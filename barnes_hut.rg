import "regent"
require("barnes_hut_io")
require("quad_tree")

local assert = regentlib.assert
local c = regentlib.c
local pow = regentlib.pow(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")

local QuadTreeSizer = require("quad_tree_sizer")

local gee = 100.0
local delta = 0.1
local theta = 0.5
local elimination_threshold = 8.0
local elimination_quantity = 4

task update_boundaries(bodies : region(body), boundaries : region(boundary))
  where
  reads(bodies.{mass_x, mass_y}),
  reads(boundaries),
  reduces min(boundaries.{min_x, min_y}),
  reduces max(boundaries.{max_x, max_y})
do
  for body in bodies do
    boundaries[0].min_x min = min(body.mass_x, boundaries[0].min_x)
    boundaries[0].min_y min = min(body.mass_y, boundaries[0].min_y)
    boundaries[0].max_x max = max(body.mass_x, boundaries[0].max_x)
    boundaries[0].max_y max = max(body.mass_y, boundaries[0].max_y)
  end
end

task assign_sectors(bodies : region(body), min_x : double, min_y : double, size : double, sector_precision : uint)
  where
  reads(bodies.{mass_x, mass_y, sector}),
  writes(bodies.sector)
do
  for body in bodies do
    var sector_x : int64 = cmath.floor((body.mass_x - min_x) / (size / sector_precision))
    if (sector_x >= sector_precision) then
      sector_x = sector_x - 1
    end

    var sector_y: int64 = cmath.floor((body.mass_y - min_y) / (size / sector_precision))
    if (sector_y >= sector_precision) then
      sector_y = sector_y - 1
    end

    body.sector = sector_x + sector_y * sector_precision
  end
end

task size_quad(bodies : region(body), quad_size : region(uint), min_x : double, min_y : double, size : double, sector_precision : uint, leaf_size : uint, sector : int1d)
  where reads(bodies.{mass_x, mass_y, index}),
  writes (quad_size)
do
  var chunk = create_quad_chunk(2048)
  var sector_x : int64 = sector % sector_precision
  var sector_y : int64 = cmath.floor(sector / sector_precision)
  var center_x = min_x + (sector_x + 0.5) * size / sector_precision
  var center_y = min_y + (sector_y + 0.5) * size / sector_precision

  var root = init_placeholder(chunk)
  root.center_x = center_x
  root.center_y = center_y
  root.size = size / sector_precision
  root.type = 2

  for body in bodies do
    var body_quad = init_placeholder(chunk)
    body_quad.mass_x = body.mass_x
    body_quad.mass_y = body.mass_y
    body_quad.type = 1
    add_placeholder(root, body_quad, chunk, leaf_size)
  end

  quad_size[sector] = count(chunk, true)
end

task update_body_positions(bodies : region(body), quads : region(ispace(int1d), quad), root_index : uint)
where
  reads writes(bodies),
  reads(quads)
do
  var traverse_list : int1d[1024]
  for body in bodies do
    traverse_list[0] = root_index
    var traverse_index = 0
    var force_x : double = 0.0f
    var force_y : double = 0.0f

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1
      -- c.printf("traverse index %d quad %d\n", traverse_index + 1, cur_index)

      if quads[cur_index].type == 2 then
        var dist = sqrt((body.mass_x - quads[cur_index].mass_x) * (body.mass_x - quads[cur_index].mass_x) + (body.mass_y - quads[cur_index].mass_y) * (body.mass_y - quads[cur_index].mass_y))
        if dist > 1.0 then
          if quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = quads[cur_index].sw
            end

            if quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = quads[cur_index].nw
            end

            if quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = quads[cur_index].se
            end

            if quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * quads[cur_index].mass / (dist * dist)
            var xn = (quads[cur_index].mass_x - body.mass_x) / dist
            var yn = (quads[cur_index].mass_y - body.mass_y) / dist
            -- c.printf("Heuristic updating body %d: d_force_x %f d_force_y %f\n", body.index, d_force_x, d_force_x)

            force_x += d_force * xn
            force_y += d_force * yn
          end
        end
      else
        while cur_index ~= -1 do
          if quads[cur_index].index ~= body.index then
            var dist = sqrt((body.mass_x - quads[cur_index].mass_x) * (body.mass_x - quads[cur_index].mass_x) + (body.mass_y - quads[cur_index].mass_y) * (body.mass_y - quads[cur_index].mass_y))
            if dist > 1.0 then
              var d_force = gee * body.mass * quads[cur_index].mass / (dist * dist)
              var xn = (quads[cur_index].mass_x - body.mass_x) / dist
              var yn = (quads[cur_index].mass_y - body.mass_y) / dist
              force_x += d_force * xn
              force_y += d_force * yn
              -- c.printf("Updating body %d: d_force_x %f d_force_y %f d_force %f body_mass %f cur_mass %f dist %f xn %f yn %f\n", body.index, d_force_x, d_force_x, d_force, body.mass, quads[cur_index].mass, dist, xn, yn)
              end
            end
          cur_index = quads[cur_index].next_in_leaf
        end
      end
    end

    -- c.printf("Updating body %d: force_x %f force_y %f\n", body.index, force_x, force_y)

    body.mass_x = body.mass_x + body.speed_x * delta
    body.mass_y = body.mass_y + body.speed_y * delta
    body.speed_x = body.speed_x + force_x / body.mass * delta
    body.speed_y = body.speed_y + force_y / body.mass * delta
  end
end

task eliminate_outliers(bodies : region(body), quad_size : region(uint), root_mass_x : double, root_mass_y : double, root_mass : double, sector_size : double, sector : int1d)
where
  reads(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(quad_size),
  writes(bodies.eliminated)
do
  if quad_size[sector] < elimination_quantity then
    for body in bodies do
      var dx = root_mass_x - body.mass_x
      var dy = root_mass_y - body.mass_y
      var d = sqrt(dx * dx + dy * dy)
      if d > elimination_threshold * sector_size then
        var nx = dx / d
        var ny = dy / d
        var relative_speed = body.speed_x * nx + body.speed_y * ny
        if relative_speed < 0 then
          var escape_speed = sqrt(2 * gee * root_mass / d)
          if relative_speed < -2 * escape_speed then
            body.eliminated = 1
          end
        end
      end
    end
  end
end

task main()
  var ts_start = c.legion_get_current_time_in_micros()
  var conf = parse_input_args()
  
  var num_bodies = get_number_of_bodies(conf)
  c.printf("Loading %d bodies\n", num_bodies)
  var all_bodies = region(ispace(ptr, num_bodies), body)

  load_bodies(all_bodies, conf, num_bodies)

  if conf.csv_dir_set then
    print_bodies_csv_initial(all_bodies, conf)
  end

  var boundaries = region(ispace(ptr, 1), boundary)
  
  if conf.svg_dir_set then
    boundaries[0] = { min_x = all_bodies[0].mass_x, min_y = all_bodies[0].mass_y, max_x = all_bodies[0].mass_x, max_y = all_bodies[0].mass_y }
    update_boundaries(all_bodies, boundaries)
    print_bodies_svg(all_bodies, boundaries, conf, 0)
  end

  var sector_precision : uint = pow(2, conf.N)

  var quad_range_space = ispace(int1d, sector_precision * sector_precision + 1)
  var quad_ranges = region(quad_range_space, rect1d)
  var quad_sizes = region(ispace(ptr, sector_precision * sector_precision), uint)
  var sector_index = ispace(int1d, sector_precision * sector_precision)
  var sector_quad_sizes = partition(equal, quad_sizes, sector_index)

  var num_quads = min(num_bodies * 2, 2400000)
  for i=0,conf.N do
    num_quads += pow(4, i)
  end
  num_quads += sector_precision*sector_precision

  var quads = region(ispace(int1d, num_quads), quad)

  for t=0,conf.time_steps do
      var iter_start = c.legion_get_current_time_in_micros()
      
      var elimination_partition = partition(all_bodies.eliminated, ispace(int1d, 2))
      var bodies = elimination_partition[0]

      boundaries[0] = { min_x = bodies[0].mass_x, min_y = bodies[0].mass_y, max_x = bodies[0].mass_x, max_y = bodies[0].mass_y }
      
      var body_partition_index = ispace(ptr, conf.parallelism * 2)

      var bodies_partition = partition(equal, bodies, body_partition_index)
      for i in body_partition_index do
        update_boundaries(bodies_partition[i], boundaries)
      end

      var min_x = boundaries[0].min_x
      var min_y = boundaries[0].min_y
      var size_x = boundaries[0].max_x - min_x
      var size_y = boundaries[0].max_y - min_y
      var size = max(size_x, size_y)

      __demand(__parallel)
      for i in body_partition_index do
        assign_sectors(bodies_partition[i], min_x, min_y, size, sector_precision)
      end
      
      var bodies_by_sector = partition(bodies.sector, sector_index)

      if conf.fixed_partition_size == -1 then
        __demand(__parallel)
        for i in sector_index do
          size_quad(bodies_by_sector[i], sector_quad_sizes[i], min_x, min_y, size, sector_precision, conf.leaf_size, i)
        end
      else
        for i in sector_index do
          quad_sizes[i] = conf.fixed_partition_size
        end
      end

      var offset = 0
      for i=0,sector_precision*sector_precision do
        quad_ranges[i] = rect1d({offset, offset + quad_sizes[i] - 1})
        offset += quad_sizes[i]
      end

      quad_ranges[sector_precision * sector_precision] = rect1d({offset, num_quads - 1})

      fill(quads.{nw, sw, ne, se, next_in_leaf}, -1)
      fill(quads.{mass_x, mass_y, mass, total, type}, 0)

      var quad_range_by_sector = partition(equal, quad_ranges, quad_range_space) 
      var quads_by_sector = image(quads, quad_range_by_sector, quad_ranges)

      var quads_by_sector_colors = quads_by_sector.colors
      var quads_by_sector_disjoint = dynamic_cast(partition(disjoint, quads, quads_by_sector_colors), quads_by_sector)

      __demand(__parallel)
      for i in sector_index do
        build_quad(bodies_by_sector[i], quads_by_sector_disjoint[i], quad_ranges, min_x, min_y, size, sector_precision, conf.leaf_size, i)
      end

      var to_merge : int[64][64]
      for i=0,sector_precision do
        for j=0,sector_precision do
          var quad_range = quad_ranges[i + j*sector_precision]
          var root_index = quad_range.lo
          if quads[root_index].total > 0 then
            to_merge[i][j] = root_index
          else
            to_merge[i][j] = -1
          end
        end
      end

      var allocation_index = num_quads - 1
      var level = sector_precision
      while level > 1 do
        var next_level = level / 2
        for i=0,next_level do
          for j=0,next_level do
            quads[allocation_index].size = size / next_level
            quads[allocation_index].center_x = min_x + size / next_level * (i + 0.5)
            quads[allocation_index].center_y = min_y + size / next_level * (j + 0.5)
            quads[allocation_index].type = 2

            quads[allocation_index].mass = 0
            quads[allocation_index].mass_x = 0
            quads[allocation_index].mass_y = 0
            quads[allocation_index].total = 0

            if to_merge[2*i][2*j] ~= -1 then
              quads[allocation_index].sw = to_merge[2*i][2*j]
              quads[allocation_index].mass += quads[to_merge[2*i][2*j]].mass
              quads[allocation_index].mass_x += quads[to_merge[2*i][2*j]].mass_x * quads[to_merge[2*i][2*j]].mass
              quads[allocation_index].mass_y += quads[to_merge[2*i][2*j]].mass_y * quads[to_merge[2*i][2*j]].mass
              quads[allocation_index].total += quads[to_merge[2*i][2*j]].total
            end

            if to_merge[2*i][2*j+1] ~= -1 then
              quads[allocation_index].nw = to_merge[2*i][2*j+1]
              quads[allocation_index].mass += quads[to_merge[2*i][2*j+1]].mass
              quads[allocation_index].mass_x += quads[to_merge[2*i][2*j+1]].mass_x * quads[to_merge[2*i][2*j+1]].mass
              quads[allocation_index].mass_y += quads[to_merge[2*i][2*j+1]].mass_y * quads[to_merge[2*i][2*j+1]].mass
              quads[allocation_index].total += quads[to_merge[2*i][2*j+1]].total
            end

            if to_merge[2*i+1][2*j] ~= -1 then
              quads[allocation_index].se = to_merge[2*i+1][2*j]
              quads[allocation_index].mass += quads[to_merge[2*i+1][2*j]].mass
              quads[allocation_index].mass_x += quads[to_merge[2*i+1][2*j]].mass_x * quads[to_merge[2*i+1][2*j]].mass
              quads[allocation_index].mass_y += quads[to_merge[2*i+1][2*j]].mass_y * quads[to_merge[2*i+1][2*j]].mass
              quads[allocation_index].total += quads[to_merge[2*i+1][2*j]].total
            end

            if to_merge[2*i+1][2*j+1] ~= -1 then
              quads[allocation_index].ne = to_merge[2*i+1][2*j+1]
              quads[allocation_index].mass += quads[to_merge[2*i+1][2*j+1]].mass
              quads[allocation_index].mass_x += quads[to_merge[2*i+1][2*j+1]].mass_x * quads[to_merge[2*i+1][2*j+1]].mass
              quads[allocation_index].mass_y += quads[to_merge[2*i+1][2*j+1]].mass_y * quads[to_merge[2*i+1][2*j+1]].mass
              quads[allocation_index].total += quads[to_merge[2*i+1][2*j+1]].total
            end

            if quads[allocation_index].total > 0 then
              quads[allocation_index].mass_x = quads[allocation_index].mass_x / quads[allocation_index].mass
              quads[allocation_index].mass_y = quads[allocation_index].mass_y / quads[allocation_index].mass
              to_merge[i][j] = allocation_index
            else
              to_merge[i][j] = -1
            end

            allocation_index = allocation_index - 1
          end
        end
        level = next_level
      end

       -- for i in quads_index do
        -- if quads[i].total ~= 0 then
          -- c.printf("%d Quad index: %d, type %d mass_x %f, mass_y %f, mass %f, center_x %f, center_y %f, size %f, total %d, sw %d, nw %d, se %d, ne %d\n", i, quads[i].index, quads[i].type, quads[i].mass_x, quads[i].mass_y, quads[i].mass, quads[i].center_x, quads[i].center_y, quads[i].size, quads[i].total, quads[i].sw, quads[i].nw, quads[i].se, quads[i].ne)
        -- end
       -- end

      -- var i = allocation_index + 1
      -- c.printf("\n%d Root index: %d, type %d mass_x %f, mass_y %f, mass %f, center_x %f, center_y %f, size %f, total %d, sw %d, nw %d, se %d, ne %d\n", i, quads[i].index, quads[i].type, quads[i].mass_x, quads[i].mass_y, quads[i].mass, quads[i].center_x, quads[i].center_y, quads[i].size, quads[i].total, quads[i].sw, quads[i].nw, quads[i].se, quads[i].ne)

      var root_index = allocation_index + 1
      __demand(__parallel)
      for i in body_partition_index do
        update_body_positions(bodies_partition[i], quads, root_index)
      end

      var root = quads[root_index]
      var root_mass_x = root.mass_x
      var root_mass_y = root.mass_y
      var root_mass = root.mass

      __demand(__parallel)
      for x=0,sector_precision do
        eliminate_outliers(bodies_by_sector[x], sector_quad_sizes[x], root_mass_x, root_mass_y, root_mass, size, x)
      end

      var start_index = sector_precision * (sector_precision - 1)
      var end_index = sector_precision * sector_precision - 1

      __demand(__parallel)
      for x=start_index,end_index do
        eliminate_outliers(bodies_by_sector[x], sector_quad_sizes[x], root_mass_x, root_mass_y, root_mass, size, x)
      end

      start_index = sector_precision
      end_index = sector_precision * (sector_precision - 1)

      -- __demand(__parallel)
      for y=start_index,end_index,sector_precision do
        eliminate_outliers(bodies_by_sector[y], sector_quad_sizes[y], root_mass_x, root_mass_y, root_mass, size, y)
      end

      start_index = sector_precision + sector_precision - 1
      end_index = sector_precision * (sector_precision - 1) + sector_precision - 1

      -- __demand(__parallel)
      for y=start_index,end_index,sector_precision do
        eliminate_outliers(bodies_by_sector[y], sector_quad_sizes[y], root_mass_x, root_mass_y, root_mass, size, y)
      end

      var iter_end = c.legion_get_current_time_in_micros()
      c.printf("Iteration time: %d ms\n", (iter_end - iter_start) / 1000)

      if conf.csv_dir_set then
        print_bodies_csv_update(bodies, conf, t+1)
      end

      if conf.svg_dir_set then
        print_bodies_svg(bodies, boundaries, conf, t+1)
      end
  end
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("Total time: %d ms\n", (ts_end - ts_start) / 1000)
end
regentlib.start(main)
