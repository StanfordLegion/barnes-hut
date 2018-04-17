import "regent"
require("barnes_hut_io")
require("quad_tree")

local assert = regentlib.assert
local c = regentlib.c
local pow = regentlib.pow(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")

local gee = 100.0
local delta = 0.1
local theta = 0.5
local elimination_threshold = 8.0
local elimination_quantity = 4

task init_boundaries(bodies : region(body), boundaries : region(boundary))
  where
  reads(bodies.{mass_x, mass_y}),
  writes(boundaries)
do
  boundaries[0] = { min_x = bodies[0].mass_x, min_y = bodies[0].mass_y, max_x = bodies[0].mass_x, max_y = bodies[0].mass_y }
end

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

task assign_sectors(bodies : region(body), boundaries : region(boundary), sector_precision : uint)
  where
    reads(bodies.{mass_x, mass_y}),
    writes(bodies.sector),
    reads(boundaries)
do
  var min_x = boundaries[0].min_x
  var min_y = boundaries[0].min_y
  var size_x = boundaries[0].max_x - min_x
  var size_y = boundaries[0].max_y - min_y
  var size = max(size_x, size_y)
  var sector_size = size / sector_precision

  for body in bodies do
    var sector_x : int64 = cmath.floor((body.mass_x - min_x) / sector_size)
    if (sector_x >= sector_precision) then
      sector_x = sector_x - 1
    end

    var sector_y: int64 = cmath.floor((body.mass_y - min_y) / sector_size)
    if (sector_y >= sector_precision) then
      sector_y = sector_y - 1
    end

    var sector = sector_x + sector_y * sector_precision
    body.sector = sector
  end
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
        var dist_x = quads[cur_index].mass_x - body.mass_x
        var dist_y = quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if quads[cur_index].size / dist >= theta then
            -- assert(traverse_index < 1020, "possible traverse list overflow")
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
            var d_force = gee * body.mass * quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if quads[cur_index].index ~= body.index then
            var dist_x = quads[cur_index].mass_x - body.mass_x
            var dist_y = quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
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

task eliminate_outliers(bodies : region(body), quads : region(ispace(int1d), quad), root_index : uint, sector_precision : double)
where
  reads(bodies.{mass_x, mass_y, speed_x, speed_y}),
  writes(bodies.eliminated),
  reads(quads.{mass_x, mass_y, mass, size})
do
  var root_mass_x = quads[root_index].mass_x
  var root_mass_y = quads[root_index].mass_y
  var root_mass = quads[root_index].mass
  var sector_size = quads[root_index].size / sector_precision

  if bodies.ispace.volume < elimination_quantity then
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

task merge_tree(quads : region(ispace(int1d), quad), quad_ranges : region(ispace(int1d), rect1d), boundaries : region(boundary), sector_precision : uint)
where
  reads(quads),
  reads(quad_ranges),
  writes(quads),
  reads(boundaries)
do
  var min_x = boundaries[0].min_x
  var min_y = boundaries[0].min_y
  var size_x = boundaries[0].max_x - min_x
  var size_y = boundaries[0].max_y - min_y
  var size = max(size_x, size_y)

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

  var quad_range = quad_ranges[sector_precision * sector_precision]
  var allocation_index = quad_range.hi
  var level = sector_precision
  while level > 1 do
    var next_level = level / 2
    for i=0,next_level do
      for j=0,next_level do
        var quad_size = size / next_level
        quads[allocation_index].size = quad_size
        quads[allocation_index].center_x = min_x + quad_size * (i + 0.5)
        quads[allocation_index].center_y = min_y + quad_size * (j + 0.5)
        quads[allocation_index].type = 2

        var mass = 0.0f
        var mass_x = 0.0f
        var mass_y = 0.0f
        var total = 0

        if to_merge[2*i][2*j] ~= -1 then
          quads[allocation_index].sw = to_merge[2*i][2*j]
          var to_merge_mass = quads[to_merge[2*i][2*j]].mass
          mass += to_merge_mass
          mass_x += quads[to_merge[2*i][2*j]].mass_x * to_merge_mass
          mass_y += quads[to_merge[2*i][2*j]].mass_y * to_merge_mass
          total += quads[to_merge[2*i][2*j]].total
        end

        if to_merge[2*i][2*j+1] ~= -1 then
          quads[allocation_index].nw = to_merge[2*i][2*j+1]
          var to_merge_mass = quads[to_merge[2*i][2*j+1]].mass
          mass += to_merge_mass
          mass_x += quads[to_merge[2*i][2*j+1]].mass_x * to_merge_mass
          mass_y += quads[to_merge[2*i][2*j+1]].mass_y * to_merge_mass
          total += quads[to_merge[2*i][2*j+1]].total
        end

        if to_merge[2*i+1][2*j] ~= -1 then
          quads[allocation_index].se = to_merge[2*i+1][2*j]
          var to_merge_mass = quads[to_merge[2*i+1][2*j]].mass
          mass += to_merge_mass
          mass_x += quads[to_merge[2*i+1][2*j]].mass_x * to_merge_mass
          mass_y += quads[to_merge[2*i+1][2*j]].mass_y * to_merge_mass
          total += quads[to_merge[2*i+1][2*j]].total
        end

        if to_merge[2*i+1][2*j+1] ~= -1 then
          quads[allocation_index].ne = to_merge[2*i+1][2*j+1]
          var to_merge_mass = quads[to_merge[2*i+1][2*j+1]].mass
          mass += to_merge_mass
          mass_x += quads[to_merge[2*i+1][2*j+1]].mass_x * to_merge_mass
          mass_y += quads[to_merge[2*i+1][2*j+1]].mass_y * to_merge_mass
          total += quads[to_merge[2*i+1][2*j+1]].total
        end

        if total > 0 then
          quads[allocation_index].mass = mass
          quads[allocation_index].mass_x = mass_x / mass
          quads[allocation_index].mass_y = mass_y / mass
          quads[allocation_index].total = total
          to_merge[i][j] = allocation_index
        else
          to_merge[i][j] = -1
        end

        allocation_index = allocation_index - 1
      end
    end
    level = next_level
  end  
end

task calculate_quad_ranges(bodies : region(body), bodies_by_sector : partition(disjoint, bodies, ispace(int1d)), quad_ranges : region(ispace(int1d), rect1d), num_quads : uint, sector_precision : uint)
where
  reads (bodies),
  writes (quad_ranges)
do
    var offset = 0
    for i=0,sector_precision*sector_precision do
      var current = bodies_by_sector[i]
      var quad_size_estimate = current.ispace.volume * 12 / 5
      quad_ranges[i] = rect1d({offset, offset + quad_size_estimate})
      offset += quad_size_estimate + 1
    end

    quad_ranges[sector_precision * sector_precision] = rect1d({offset, num_quads - 1})
end

task run_iteration(all_bodies : region(body), quads : region(ispace(int1d), quad), quad_ranges : region(ispace(int1d), rect1d), boundaries : region(boundary), num_quads : uint, sector_precision : uint, merge_quads : uint, conf : Config)
where
  reads writes(all_bodies),
  reads writes(quads),
  reads writes(quad_ranges),
  reads writes(boundaries)
do
  var quad_range_space = ispace(int1d, sector_precision * sector_precision + 1)
  var sector_index = ispace(int1d, sector_precision * sector_precision)
  var elimination_partition = partition(all_bodies.eliminated, ispace(int1d, 2))
  var bodies = elimination_partition[0]

  init_boundaries(bodies, boundaries)      
  var body_partition_index = ispace(ptr, conf.parallelism * 2)

  var bodies_partition = partition(equal, bodies, body_partition_index)
  for i in body_partition_index do
    update_boundaries(bodies_partition[i], boundaries)
  end

  __demand(__parallel)
  for i in body_partition_index do
    assign_sectors(bodies_partition[i], boundaries, sector_precision)
  end
      
  var bodies_by_sector = partition(bodies.sector, sector_index)

  calculate_quad_ranges(bodies, bodies_by_sector, quad_ranges, num_quads, sector_precision)

  fill(quads.{nw, sw, ne, se, next_in_leaf}, -1)
  fill(quads.{mass_x, mass_y, mass, total, type}, 0)

  var quad_range_by_sector = partition(equal, quad_ranges, quad_range_space)
  var quads_by_sector = image(quads, quad_range_by_sector, quad_ranges)

  var quads_by_sector_colors = quads_by_sector.colors
  var quads_by_sector_disjoint = dynamic_cast(partition(disjoint, quads, quads_by_sector_colors), quads_by_sector)

  __demand(__parallel)
  for i in sector_index do
    var current = bodies_by_sector[i]
    build_quad(bodies_by_sector[i], quads_by_sector_disjoint[i], quad_ranges, boundaries, sector_precision, conf.leaf_size, conf.max_depth, i)
  end

  merge_tree(quads, quad_ranges, boundaries, sector_precision)

  var root_index = num_quads - merge_quads
  __demand(__parallel)
  for i in body_partition_index do
    update_body_positions(bodies_partition[i], quads, root_index)
  end

  __demand(__parallel)
  for x=0,sector_precision do
    eliminate_outliers(bodies_by_sector[x], quads, root_index, sector_precision)
  end

  var start_index = sector_precision * (sector_precision - 1)
  var end_index = sector_precision * sector_precision - 1

  __demand(__parallel)
  for x=start_index,end_index do
    eliminate_outliers(bodies_by_sector[x], quads, root_index, sector_precision)
  end

  start_index = sector_precision
  end_index = sector_precision * (sector_precision - 1)

  -- __demand(__parallel)
  for y=start_index,end_index,sector_precision do
    eliminate_outliers(bodies_by_sector[y], quads, root_index, sector_precision)
  end

  start_index = sector_precision + sector_precision - 1
  end_index = sector_precision * (sector_precision - 1) + sector_precision - 1

  -- __demand(__parallel)
  for y=start_index,end_index,sector_precision do
    eliminate_outliers(bodies_by_sector[y], quads, root_index, sector_precision)
  end

  __delete(elimination_partition)
  __delete(bodies_partition)
  __delete(bodies_by_sector)
  __delete(quad_range_by_sector)
  __delete(quads_by_sector)
  __delete(quads_by_sector_disjoint)
end

task main()
  var ts_start = c.legion_get_current_time_in_micros()
  var conf = parse_input_args()
  
  var num_bodies = get_number_of_bodies(conf)
  -- c.printf("Loading %d bodies\n", num_bodies)
  var all_bodies = region(ispace(ptr, num_bodies), body)

  load_bodies(all_bodies, conf, num_bodies)

  -- if conf.csv_dir_set then
    -- print_bodies_csv_initial(all_bodies, conf)
  -- end

  var boundaries = region(ispace(ptr, 1), boundary)
  
--  if conf.svg_dir_set then
--    boundaries[0] = { min_x = all_bodies[0].mass_x, min_y = all_bodies[0].mass_y, max_x = all_bodies[0].mass_x, max_y = all_bodies[0].mass_y }
--    update_boundaries(all_bodies, boundaries)
--    print_bodies_svg(all_bodies, boundaries, conf, 0)
--  end

  var sector_precision : uint = pow(2, conf.N)

  var quad_range_space = ispace(int1d, sector_precision * sector_precision + 1)
  var quad_ranges = region(quad_range_space, rect1d)

  var merge_quads = 0
  for i=0,conf.N do
    merge_quads += pow(4, i)
  end
  var num_quads = num_bodies * 12 / 5 + merge_quads + sector_precision*sector_precision

  var quads = region(ispace(int1d, num_quads), quad)

  for t=0,conf.time_steps do
    run_iteration(all_bodies, quads, quad_ranges, boundaries, num_quads, sector_precision, merge_quads, conf)
    -- if conf.csv_dir_set then
      -- print_bodies_csv_update(bodies, conf, t+1)
    -- end

    -- if conf.svg_dir_set then
      -- print_bodies_svg(bodies, boundaries, conf, t+1)
    -- end
  end

  __fence(__execution, __block)
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("%d\n", (ts_end - ts_start) / 1000)
end

if os.getenv('SAVEOBJ') == '1' then
  local exe = os.getenv('OBJNAME') or "barnes_hut"
  regentlib.saveobj(main, exe, "executable")
else
  regentlib.start(main)
end
