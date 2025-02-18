import "regent"
require("barnes_hut_io")
require("quad_tree")

local assert = regentlib.assert
local c = regentlib.c
local pow = regentlib.pow(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")

local hdf5 = terralib.includec(os.getenv("HDF_HEADER") or "hdf5.h")
hdf5.H5F_ACC_TRUNC = 2
hdf5.H5T_STD_I32LE = hdf5.H5T_STD_I32LE_g
hdf5.H5T_STD_I64LE = hdf5.H5T_STD_I64LE_g
hdf5.H5T_IEEE_F64LE = hdf5.H5T_IEEE_F64LE_g
hdf5.H5P_DEFAULT = 0

local gee = 100.0
local delta = 0.1
local theta = 0.5
local sector_precision = 8

local launcher = require("std/launcher")
local cbarnes_hut = launcher.build_library("barnes_hut")

task update_boundaries_mass_x_max(bodies : region(body))
  where
  reads(bodies.mass_x)
do
  var mass = -math.huge
  for body in bodies do
    mass max= body.mass_x
  end
  return mass
end

task update_boundaries_mass_x_min(bodies : region(body))
  where
  reads(bodies.mass_x)
do
  var mass = math.huge
  for body in bodies do
    mass min= body.mass_x
  end
  return mass
end

task update_boundaries_mass_y_max(bodies : region(body))
  where
  reads(bodies.mass_y)
do
  var mass = -math.huge
  for body in bodies do
    mass max= body.mass_y
  end
  return mass
end

task update_boundaries_mass_y_min(bodies : region(body))
  where
  reads(bodies.mass_y)
do
  var mass = math.huge
  for body in bodies do
    mass min= body.mass_y
  end
  return mass
end

task assign_sectors(bodies : region(body), min_x : double, min_y : double, max_x : double, max_y : double)
  where
    reads(bodies.{mass_x, mass_y}),
    writes(bodies.sector)
do
  var size_x = max_x - min_x
  var size_y = max_y - min_y
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

task update_body_force_root(bodies : region(body), roots : region(ispace(int1d), quad), sector : int1d)
where
  reads (bodies.{mass_x, mass_y, mass}),
  reduces +(bodies.{force_x, force_y}),
  reads(roots.{mass_x, mass_y, mass})
do
  for i=0,sector_precision*sector_precision do
    if i ~= [int](sector) and i ~= [int](sector) - 1 and i ~= [int](sector) + 1 and i ~= [int](sector) - sector_precision and i ~= [int](sector) + sector_precision and i ~= [int](sector) - sector_precision - 1 and i ~= [int](sector) - sector_precision + 1 and i ~= [int](sector) + sector_precision - 1 and i ~= [int](sector) + sector_precision + 1 then
      for body in bodies do
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
          body.force_x += d_force * dist_x
          body.force_y += d_force * dist_y
        end
      end
    end
  end
end

task update_body_force(bodies : region(body), quads : region(ispace(int1d), quad), quad_range : region(ispace(int1d), rect1d), sector: int1d)
where
  reads (bodies.{mass_x, mass_y, mass, index}),
  reduces +(bodies.{force_x, force_y}),
  reads(quads),
  reads(quad_range)
do
  var sector_quad_range = quad_range[sector]
  var root_index = sector_quad_range.lo
  var traverse_list : int1d[1024]
  for body in bodies do
    traverse_list[0] = root_index
    var traverse_index = 0
    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1
      if quads[cur_index].type == 2 then
        var dist_x = quads[cur_index].mass_x - body.mass_x
        var dist_y = quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
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
            var d_force = gee * body.mass * quads[cur_index].mass / (dist * dist * dist)
            body.force_x += d_force * dist_x
            body.force_y += d_force * dist_y
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
              body.force_x += d_force * dist_x
              body.force_y += d_force * dist_y
            end
          end
          cur_index = quads[cur_index].next_in_leaf
        end
      end
    end
  end
end

task update_body_mass(bodies : region(body))
where
  reads(bodies.{speed_x, speed_y}),
  reduces +(bodies.{mass_x, mass_y})
do
  for body in bodies do
    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
  end
end

task update_body_speed(bodies : region(body))
where
  reads(bodies.{mass, force_x, force_y}),
  reduces +(bodies.{speed_x, speed_y})
do
  for body in bodies do
    body.speed_x += (body.force_x / body.mass * delta)
    body.speed_y += (body.force_y / body.mass * delta)
  end
end

task calculate_quad_ranges(bodies : region(body), bodies_by_sector : partition(disjoint, bodies, ispace(int1d)), quad_ranges : region(ispace(int1d), rect1d), num_quads : uint)
where
  reads(bodies),
  writes(quad_ranges)
do
    var offset = 0
    for i=0,sector_precision*sector_precision do
      var current = bodies_by_sector[i]
      var quad_size_estimate = current.ispace.volume * 2
      quad_ranges[i] = rect1d({offset, offset + quad_size_estimate})
      offset += quad_size_estimate + 1
    end
end

task run_iteration(bodies : region(body), roots : region(ispace(int1d), quad), quads : region(ispace(int1d), quad), quad_ranges : region(ispace(int1d), rect1d), num_quads : uint, conf : Config, iteration : int)
where
  roots * quads,
  reads writes(bodies),
  reads writes(roots),
  reads writes(quads),
  reads writes(quad_ranges)
do
  var sector_space = ispace(int1d, sector_precision * sector_precision)
  var min_x = math.huge
  var min_y = math.huge
  var max_x = -math.huge
  var max_y = -math.huge
  
  do
  var body_partition_index = ispace(ptr, conf.parallelism * 2)
  var bodies_partition = partition(equal, bodies, body_partition_index)
  
  __demand(__index_launch)
  for i in body_partition_index do
    min_x min= update_boundaries_mass_x_min(bodies_partition[i])
  end
  
  __demand(__index_launch)
  for i in body_partition_index do
    max_x max= update_boundaries_mass_x_max(bodies_partition[i])
  end
  
  __demand(__index_launch)
  for i in body_partition_index do
    min_y min= update_boundaries_mass_y_min(bodies_partition[i])
  end
  
  __demand(__index_launch)
  for i in body_partition_index do
    max_y max= update_boundaries_mass_y_max(bodies_partition[i])
  end

  __demand(__index_launch)
  for i in body_partition_index do
    assign_sectors(bodies_partition[i], min_x, min_y, max_x, max_y)
  end
  end -- body_partition
      
  var bodies_by_sector = partition(bodies.sector, sector_space)
  fill(bodies.{force_x, force_y}, 0)
  calculate_quad_ranges(bodies, bodies_by_sector, quad_ranges, num_quads)
  fill(quads.{nw, sw, ne, se, next_in_leaf}, -1)
  fill(quads.{mass_x, mass_y, mass, total, type}, 0)

  var quad_range_by_sector = partition(equal, quad_ranges, sector_space)
  var roots_by_sector = partition(equal, roots, sector_space)
  var quads_by_sector = image(disjoint, quads, quad_range_by_sector, quad_ranges)
  var quads_by_sector_colors = quads_by_sector.colors

  __demand(__index_launch)
  for i in sector_space do
    build_quad(bodies_by_sector[i], roots_by_sector[i], quads_by_sector[i], quad_range_by_sector[i], min_x, min_y, max_x, max_y, sector_precision, conf.leaf_size, conf.max_depth, i)
  end

  __demand(__index_launch)
  for i in sector_space do
    update_body_force_root(bodies_by_sector[i], roots, i)
  end
  
  __demand(__index_launch)
  for i in sector_space do
    update_body_force(bodies_by_sector[i], quads_by_sector[i], quad_range_by_sector[i], i)
  end

  __demand(__index_launch)
  for i=0,sector_precision*sector_precision-1 do
    update_body_force(bodies_by_sector[i], quads_by_sector[i + 1], quad_range_by_sector[i + 1], i + 1)
  end
  
  __demand(__index_launch)
  for i=1,sector_precision*sector_precision do
    update_body_force(bodies_by_sector[i], quads_by_sector[i - 1], quad_range_by_sector[i - 1], i - 1)
  end

  __demand(__index_launch)
  for i=0,sector_precision*(sector_precision-1) do
    update_body_force(bodies_by_sector[i], quads_by_sector[i + sector_precision], quad_range_by_sector[i + sector_precision], i + sector_precision)
  end

  __demand(__index_launch)
  for i=sector_precision,sector_precision*sector_precision do
    update_body_force(bodies_by_sector[i], quads_by_sector[i - sector_precision], quad_range_by_sector[i - sector_precision], i - sector_precision)
  end

  __demand(__index_launch)
  for i=0,sector_precision*(sector_precision-1)-1 do
    update_body_force(bodies_by_sector[i], quads_by_sector[i + sector_precision + 1], quad_range_by_sector[i + sector_precision + 1], i + sector_precision + 1)
  end

  __demand(__index_launch)
  for i=0,sector_precision*(sector_precision-1) do
    update_body_force(bodies_by_sector[i], quads_by_sector[i + sector_precision - 1], quad_range_by_sector[i + sector_precision - 1], i + sector_precision - 1)
  end

  __demand(__index_launch)
  for i=sector_precision+1,sector_precision*sector_precision do
    update_body_force(bodies_by_sector[i], quads_by_sector[i - sector_precision - 1], quad_range_by_sector[i - sector_precision - 1], i - sector_precision - 1)
  end

  __demand(__index_launch)
  for i=sector_precision,sector_precision*sector_precision do
    update_body_force(bodies_by_sector[i], quads_by_sector[i - sector_precision + 1], quad_range_by_sector[i - sector_precision + 1], i - sector_precision + 1)
  end

  __demand(__index_launch)
  for i in sector_space do
    update_body_mass(bodies_by_sector[i])
  end

  __demand(__index_launch)
  for i in sector_space do
    update_body_speed(bodies_by_sector[i])
  end
end

task main()
  var ts_start = c.legion_get_current_time_in_micros()
  var conf = parse_input_args()

  var all_bodies = region(ispace(ptr, conf.num_bodies), body)

  do
    var input = region(ispace(ptr, conf.num_bodies), body)
    attach(hdf5, input.{mass, mass_x, mass_y, speed_x, speed_y, index}, conf.input_file, regentlib.file_read_only)
    acquire(input.{mass, mass_x, mass_y, speed_x, speed_y, index})
    copy(input.{mass, mass_x, mass_y, speed_x, speed_y, index}, all_bodies.{mass, mass_x, mass_y, speed_x, speed_y, index})
    release(input.{mass, mass_x, mass_y, speed_x, speed_y, index})
    detach(hdf5, input.{mass, mass_x, mass_y, speed_x, speed_y, index})
  end

  var sector_space = ispace(int1d, sector_precision * sector_precision)
  var roots = region(sector_space, quad)
  var quad_ranges = region(sector_space, rect1d)

  var num_quads = conf.num_bodies * 2 + sector_precision*sector_precision

  var quads = region(ispace(int1d, num_quads), quad)

  for t=0,conf.time_steps do
    run_iteration(all_bodies, roots, quads, quad_ranges, num_quads, conf, t)
  end

  __fence(__execution, __block)
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("%d\n", (ts_end - ts_start) / 1000)
end

launcher.launch(main, "barnes_hut", cbarnes_hut.register_mappers, {"-lbarnes_hut", "-lm"})
