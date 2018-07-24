import "regent"
require("barnes_hut_io")
require("quad_tree")

local assert = regentlib.assert
local c = regentlib.c
local pow = regentlib.pow(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")

local hdf5 = terralib.includec("hdf5.h", {"-I", "/share/software/user/open/hdf5/1.10.2/include", "-I", "/share/software/user/open/openmpi/2.0.2/include"})
hdf5.H5F_ACC_TRUNC = 2
hdf5.H5T_STD_I32LE = hdf5.H5T_STD_I32LE_g
hdf5.H5T_STD_I64LE = hdf5.H5T_STD_I64LE_g
hdf5.H5T_IEEE_F64LE = hdf5.H5T_IEEE_F64LE_g
hdf5.H5P_DEFAULT = 0

local gee = 100.0
local delta = 0.1
local theta = 0.5
local sector_precision = 8

local cbarnes_hut
do
  assert(os.getenv('LG_RT_DIR') ~= nil, "$LG_RT_DIR should be set!")
  local root_dir = arg[0]:match(".*/") or "./"
  local runtime_dir = os.getenv("LG_RT_DIR") .. "/"
  local barnes_hut_cc = root_dir .. "barnes_hut.cc"
  local barnes_hut_so
  if os.getenv('SAVEOBJ') == '1' then
    barnes_hut_so = root_dir .. "libbarnes_hut.so"
  else
    barnes_hut_so = os.tmpname() .. ".so" -- root_dir .. "mapper.so"
  end
  local cxx = os.getenv('CXX') or 'c++'

  local cxx_flags = os.getenv('CC_FLAGS') or ''
  cxx_flags = cxx_flags .. " -O2 -Wall -Werror"
  if os.execute('test "$(uname)" = Darwin') == 0 then
    cxx_flags =
      (cxx_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cxx_flags = cxx_flags .. " -shared -fPIC"
  end

  local cmd = (cxx .. " " .. cxx_flags .. " -I " .. runtime_dir .. " " ..
                 barnes_hut_cc .. " -o " .. barnes_hut_so)
  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. barnes_hut_cc)
    assert(false, "")
  end
  terralib.linklibrary(barnes_hut_so)
  cbarnes_hut =
    terralib.includec("barnes_hut.h", {"-I", root_dir, "-I", runtime_dir})
end

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

task assign_sectors_first(bodies : region(body), min_x : double, min_y : double, max_x : double, max_y : double)
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

task update_bodies_first(bodies : region(body), roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  right_quads : region(ispace(int1d), quad), right_quad_range : region(ispace(int1d), rect1d),
  down_quads : region(ispace(int1d), quad), down_quad_range : region(ispace(int1d), rect1d),
  down_right_quads : region(ispace(int1d), quad), down_right_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(right_quads),
  reads(right_quad_range),
  reads(down_quads),
  reads(down_quad_range),
  reads(down_right_quads),
  reads(down_right_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= 0 and i ~= 1 and i ~= sector_precision and i ~= sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[0]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = right_quad_range[1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if right_quads[cur_index].type == 2 then
        var dist_x = right_quads[cur_index].mass_x - body.mass_x
        var dist_y = right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].sw
            end

            if right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].nw
            end

            if right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].se
            end

            if right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if right_quads[cur_index].index ~= body.index then
            var dist_x = right_quads[cur_index].mass_x - body.mass_x
            var dist_y = right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_quad_range[sector_precision]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_quads[cur_index].type == 2 then
        var dist_x = down_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].sw
            end

            if down_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].nw
            end

            if down_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].se
            end

            if down_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_quads[cur_index].index ~= body.index then
            var dist_x = down_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_right_quad_range[sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_right_quads[cur_index].type == 2 then
        var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].sw
            end

            if down_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].nw
            end

            if down_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].se
            end

            if down_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_right_quads[cur_index].index ~= body.index then
            var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_right_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_first_row(bodies : region(body), sector: int1d, roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  left_quads : region(ispace(int1d), quad), left_quad_range : region(ispace(int1d), rect1d),
  right_quads : region(ispace(int1d), quad), right_quad_range : region(ispace(int1d), rect1d),
  down_quads : region(ispace(int1d), quad), down_quad_range : region(ispace(int1d), rect1d),
  down_left_quads : region(ispace(int1d), quad), down_left_quad_range : region(ispace(int1d), rect1d),
  down_right_quads : region(ispace(int1d), quad), down_right_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(left_quads),
  reads(left_quad_range),
  reads(right_quads),
  reads(right_quad_range),
  reads(down_quads),
  reads(down_quad_range),
  reads(down_left_quads),
  reads(down_left_quad_range),
  reads(down_right_quads),
  reads(down_right_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= [int](sector) and i ~= [int](sector) - 1 and i ~= [int](sector) + 1 and i ~= [int](sector) + sector_precision and i ~= [int](sector) + sector_precision - 1 and i ~= [int](sector) + sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = left_quad_range[sector - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if left_quads[cur_index].type == 2 then
        var dist_x = left_quads[cur_index].mass_x - body.mass_x
        var dist_y = left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].sw
            end

            if left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].nw
            end

            if left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].se
            end

            if left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if left_quads[cur_index].index ~= body.index then
            var dist_x = left_quads[cur_index].mass_x - body.mass_x
            var dist_y = left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = right_quad_range[sector + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if right_quads[cur_index].type == 2 then
        var dist_x = right_quads[cur_index].mass_x - body.mass_x
        var dist_y = right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].sw
            end

            if right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].nw
            end

            if right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].se
            end

            if right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if right_quads[cur_index].index ~= body.index then
            var dist_x = right_quads[cur_index].mass_x - body.mass_x
            var dist_y = right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_quad_range[sector + sector_precision]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_quads[cur_index].type == 2 then
        var dist_x = down_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].sw
            end

            if down_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].nw
            end

            if down_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].se
            end

            if down_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_quads[cur_index].index ~= body.index then
            var dist_x = down_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_left_quad_range[sector + sector_precision - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_left_quads[cur_index].type == 2 then
        var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].sw
            end

            if down_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].nw
            end

            if down_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].se
            end

            if down_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_left_quads[cur_index].index ~= body.index then
            var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_right_quad_range[sector + sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_right_quads[cur_index].type == 2 then
        var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].sw
            end

            if down_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].nw
            end

            if down_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].se
            end

            if down_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_right_quads[cur_index].index ~= body.index then
            var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_right_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_middle_first(bodies : region(body), roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  right_quads : region(ispace(int1d), quad), right_quad_range : region(ispace(int1d), rect1d),
  up_quads : region(ispace(int1d), quad), up_quad_range : region(ispace(int1d), rect1d),
  up_right_quads : region(ispace(int1d), quad), up_right_quad_range : region(ispace(int1d), rect1d),
  down_quads : region(ispace(int1d), quad), down_quad_range : region(ispace(int1d), rect1d),
  down_right_quads : region(ispace(int1d), quad), down_right_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(right_quads),
  reads(right_quad_range),
  reads(up_quads),
  reads(up_quad_range),
  reads(up_right_quads),
  reads(up_right_quad_range),
  reads(down_quads),
  reads(down_quad_range),
  reads(down_right_quads),
  reads(down_right_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= 0 and i ~= 1 and i ~= sector_precision and i ~= sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector_precision]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = right_quad_range[sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if right_quads[cur_index].type == 2 then
        var dist_x = right_quads[cur_index].mass_x - body.mass_x
        var dist_y = right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].sw
            end

            if right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].nw
            end

            if right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].se
            end

            if right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if right_quads[cur_index].index ~= body.index then
            var dist_x = right_quads[cur_index].mass_x - body.mass_x
            var dist_y = right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_quad_range[0]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_quads[cur_index].type == 2 then
        var dist_x = up_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].sw
            end

            if up_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].nw
            end

            if up_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].se
            end

            if up_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_quads[cur_index].index ~= body.index then
            var dist_x = up_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_right_quad_range[1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_right_quads[cur_index].type == 2 then
        var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].sw
            end

            if up_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].nw
            end

            if up_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].se
            end

            if up_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_right_quads[cur_index].index ~= body.index then
            var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_quad_range[sector_precision * 2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_quads[cur_index].type == 2 then
        var dist_x = down_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].sw
            end

            if down_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].nw
            end

            if down_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].se
            end

            if down_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_quads[cur_index].index ~= body.index then
            var dist_x = down_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_right_quad_range[sector_precision * 2 + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_right_quads[cur_index].type == 2 then
        var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].sw
            end

            if down_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].nw
            end

            if down_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].se
            end

            if down_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_right_quads[cur_index].index ~= body.index then
            var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_right_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_middle(bodies : region(body), sector: int1d, roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  left_quads : region(ispace(int1d), quad), left_quad_range : region(ispace(int1d), rect1d),
  right_quads : region(ispace(int1d), quad), right_quad_range : region(ispace(int1d), rect1d),
  up_quads : region(ispace(int1d), quad), up_quad_range : region(ispace(int1d), rect1d),
  up_left_quads : region(ispace(int1d), quad), up_left_quad_range : region(ispace(int1d), rect1d),
  up_right_quads : region(ispace(int1d), quad), up_right_quad_range : region(ispace(int1d), rect1d),
  down_quads : region(ispace(int1d), quad), down_quad_range : region(ispace(int1d), rect1d),
  down_left_quads : region(ispace(int1d), quad), down_left_quad_range : region(ispace(int1d), rect1d),
  down_right_quads : region(ispace(int1d), quad), down_right_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(left_quads),
  reads(left_quad_range),
  reads(right_quads),
  reads(right_quad_range),
  reads(up_quads),
  reads(up_quad_range),
  reads(up_left_quads),
  reads(up_left_quad_range),
  reads(up_right_quads),
  reads(up_right_quad_range),
  reads(down_quads),
  reads(down_quad_range),
  reads(down_left_quads),
  reads(down_left_quad_range),
  reads(down_right_quads),
  reads(down_right_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= [int](sector) and i ~= [int](sector) - 1 and i ~= [int](sector) + 1 and i ~= [int](sector) + sector_precision and i ~= [int](sector) + sector_precision - 1 and i ~= [int](sector) + sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = left_quad_range[sector - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if left_quads[cur_index].type == 2 then
        var dist_x = left_quads[cur_index].mass_x - body.mass_x
        var dist_y = left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].sw
            end

            if left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].nw
            end

            if left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].se
            end

            if left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if left_quads[cur_index].index ~= body.index then
            var dist_x = left_quads[cur_index].mass_x - body.mass_x
            var dist_y = left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = right_quad_range[sector + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if right_quads[cur_index].type == 2 then
        var dist_x = right_quads[cur_index].mass_x - body.mass_x
        var dist_y = right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].sw
            end

            if right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].nw
            end

            if right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].se
            end

            if right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if right_quads[cur_index].index ~= body.index then
            var dist_x = right_quads[cur_index].mass_x - body.mass_x
            var dist_y = right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_quad_range[sector - sector_precision]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_quads[cur_index].type == 2 then
        var dist_x = up_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].sw
            end

            if up_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].nw
            end

            if up_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].se
            end

            if up_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_quads[cur_index].index ~= body.index then
            var dist_x = up_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_left_quad_range[sector - sector_precision - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_left_quads[cur_index].type == 2 then
        var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].sw
            end

            if up_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].nw
            end

            if up_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].se
            end

            if up_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_left_quads[cur_index].index ~= body.index then
            var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_right_quad_range[sector - sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_right_quads[cur_index].type == 2 then
        var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].sw
            end

            if up_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].nw
            end

            if up_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].se
            end

            if up_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_right_quads[cur_index].index ~= body.index then
            var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_quad_range[sector + sector_precision]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_quads[cur_index].type == 2 then
        var dist_x = down_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].sw
            end

            if down_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].nw
            end

            if down_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].se
            end

            if down_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_quads[cur_index].index ~= body.index then
            var dist_x = down_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_left_quad_range[sector + sector_precision - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_left_quads[cur_index].type == 2 then
        var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].sw
            end

            if down_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].nw
            end

            if down_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].se
            end

            if down_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_left_quads[cur_index].index ~= body.index then
            var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_right_quad_range[sector + sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_right_quads[cur_index].type == 2 then
        var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].sw
            end

            if down_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].nw
            end

            if down_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].se
            end

            if down_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_right_quads[cur_index].index ~= body.index then
            var dist_x = down_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_right_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_middle_last(bodies : region(body), roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  left_quads : region(ispace(int1d), quad), left_quad_range : region(ispace(int1d), rect1d),
  up_quads : region(ispace(int1d), quad), up_quad_range : region(ispace(int1d), rect1d),
  up_left_quads : region(ispace(int1d), quad), up_left_quad_range : region(ispace(int1d), rect1d),
  down_quads : region(ispace(int1d), quad), down_quad_range : region(ispace(int1d), rect1d),
  down_left_quads : region(ispace(int1d), quad), down_left_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(left_quads),
  reads(left_quad_range),
  reads(up_quads),
  reads(up_quad_range),
  reads(up_left_quads),
  reads(up_left_quad_range),
  reads(down_quads),
  reads(down_quad_range),
  reads(down_left_quads),
  reads(down_left_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= 0 and i ~= 1 and i ~= sector_precision and i ~= sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector_precision*(sector_precision-1) - 1]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = left_quad_range[sector_precision*(sector_precision-1) - 2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if left_quads[cur_index].type == 2 then
        var dist_x = left_quads[cur_index].mass_x - body.mass_x
        var dist_y = left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].sw
            end

            if left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].nw
            end

            if left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].se
            end

            if left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if left_quads[cur_index].index ~= body.index then
            var dist_x = left_quads[cur_index].mass_x - body.mass_x
            var dist_y = left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_quad_range[sector_precision*(sector_precision-2) - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_quads[cur_index].type == 2 then
        var dist_x = up_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].sw
            end

            if up_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].nw
            end

            if up_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].se
            end

            if up_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_quads[cur_index].index ~= body.index then
            var dist_x = up_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_left_quad_range[sector_precision*(sector_precision-2) - 2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_left_quads[cur_index].type == 2 then
        var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].sw
            end

            if up_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].nw
            end

            if up_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].se
            end

            if up_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_left_quads[cur_index].index ~= body.index then
            var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_quad_range[sector_precision * sector_precision - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_quads[cur_index].type == 2 then
        var dist_x = down_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].sw
            end

            if down_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].nw
            end

            if down_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].se
            end

            if down_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_quads[cur_index].index ~= body.index then
            var dist_x = down_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = down_left_quad_range[sector_precision * sector_precision - 2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if down_left_quads[cur_index].type == 2 then
        var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if down_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if down_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].sw
            end

            if down_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].nw
            end

            if down_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].se
            end

            if down_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = down_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if down_left_quads[cur_index].index ~= body.index then
            var dist_x = down_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = down_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * down_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = down_left_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_last_row(bodies : region(body), sector: int1d, roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  left_quads : region(ispace(int1d), quad), left_quad_range : region(ispace(int1d), rect1d),
  right_quads : region(ispace(int1d), quad), right_quad_range : region(ispace(int1d), rect1d),
  up_quads : region(ispace(int1d), quad), up_quad_range : region(ispace(int1d), rect1d),
  up_left_quads : region(ispace(int1d), quad), up_left_quad_range : region(ispace(int1d), rect1d),
  up_right_quads : region(ispace(int1d), quad), up_right_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(left_quads),
  reads(left_quad_range),
  reads(right_quads),
  reads(right_quad_range),
  reads(up_quads),
  reads(up_quad_range),
  reads(up_left_quads),
  reads(up_left_quad_range),
  reads(up_right_quads),
  reads(up_right_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= [int](sector) and i ~= [int](sector) - 1 and i ~= [int](sector) + 1 and i ~= [int](sector) + sector_precision and i ~= [int](sector) + sector_precision - 1 and i ~= [int](sector) + sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = left_quad_range[sector - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if left_quads[cur_index].type == 2 then
        var dist_x = left_quads[cur_index].mass_x - body.mass_x
        var dist_y = left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].sw
            end

            if left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].nw
            end

            if left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].se
            end

            if left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if left_quads[cur_index].index ~= body.index then
            var dist_x = left_quads[cur_index].mass_x - body.mass_x
            var dist_y = left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = right_quad_range[sector + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if right_quads[cur_index].type == 2 then
        var dist_x = right_quads[cur_index].mass_x - body.mass_x
        var dist_y = right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].sw
            end

            if right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].nw
            end

            if right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].se
            end

            if right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if right_quads[cur_index].index ~= body.index then
            var dist_x = right_quads[cur_index].mass_x - body.mass_x
            var dist_y = right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = right_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_quad_range[sector - sector_precision]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_quads[cur_index].type == 2 then
        var dist_x = up_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].sw
            end

            if up_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].nw
            end

            if up_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].se
            end

            if up_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_quads[cur_index].index ~= body.index then
            var dist_x = up_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_left_quad_range[sector - sector_precision - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_left_quads[cur_index].type == 2 then
        var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].sw
            end

            if up_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].nw
            end

            if up_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].se
            end

            if up_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_left_quads[cur_index].index ~= body.index then
            var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_right_quad_range[sector - sector_precision + 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_right_quads[cur_index].type == 2 then
        var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_right_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_right_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].sw
            end

            if up_right_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].nw
            end

            if up_right_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].se
            end

            if up_right_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_right_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_right_quads[cur_index].index ~= body.index then
            var dist_x = up_right_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_right_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_right_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_right_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task update_bodies_last(bodies : region(body), roots : region(ispace(int1d), quad),
  self_quads : region(ispace(int1d), quad), self_quad_range : region(ispace(int1d), rect1d),
  left_quads : region(ispace(int1d), quad), left_quad_range : region(ispace(int1d), rect1d),
  up_quads : region(ispace(int1d), quad), up_quad_range : region(ispace(int1d), rect1d),
  up_left_quads : region(ispace(int1d), quad), up_left_quad_range : region(ispace(int1d), rect1d))
where
  reads (bodies.{mass_x, mass_y, mass, speed_x, speed_y, index}),
  reduces +(bodies.{mass_x, mass_y, speed_x, speed_y}),
  reads(roots.{mass_x, mass_y, mass}),
  reads(self_quads),
  reads(self_quad_range),
  reads(left_quads),
  reads(left_quad_range),
  reads(up_quads),
  reads(up_quad_range),
  reads(up_left_quads),
  reads(up_left_quad_range)
do
  var traverse_list : int1d[1024]

  for body in bodies do
    var force_x = 0.0
    var force_y = 0.0

    for i=0,sector_precision*sector_precision do
      if i ~= 0 and i ~= 1 and i ~= sector_precision and i ~= sector_precision + 1 then
        var dist_x = roots[i].mass_x - body.mass_x
        var dist_y = roots[i].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        var d_force = gee * body.mass * roots[i].mass / (dist * dist * dist)
        force_x += d_force * dist_x
        force_y += d_force * dist_y
      end
    end

    var quad_range = self_quad_range[sector_precision*sector_precision-1]
    traverse_list[0] = quad_range.lo
    var traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if self_quads[cur_index].type == 2 then
        var dist_x = self_quads[cur_index].mass_x - body.mass_x
        var dist_y = self_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if self_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if self_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].sw
            end

            if self_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].nw
            end

            if self_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].se
            end

            if self_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = self_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if self_quads[cur_index].index ~= body.index then
            var dist_x = self_quads[cur_index].mass_x - body.mass_x
            var dist_y = self_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * self_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = self_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = left_quad_range[sector_precision*sector_precision-2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if left_quads[cur_index].type == 2 then
        var dist_x = left_quads[cur_index].mass_x - body.mass_x
        var dist_y = left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].sw
            end

            if left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].nw
            end

            if left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].se
            end

            if left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if left_quads[cur_index].index ~= body.index then
            var dist_x = left_quads[cur_index].mass_x - body.mass_x
            var dist_y = left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = left_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_quad_range[sector_precision*(sector_precision-1) - 1]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_quads[cur_index].type == 2 then
        var dist_x = up_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].sw
            end

            if up_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].nw
            end

            if up_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].se
            end

            if up_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_quads[cur_index].index ~= body.index then
            var dist_x = up_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_quads[cur_index].next_in_leaf
        end
      end
    end

    quad_range = up_left_quad_range[sector_precision*(sector_precision-1) - 2]
    traverse_list[0] = quad_range.lo
    traverse_index = 0

    while traverse_index >= 0 do
      var cur_index : int = traverse_list[traverse_index]
      traverse_index = traverse_index - 1

      if up_left_quads[cur_index].type == 2 then
        var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
        var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
        var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
        if dist > 1.0 then
          if up_left_quads[cur_index].size / dist >= theta then
            assert(traverse_index < 1020, "possible traverse list overflow")
            if up_left_quads[cur_index].sw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].sw
            end

            if up_left_quads[cur_index].nw ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].nw
            end

            if up_left_quads[cur_index].se ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].se
            end

            if up_left_quads[cur_index].ne ~= -1 then
              traverse_index += 1
              traverse_list[traverse_index] = up_left_quads[cur_index].ne
            end
          else
            var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
            force_x += d_force * dist_x
            force_y += d_force * dist_y
          end
        end
      else
        while cur_index ~= -1 do
          if up_left_quads[cur_index].index ~= body.index then
            var dist_x = up_left_quads[cur_index].mass_x - body.mass_x
            var dist_y = up_left_quads[cur_index].mass_y - body.mass_y
            var dist = sqrt(dist_x * dist_x + dist_y * dist_y)
            if dist > 1.0 then
              var d_force = gee * body.mass * up_left_quads[cur_index].mass / (dist * dist * dist)
              force_x += d_force * dist_x
              force_y += d_force * dist_y
            end
          end

          cur_index = up_left_quads[cur_index].next_in_leaf
        end
      end
    end

    body.mass_x += body.speed_x * delta
    body.mass_y += body.speed_y * delta
    body.speed_x += (force_x / body.mass * delta)
    body.speed_y += (force_y / body.mass * delta)
  end
end

task calculate_quad_ranges(bodies : region(body), bodies_by_sector : partition(disjoint, bodies, ispace(int1d)), quad_ranges : region(ispace(int1d), rect1d), num_quads : uint)
where
  writes(quad_ranges)
do
    var offset = 0
    for i=0,sector_precision*sector_precision do
      var current = bodies_by_sector[i]
      var quad_size_estimate = current.ispace.volume * 12 / 5
      quad_ranges[i] = rect1d({offset, offset + quad_size_estimate})
      offset += quad_size_estimate + 1
    end
end

task run_iteration(bodies : region(body), roots : region(ispace(int1d), quad), quads : region(ispace(int1d), quad), quad_ranges : region(ispace(int1d), rect1d), num_quads : uint, conf : Config, iteration : int)
where
  roots * quads,
  reads writes(bodies),
  reads writes(roots.{mass_x, mass_y, mass}),
  reads writes(quads),
  reads writes(quad_ranges)
do
  var sector_space = ispace(int1d, sector_precision * sector_precision)

  var min_x = math.huge
  var min_y = math.huge
  var max_x = -math.huge
  var max_y = -math.huge

  if iteration == 0 then
    var body_partition_index = ispace(ptr, conf.parallelism * 2)
    var bodies_partition = partition(equal, bodies, body_partition_index)

    __demand(__parallel)
    for i in body_partition_index do
      min_x min= update_boundaries_mass_x_min(bodies_partition[i])
    end

    __demand(__parallel)
    for i in body_partition_index do
      max_x max= update_boundaries_mass_x_max(bodies_partition[i])
    end

    __demand(__parallel)
    for i in body_partition_index do
      min_y min= update_boundaries_mass_y_min(bodies_partition[i])
    end
    
    __demand(__parallel)
    for i in body_partition_index do
      max_y max= update_boundaries_mass_y_max(bodies_partition[i])
    end

    __demand(__parallel)
    for i in body_partition_index do
      assign_sectors_first(bodies_partition[i], min_x, min_y, max_x, max_y)
    end
    
    __delete(bodies_partition)
  else 
    var bodies_partition = partition(bodies.sector, sector_space)

    __demand(__parallel)
    for i in sector_space do
      min_x min= update_boundaries_mass_x_min(bodies_partition[i])
    end

    __demand(__parallel)
    for i in sector_space do
      max_x max= update_boundaries_mass_x_max(bodies_partition[i])
    end

    __demand(__parallel)
    for i in sector_space do
      min_y min= update_boundaries_mass_y_min(bodies_partition[i])
    end
    
    __demand(__parallel)
    for i in sector_space do
      max_y max= update_boundaries_mass_y_max(bodies_partition[i])
    end

    __demand(__parallel)
    for i in sector_space do
      assign_sectors(bodies_partition[i], min_x, min_y, max_x, max_y)
    end
    
    __delete(bodies_partition)
  end
      
  var bodies_by_sector = partition(bodies.sector, sector_space)

  calculate_quad_ranges(bodies, bodies_by_sector, quad_ranges, num_quads)

  fill(quads.{nw, sw, ne, se, next_in_leaf}, -1)
  fill(quads.{mass_x, mass_y, mass, total, type}, 0)

  var quad_range_by_sector = partition(equal, quad_ranges, sector_space)
  var roots_by_sector = partition(equal, roots, sector_space)
  var quads_by_sector = image(quads, quad_range_by_sector, quad_ranges)

  var quads_by_sector_colors = quads_by_sector.colors
  var quads_by_sector_disjoint = dynamic_cast(partition(disjoint, quads, quads_by_sector_colors), quads_by_sector)

  __demand(__parallel)
  for i in sector_space do
    build_quad(bodies_by_sector[i], roots_by_sector[i], quads_by_sector_disjoint[i], quad_range_by_sector[i], min_x, min_y, max_x, max_y, sector_precision, conf.leaf_size, conf.max_depth, i)
  end

  update_bodies_first(bodies_by_sector[0], roots,
    quads_by_sector_disjoint[0], quad_range_by_sector[0],
    quads_by_sector_disjoint[1], quad_range_by_sector[1], 
    quads_by_sector_disjoint[sector_precision], quad_range_by_sector[sector_precision],
    quads_by_sector_disjoint[sector_precision + 1], quad_range_by_sector[sector_precision + 1])

  __demand(__parallel)
  for i=1,sector_precision do
    update_bodies_first_row(bodies_by_sector[i], i, roots,
      quads_by_sector_disjoint[i], quad_range_by_sector[i],
      quads_by_sector_disjoint[i - 1], quad_range_by_sector[i - 1],
      quads_by_sector_disjoint[i + 1], quad_range_by_sector[i + 1], 
      quads_by_sector_disjoint[i + sector_precision], quad_range_by_sector[i + sector_precision],
      quads_by_sector_disjoint[i + sector_precision - 1], quad_range_by_sector[i + sector_precision - 1],
      quads_by_sector_disjoint[i + sector_precision + 1], quad_range_by_sector[i + sector_precision + 1])
  end

  update_bodies_middle_first(bodies_by_sector[sector_precision], roots,
    quads_by_sector_disjoint[sector_precision], quad_range_by_sector[sector_precision],
    quads_by_sector_disjoint[sector_precision + 1], quad_range_by_sector[sector_precision + 1],
    quads_by_sector_disjoint[0], quad_range_by_sector[0],
    quads_by_sector_disjoint[1], quad_range_by_sector[1],
    quads_by_sector_disjoint[sector_precision*2], quad_range_by_sector[sector_precision*2],
    quads_by_sector_disjoint[sector_precision*2 + 1], quad_range_by_sector[sector_precision*2 + 1])

  __demand(__parallel)
  for i=sector_precision+1,sector_precision*(sector_precision-1)-1 do
    update_bodies_middle(bodies_by_sector[i], i, roots,
      quads_by_sector_disjoint[i], quad_range_by_sector[i],
      quads_by_sector_disjoint[i - 1], quad_range_by_sector[i - 1],
      quads_by_sector_disjoint[i + 1], quad_range_by_sector[i + 1],
      quads_by_sector_disjoint[i - sector_precision], quad_range_by_sector[i - sector_precision],
      quads_by_sector_disjoint[i - sector_precision - 1], quad_range_by_sector[i - sector_precision - 1],
      quads_by_sector_disjoint[i - sector_precision + 1], quad_range_by_sector[i - sector_precision + 1],
      quads_by_sector_disjoint[i + sector_precision], quad_range_by_sector[i + sector_precision],
      quads_by_sector_disjoint[i + sector_precision - 1], quad_range_by_sector[i + sector_precision - 1],
      quads_by_sector_disjoint[i + sector_precision + 1], quad_range_by_sector[i + sector_precision + 1])
  end

  update_bodies_middle_last(bodies_by_sector[sector_precision*(sector_precision-1)-1], roots,
    quads_by_sector_disjoint[sector_precision*(sector_precision-1)-1], quad_range_by_sector[sector_precision*(sector_precision-1)-1],
    quads_by_sector_disjoint[sector_precision*(sector_precision-1)-2], quad_range_by_sector[sector_precision*(sector_precision-1)-2],
    quads_by_sector_disjoint[sector_precision*(sector_precision-2)-1], quad_range_by_sector[sector_precision*(sector_precision-2)-1],
    quads_by_sector_disjoint[sector_precision*(sector_precision-2)-2], quad_range_by_sector[sector_precision*(sector_precision-2)-2],
    quads_by_sector_disjoint[sector_precision*sector_precision-1], quad_range_by_sector[sector_precision*sector_precision-1],
    quads_by_sector_disjoint[sector_precision*sector_precision-2], quad_range_by_sector[sector_precision*sector_precision-2])

  __demand(__parallel)
  for i=sector_precision*(sector_precision-1),sector_precision*sector_precision-1 do
    update_bodies_last_row(bodies_by_sector[i], i, roots,
      quads_by_sector_disjoint[i], quad_range_by_sector[i],
      quads_by_sector_disjoint[i - 1], quad_range_by_sector[i - 1],
      quads_by_sector_disjoint[i + 1], quad_range_by_sector[i + 1], 
      quads_by_sector_disjoint[i - sector_precision], quad_range_by_sector[i - sector_precision],
      quads_by_sector_disjoint[i - sector_precision - 1], quad_range_by_sector[i - sector_precision - 1],
      quads_by_sector_disjoint[i - sector_precision + 1], quad_range_by_sector[i - sector_precision + 1])
  end

  update_bodies_last(bodies_by_sector[sector_precision*sector_precision-1], roots,
    quads_by_sector_disjoint[sector_precision*sector_precision-1], quad_range_by_sector[sector_precision*sector_precision-1],
    quads_by_sector_disjoint[sector_precision*sector_precision-2], quad_range_by_sector[sector_precision*sector_precision-2], 
    quads_by_sector_disjoint[sector_precision*(sector_precision-1) - 1], quad_range_by_sector[sector_precision*(sector_precision-1) - 1],
    quads_by_sector_disjoint[sector_precision*(sector_precision-1) - 2], quad_range_by_sector[sector_precision*(sector_precision-1) - 2])

  __delete(bodies_by_sector)
  __delete(quad_range_by_sector)
  __delete(quads_by_sector)
  __delete(quads_by_sector_disjoint)
end

task main()
  var ts_start = c.legion_get_current_time_in_micros()
  var conf = parse_input_args()

  var input = region(ispace(ptr, conf.num_bodies), body)  
  var all_bodies = region(ispace(ptr, conf.num_bodies), body)

  attach(hdf5, input.{mass, mass_x, mass_y, speed_x, speed_y, index}, conf.input_file, regentlib.file_read_only)
  acquire(input.{mass, mass_x, mass_y, speed_x, speed_y, index})
  copy(input.{mass, mass_x, mass_y, speed_x, speed_y, index}, all_bodies.{mass, mass_x, mass_y, speed_x, speed_y, index})
  release(input.{mass, mass_x, mass_y, speed_x, speed_y, index})
  detach(hdf5, input.{mass, mass_x, mass_y, speed_x, speed_y, index})
  __delete(input)

  var sector_space = ispace(int1d, sector_precision * sector_precision)
  var roots = region(sector_space, quad)
  var quad_ranges = region(sector_space, rect1d)

  var num_quads = conf.num_bodies * 12 / 5 + sector_precision*sector_precision

  var quads = region(ispace(int1d, num_quads), quad)

  for t=0,conf.time_steps do
    run_iteration(all_bodies, roots, quads, quad_ranges, num_quads, conf, t)
  end

  __fence(__execution, __block)
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("%d\n", (ts_end - ts_start) / 1000)
end

if os.getenv('SAVEOBJ') == '1' then
  local root_dir = arg[0]:match(".*/") or "./"
  local out_dir = (os.getenv('OBJNAME') and os.getenv('OBJNAME'):match('.*/')) or root_dir
  local link_flags = terralib.newlist({"-L" .. out_dir, "-lm", "-lpmi2", "-lbarnes_hut"})
  if os.getenv('STANDALONE') == '1' then
    os.execute('cp ' .. os.getenv('LG_RT_DIR') .. '/../bindings/regent/libregent.so ' .. out_dir)
  end
  local exe = os.getenv('OBJNAME') or "barnes_hut"
  regentlib.saveobj(main, exe, "executable", cbarnes_hut.register_mappers, link_flags)
else
  regentlib.start(main, cbarnes_hut.register_mappers)
end
