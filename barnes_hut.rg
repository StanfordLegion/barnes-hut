import "regent"

local c = regentlib.c
local cos = regentlib.cos(float)
local sin = regentlib.sin(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")
local cstring = terralib.includec("string.h")
local std = terralib.includec("stdlib.h")

rawset(_G, "drand48", std.drand48)
rawset(_G, "srand48", std.srand48)

local gee = 100
local sector_precision = 16

struct Config {
  num_bodies : uint,
  random_seed : uint,
  iterations: uint
}

fspace sector {
  x: uint,
  y: uint
}

fspace body {
  x: double,
  y: double,
  x_speed: double,
  y_speed: double,
  mass: double,
  index: uint,
  sector: int2d
}

fspace boundary {
  min_x: double,
  min_y: double,
  max_x: double,
  max_y: double
}

fspace quad(r : region(quad(wild))) {
  mass_x: double,
  mass_y: double,
  mass: double,
  center_x: double,
  center_y: double,
  size: double,
  total: uint,
  type: uint,
  ne: ptr(quad(wild), r),
  nw: ptr(quad(wild), r),
  se: ptr(quad(wild), r),
  sw: ptr(quad(wild), r)
}

struct quad_str {
  mass_x: double,
  mass_y: double,
  mass: double,
  total: uint,
  ne: &quad_str,
  nw: &quad_str,
  se: &quad_str,
  sw: &quad_str,

  center_x: double,
  center_y: double,
  size: double,

  type: uint
}

terra parse_input_args(conf : Config)
  var args = c.legion_runtime_get_input_args()
  for i = 0, args.argc do
    if cstring.strcmp(args.argv[i], "-b") == 0 then
      i = i + 1
      conf.num_bodies = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-s") == 0 then
      i = i + 1
      conf.random_seed = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      conf.iterations = std.atoi(args.argv[i])
    end
  end
  return conf
end

task init_black_hole(bodies : region(ispace(ptr), body), mass : uint, cx : double, cy : double, sx : double, sy : double, index: uint)
  where writes(bodies)
do
  bodies[index] = { x = cx, y = cy, x_speed = sx, y_speed = sy, mass = mass, index = index, sector = {x = 0, y = 0} }
end

task init_star(bodies : region(ispace(ptr), body), num : uint, max_radius : double, cx : double, cy : double, sx : double, sy : double, index: uint)
  where writes(bodies)
do
  var total_m = 1.5 * num
  var cube_max_radius = max_radius * max_radius * max_radius
  
  var angle = drand48() * 2 * cmath.M_PI
  var radius = 25 + max_radius * drand48()
  var x_star = cx + radius * sin(angle)
  var y_star = cy + radius * cos(angle)
  var speed = sqrt(gee * num / radius + gee * total_m * radius * radius / cube_max_radius)
  var x_speed_star = sx + speed * sin(angle + cmath.M_PI / 2)
  var y_speed_star = sy + speed * cos(angle + cmath.M_PI / 2)
  var mass_star = 1.0 + drand48()
  bodies[index] = { x = x_star, y = y_star, x_speed = x_speed_star, y_speed = y_speed_star, mass = mass_star, index = index, sector = {x = 0, y = 0} }
end

task init_2_galaxies(bodies : region(body), conf : Config)
  where writes(bodies)
do
  srand48(conf.random_seed)

  var bodies_partition = partition(equal, bodies, ispace(ptr, conf.num_bodies))

  var num1 = conf.num_bodies / 8
  init_black_hole(bodies_partition[0], num1, 0, 0, 0, 0, 0)
  
  __demand(__parallel)
  for i = 1, num1 do
    init_star(bodies_partition[i], num1, 300, 0, 0, 0, 0, i)
  end

  var num2 = conf.num_bodies / 8 * 7
  init_black_hole(bodies_partition[num1], num2, -1800, -1200, 0, 0, num1)
  
  __demand(__parallel)
  for i = num1 + 1, num1 + num2 do
    init_star(bodies_partition[i], num2, 350, -1800, -1200, 0, 0, i)
  end

end

task print_bodies_initial(bodies : region(body))
  where reads(bodies)
do
  c.printf("Initial bodies:\n") 
  for body in bodies do
    c.printf("%d: x: %f, y: %f, x_speed: %f, y_speed: %f, mass: %f\n",
    body.index, body.x, body.y, body.x_speed, body.y_speed, body.mass) 
  end
  c.printf("\n") 
end

task add_node(quads : region(quad(wild)), cur: ptr(quad(wild), quads), body: ptr(quad(wild), quads), last_used: uint): uint
where
  reads(quads),
  writes(quads)
do
  c.printf("Inserting cur: %d, body: %d, last_used: %d\n", cur, body, last_used)
  var half_size = cur.size / 2
  var new_last_used = last_used
  if body.mass_x <= cur.center_x then
    if body.mass_y <= cur.center_y then
      if cur.sw == dynamic_cast(ptr(quad(quads), quads), 0) then
        cur.sw = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.sw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x - half_size / 2
        new_fork.center_y = cur.center_y - half_size / 2
        new_fork.size = half_size
        cur.sw = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.sw), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.sw, body, last_used)
      end
    else
      if cur.nw == dynamic_cast(ptr(quad(quads), quads), 0) then
        cur.nw = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.nw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x - half_size / 2
        new_fork.center_y = cur.center_y + half_size / 2
        new_fork.size = half_size
        cur.nw = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.nw), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.nw, body, last_used)
      end      
    end
  else
    if body.mass_y <= cur.center_y then
      if cur.se == dynamic_cast(ptr(quad(quads), quads), 0) then
        cur.se = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.se).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x + half_size / 2
        new_fork.center_y = cur.center_y - half_size / 2
        new_fork.size = half_size
        cur.se = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.se), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.se, body, last_used)
      end
    else
      if cur.ne == dynamic_cast(ptr(quad(quads), quads), 0) then
        cur.ne = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.sw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x + half_size / 2
        new_fork.center_y = cur.center_x + half_size / 2
        new_fork.size = half_size
        cur.ne = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.ne), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.ne, body, last_used)
      end      
    end
  end  

  return new_last_used
end

task build_quad(bodies: region(body), quads: region(quad(wild)), from_x: double, from_y: double, sector_size: double, start_index: int)
  where
  reads(bodies.{x, y, mass}),
  writes(quads),
  reads(quads)
do
  var index = start_index
  var root = dynamic_cast(ptr(quad(quads), quads), index)
  root.center_x = from_x + sector_size / 2
  root.center_y = from_y + sector_size / 2
  root.size = sector_size
  root.type = 2
  index = index + 1
  for body in bodies do
    var body_quad = dynamic_cast(ptr(quad(quads), quads), index)
    body_quad.mass_x = body.x
    body_quad.mass_y = body.y
    body_quad.mass = body.mass
    body_quad.total = 1 
    body_quad.type = 1   
    index = index + 1
    index = add_node(quads, root, body_quad, index)
  end
end

task assign_sectors(bodies: region(body), min_x: double, min_y: double, size_x: double, size_y: double)
  where
  reads(bodies.{x, y, sector, index}),
  writes(bodies.sector)
do
  for body in bodies do
    var sector_x: int64 = cmath.floor((body.x - min_x) / size_x)
    if (sector_x >= sector_precision) then
      sector_x = sector_x - 1
    end

    var sector_y: int64 = cmath.floor((body.y - min_y) / size_y)
    if (sector_y >= sector_precision) then
      sector_y = sector_y - 1
    end

    body.sector = { x = sector_x , y = sector_y }

    c.printf("x: %d, y: %d\n", sector_x, sector_y)
  end
end

task update_boundaries(bodies: region(body), boundaries: region(boundary))
  where
  reads(bodies.{x, y}),
  reads(boundaries),
  reduces min(boundaries.{min_x, min_y}),
  reduces max(boundaries.{max_x, max_y})
do
  for body in bodies do
    boundaries[0].min_x min = min(body.x, boundaries[0].min_x)
    boundaries[0].min_y min = min(body.y, boundaries[0].min_y)
    boundaries[0].max_x max = max(body.x, boundaries[0].max_x)
    boundaries[0].max_y max = max(body.y, boundaries[0].max_y)
  end
end

task run_iteration(bodies : region(body), body_index : ispace(ptr))
  where
  reads(bodies),
  writes(bodies.sector)
do
  var boundaries_index = ispace(ptr, 1)
  var boundaries = region(boundaries_index, boundary) 
  boundaries[0] = { min_x = bodies[0].x, min_y = bodies[0].y, max_x = bodies[0].x, max_y = bodies[0].y }
  
  var bodies_partition = partition(equal, bodies, body_index)
  for i in body_index do
    update_boundaries(bodies_partition[i], boundaries)
  end

  c.printf("boundaries: min_x=%f min_y=%f max_x=%f max_y=%f\n", boundaries[0].min_x, boundaries[0].min_y, boundaries[0].max_x, boundaries[0].max_y)

  var size_x = (boundaries[0].max_x - boundaries[0].min_x) / sector_precision
  var size_y = (boundaries[0].max_y - boundaries[0].min_y) / sector_precision

  for i in body_index do
    assign_sectors(bodies_partition[i], boundaries[0].min_x, boundaries[0].min_y, size_x, size_y)
  end

  var sector_index = ispace(int2d, { x = sector_precision, y = sector_precision })
  var bodies_by_sector = partition(bodies.sector, sector_index)

  var max_size = sector_precision * sector_precision * 100
  var quads = region(ispace(ptr, max_size), quad(quads))
  var quad_partitions = partition(equal, quads, sector_index)

  for i in bodies_by_sector.colors do
    build_quad(bodies_by_sector[i], quad_partitions[i], 1, 1, 1, max_size / (sector_precision * sector_precision) * (i.x * sector_precision + i.y))
  end

  --[[
  c.printf("\n")
  for i in quads do
    c.printf("total: %d\n", i.total)
  end
  --]]
end

task main()
  var conf : Config
  conf.num_bodies = 16
  conf.random_seed = 213
  conf.iterations = 5

  conf = parse_input_args(conf)
  c.printf("circuit settings: bodies=%d seed=%d\n", conf.num_bodies, conf.random_seed) 

  var body_index = ispace(ptr, conf.num_bodies)
  var bodies = region(body_index, body)

  init_2_galaxies(bodies, conf)

  print_bodies_initial(bodies)

  run_iteration(bodies, body_index)
  
end
regentlib.start(main)
