import "regent"
require("barnes_hut_common")

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
local delta = 0.1
local theta = 0.5

struct Config {
  num_bodies : uint,
  random_seed : uint,
  iterations: uint
}

fspace body {
  x: double,
  y: double,
  speed_x: double,
  speed_y: double,
  mass: double,
  index: uint
}

fspace body_update {
  x: double,
  y: double,
  speed_x: double,
  speed_y: double,
  force_x: double,
  force_y: double
}

fspace boundary {
  min_x: double,
  min_y: double,
  max_x: double,
  max_y: double
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
  bodies[index] = { x = cx, y = cy, speed_x = sx, speed_y = sy, mass = mass, index = index }
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
  var speed_x_star = sx + speed * sin(angle + cmath.M_PI / 2)
  var speed_y_star = sy + speed * cos(angle + cmath.M_PI / 2)
  var mass_star = 1.0 + drand48()
  bodies[index] = { x = x_star, y = y_star, speed_x = speed_x_star, speed_y = speed_y_star, mass = mass_star, index = index }
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
    c.printf("%d: x: %f, y: %f, speed_x: %f, speed_y: %f, mass: %f\n",
    body.index, body.x, body.y, body.speed_x, body.speed_y, body.mass) 
  end
  c.printf("\n") 
end

task print_update(iteration : uint, bodies : region(body))
  where reads(bodies)
do
  c.printf("Iteration %d\n", iteration + 1)
  for body in bodies do
    c.printf("%d: x: %f, y: %f, speed_x: %f, speed_y: %f\n",
    body.index, body.x, body.y, body.speed_x, body.speed_y) 
  end
  c.printf("\n") 
end

task update_boundaries(bodies : region(body), boundaries : region(boundary))
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

task build_quad(bodies : region(body), quads : region(quad(wild)), center_x : double, center_y : double, size : double)
  where
  reads(bodies.{x, y, mass, index}),
  writes(quads),
  reads(quads)
do
  var root = dynamic_cast(ptr(quad(quads), quads), 0)
  root.center_x = center_x
  root.center_y = center_y
  root.size = size
  root.type = 2
  var index = 1
  for body in bodies do
    var body_quad = dynamic_cast(ptr(quad(quads), quads), index)
    body_quad.mass_x = body.x
    body_quad.mass_y = body.y
    body_quad.mass = body.mass
    body_quad.total = 1 
    body_quad.type = 1
    body_quad.index = body.index
    index = add_node(quads, root, body_quad, index) + 1
  end
end

task calculate_net_force(bodies : region(body), body_updates : region(body_update), quads : region(quad(wild)), index : uint, node : ptr(quad(wild), quads)): uint
  where
  reads(bodies),
  reduces +(body_updates.{force_x, force_y}),
  reads(quads)
do
  -- c.printf("index %d", index)
  if node.type == 1 and node.index == index then
    -- c.printf("nope", index)
    return 1
  end

  var body = bodies[index]
  var dist = sqrt((body.x - node.mass_x) * (body.x - node.mass_x) + (body.y - node.mass_y) * (body.y - node.mass_y))

  if node.type == 2 and node.size / dist >= theta then
    -- c.printf("threshold hit: %f dist: %f", node.size / dist, dist)
    calculate_net_force(bodies, body_updates, quads, index, node.sw)
    calculate_net_force(bodies, body_updates, quads, index, node.nw)
    calculate_net_force(bodies, body_updates, quads, index, node.se)
    calculate_net_force(bodies, body_updates, quads, index, node.ne)
    return 1
  end

  var d_force = gee * body.mass * node.mass / (dist * dist)
  var xn = (node.mass_x - body.x) / dist
  var yn = (node.mass_y - body.y) / dist
  var d_force_x = d_force * xn
  var d_force_y = d_force * yn

  body_updates[index].force_x += d_force_x
  body_updates[index].force_y += d_force_y
  return 1
end

task update_body_positions(bodies : region(body), body_updates : region(body_update), quads : region(quad(wild)))
  where
  reads(bodies),
  reads(body_updates.{force_x, force_y}),
  writes(body_updates.{x, y, speed_x, speed_y}),
  reduces +(body_updates.{force_x, force_y}),
  reads(quads)
do
  var root = dynamic_cast(ptr(quad(quads), quads), 0)
  for body in bodies do
    calculate_net_force(bodies, body_updates, quads, body.index, root)
    body_updates[body.index].x = body.x + body.speed_x * delta
    body_updates[body.index].y = body.y + body.speed_y * delta
    body_updates[body.index].speed_x = body.speed_x + body_updates[body.index].force_x / body.mass * delta
    body_updates[body.index].speed_y = body.speed_y + body_updates[body.index].force_y / body.mass * delta
  end
end

task run_iteration(bodies : region(body), body_updates : region(body_update), body_index : ispace(ptr))
  where
  reads(bodies),
  reads(body_updates.{force_x, force_y}),
  writes(body_updates.{x, y, speed_x, speed_y}),
  reduces +(body_updates.{force_x, force_y})
do
  var boundaries_index = ispace(ptr, 1)
  var boundaries = region(boundaries_index, boundary) 
  boundaries[0] = { min_x = bodies[0].x, min_y = bodies[0].y, max_x = bodies[0].x, max_y = bodies[0].y }
  
  var bodies_partition = partition(equal, bodies, body_index)
  for i in body_index do
    update_boundaries(bodies_partition[i], boundaries)
  end

  c.printf("boundaries: min_x=%f min_y=%f max_x=%f max_y=%f\n\n", boundaries[0].min_x, boundaries[0].min_y, boundaries[0].max_x, boundaries[0].max_y)

  var size_x = boundaries[0].max_x - boundaries[0].min_x
  var size_y = boundaries[0].max_y - boundaries[0].min_y
  var size = max(size_x, size_y)

  var quads = region(ispace(ptr, 1000000), quad(quads))
  fill(quads.{nw, sw, ne, se}, null(ptr(quad(quads), quads)))

  build_quad(bodies, quads, boundaries[0].min_x + size / 2, boundaries[0].min_y + size / 2, size)

  var body_updates_partition = partition(equal, body_updates, body_index)
  for i in body_index do
    update_body_positions(bodies_partition[i], body_updates_partition[i], quads)
  end
end

task main()
  var conf : Config
  conf.num_bodies = 16
  conf.random_seed = 213
  conf.iterations = 1

  conf = parse_input_args(conf)
  c.printf("settings: bodies=%d seed=%d\n\n", conf.num_bodies, conf.random_seed) 

  var body_index = ispace(ptr, conf.num_bodies)
  var bodies = region(body_index, body)
  var body_updates = region(body_index, body_update)

  init_2_galaxies(bodies, conf)

  print_bodies_initial(bodies)

  for i=0,conf.iterations do
      fill(body_updates.{force_x, force_y}, 0)
      run_iteration(bodies, body_updates, body_index)
      copy(body_updates.{x, y, speed_x, speed_y}, bodies.{x, y, speed_x, speed_y})
      print_update(i, bodies)
  end
end
regentlib.start(main)
