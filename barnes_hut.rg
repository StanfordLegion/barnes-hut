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

struct Config {
  num_bodies : uint,
  random_seed : uint,
  iterations: uint
}

fspace body {
  x: double,
  y: double,
  x_speed: double,
  y_speed: double,
  mass: double,
  index: uint
}

fspace boundary {
  minX: double,
  minY: double,
  maxX: double,
  maxY: double
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

task update_boundaries(bodies: region(ispace(ptr), body), boundaries: region(ispace(ptr), boundary))
  where
  reads(bodies.{x, y}),
  reads(boundaries),
  reduces min(boundaries.{minX, minY}),
  reduces max(boundaries.{maxX, maxY})
do
  for body in bodies do
    boundaries[0].minX min = min(body.x, boundaries[0].minX)
    boundaries[0].minY min = min(body.y, boundaries[0].minY)
    boundaries[0].maxX max = max(body.x, boundaries[0].maxX)
    boundaries[0].maxY max = max(body.y, boundaries[0].maxY)
  end
end

task init_galaxy(bodies : region(ispace(ptr), body), from : uint, num : uint, max_radius : double, cx : double, cy : double, sx : double, sy : double)
  where writes(bodies)
do
  var total_m = 1.5 * num
  var mass_black_hole = num
  var cube_max_radius = max_radius * max_radius * max_radius
  
  bodies[from] = { x = cx, y = cy, x_speed = sx, y_speed = sy, mass = mass_black_hole, index = from }

  for i = from + 1, from + num do
    var angle = drand48() * 2 * cmath.M_PI
    var radius = 25 + max_radius * drand48()
    var x_star = cx + radius * sin(angle)
    var y_star = cy + radius * cos(angle)
    var speed = sqrt(gee * mass_black_hole / radius + gee * total_m * radius * radius / cube_max_radius)
    var x_speed_star = sx + speed * sin(angle + cmath.M_PI / 2)
    var y_speed_star = sy + speed * cos(angle + cmath.M_PI / 2)
    var mass_star = 1.0 + drand48()
    bodies[i] = { x = x_star, y = y_star, x_speed = x_speed_star, y_speed = y_speed_star, mass = mass_star, index = i }
  end
end

task init_2_galaxies(bodies : region(ispace(ptr), body), conf : Config)
  where writes(bodies)
do
  srand48(conf.random_seed)

  init_galaxy(bodies, 0, conf.num_bodies / 8, 300, 0, 0, 0, 0)
  init_galaxy(bodies, conf.num_bodies / 8, conf.num_bodies / 8 * 7, 350, -1800, -1200, 0, 0)
end

task print_bodies_initial(bodies : region(ispace(ptr), body))
  where reads(bodies)
do
  c.printf("Initial bodies:\n") 
  for body in bodies do
    c.printf("%d: x: %f, y: %f, x_speed: %f, y_speed: %f, mass: %f\n",
    body.index, body.x, body.y, body.x_speed, body.y_speed, body.mass) 
  end
  c.printf("\n") 
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

  var boundaries_index = ispace(ptr, 1)
  var boundaries = region(boundaries_index, boundary) 
  boundaries[0] = { minX = bodies[0].x, minY = bodies[0].y, maxX = bodies[0].x, maxY = bodies[0].y }
  
  var bodies_partition = partition(equal, bodies, body_index)
  for i in body_index do
    update_boundaries(bodies_partition[i], boundaries)
  end

  c.printf("boundaries: minX=%f minY=%f maxX=%f maxY=%f\n", boundaries[0].minX, boundaries[0].minY, boundaries[0].maxX, boundaries[0].maxY) 
end
regentlib.start(main)