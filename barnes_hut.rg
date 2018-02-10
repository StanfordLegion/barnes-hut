import "regent"
require("quad_tree")

local assert = regentlib.assert
local c = regentlib.c
local cos = regentlib.cos(float)
local sin = regentlib.sin(float)
local sqrt = regentlib.sqrt(float)

local cmath = terralib.includec("math.h")
local cstring = terralib.includec("string.h")
local std = terralib.includec("stdlib.h")

local BarnesHitIO = require("barnes_hut_io")
local QuadTreeSizer = require("quad_tree_sizer")

rawset(_G, "drand48", std.drand48)
rawset(_G, "srand48", std.srand48)

local gee = 100
local delta = 0.1
local theta = 0.5
local N = 4
local sector_precision = 16

struct Config {
  num_bodies : uint,
  random_seed : uint,
  iterations : uint,
  verbose : bool,
  output_dir : rawstring,
  output_dir_set : bool
}

fspace body {
  {mass_x, mass_y, speed_x, speed_y, mass, force_x, force_y} : double,
  sector : int1d,
  {color, index} : uint,
}

fspace boundary {
  {min_x, min_y, max_x, max_y} : double
}

terra parse_input_args(conf : Config)
  var args = c.legion_runtime_get_input_args()
  for i = 0, args.argc do
    if cstring.strcmp(args.argv[i], "-o") == 0 then
      i = i + 1
      conf.output_dir = args.argv[i]
      conf.output_dir_set = true
    elseif cstring.strcmp(args.argv[i], "-b") == 0 then
      i = i + 1
      conf.num_bodies = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-s") == 0 then
      i = i + 1
      conf.random_seed = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      conf.iterations = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-v") == 0 then
      conf.verbose = true
    end
  end
  return conf
end

task init_black_hole(bodies : region(ispace(ptr), body), mass : uint, cx : double, cy : double, sx : double, sy : double, color : uint, index : uint)
  where writes(bodies)
do
  for body in bodies do
    body.mass_x = cx
    body.mass_y = cy
    body.speed_x = sx
    body.speed_y = sy
    body.mass = mass
    body.color = color
    body.index = index
  end
end

task init_star(bodies : region(ispace(ptr), body), num : uint, max_radius : double, cx : double, cy : double, sx : double, sy : double, color : uint, index : uint)
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

  for body in bodies do
    body.mass_x = x_star
    body.mass_y = y_star
    body.speed_x = speed_x_star
    body.speed_y = speed_y_star
    body.mass = mass_star
    body.color = color
    body.index = index
  end
end

task init_2_galaxies(bodies : region(body), conf : Config)
  where writes(bodies)
do
  srand48(conf.random_seed)

  var bodies_partition = partition(equal, bodies, ispace(ptr, conf.num_bodies))

  var num1 = conf.num_bodies / 8
  init_black_hole(bodies_partition[0], num1, 0, 0, 0, 0, 0, 0)
  
  __demand(__parallel)
  for i = 1, num1 do
    init_star(bodies_partition[i], num1, 300, 0, 0, 0, 0, 1, i)
  end

  var num2 = conf.num_bodies / 8 * 7
  init_black_hole(bodies_partition[num1], num2, -1800, -1200, 0, 0, 0, num1)
  
  __demand(__parallel)
  for i = num1 + 1, num1 + num2 do
    init_star(bodies_partition[i], num2, 350, -1800, -1200, 0, 0, 2, i)
  end
end

task print_bodies_initial(bodies : region(body))
  where reads(bodies)
do
  c.printf("Initial bodies:\n") 
  for body in bodies do
    c.printf("%d: x: %f, y: %f, speed_x: %f, speed_y: %f, mass: %f\n",
    body.index, body.mass_x, body.mass_y, body.speed_x, body.speed_y, body.mass)
  end
  c.printf("\n") 
end

task print_update(iteration : uint, bodies : region(body))
  where reads(bodies)
do
  c.printf("Iteration %d\n", iteration + 1)
  for body in bodies do
    var sector_x = body.sector % sector_precision
    var sector_y: int64 = cmath.floor(body.sector / sector_precision)

    c.printf("%d: x: %f, y: %f, speed_x: %f, speed_y: %f, sector: (%d, %d)\n",
    body.index, body.mass_x, body.mass_y, body.speed_x, body.speed_y, sector_x, sector_y)
  end
  c.printf("\n") 
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

task assign_sectors(bodies : region(body), min_x : double, min_y : double, size : double)
  where
  reads(bodies.{mass_x, mass_y, sector}),
  writes(bodies.sector)
do
  for body in bodies do
    var sector_x: int64 = cmath.floor((body.mass_x - min_x) / (size / sector_precision))
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

task size_quad(bodies : region(body), center_x : double, center_y : double, size : double, max_size: region(int))
  where reads(bodies.{mass_x, mass_y, index}),
  reads (max_size),
  reduces max(max_size)
do
  var root = create_placeholder()
  root.center_x = center_x
  root.center_y = center_y
  root.size = size
  root.type = 2

  for body in bodies do
    var body_quad = create_placeholder()
    body_quad.mass_x = body.mass_x
    body_quad.mass_y = body.mass_y
    body_quad.type = 1
    add_placeholder(root, body_quad)
  end

  max_size[0] max = max(max_size[0], count(root, true))
end

task build_quad(bodies : region(body), quads : region(quad(wild)), center_x : double, center_y : double, size : double)
  where
  reads(bodies.{mass_x, mass_y, mass, index}),
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
    assert(body_quad.type == 0, "region already allocated")
    body_quad.mass_x = body.mass_x
    body_quad.mass_y = body.mass_y
    body_quad.mass = body.mass
    body_quad.total = 1 
    body_quad.type = 1
    body_quad.index = body.index
    index = add_node(quads, root, body_quad, index) + 1
  end
end

task calculate_net_force(bodies : region(body), quads : region(quad(wild)), node : ptr(quad(wild), quads), body_ptr : ptr(body, bodies)): uint
  where
  reads(bodies),
  reduces +(bodies.{force_x, force_y}),
  reads(quads)
do
  if node.type == 1 and node.index == body_ptr.index then
    return 1
  end

  var dist = sqrt((body_ptr.mass_x - node.mass_x) * (body_ptr.mass_x - node.mass_x) + (body_ptr.mass_y - node.mass_y) * (body_ptr.mass_y - node.mass_y))
  if dist == 0 then
    return 1
  end

  if node.type == 2 and node.size / dist >= theta then
    calculate_net_force(bodies, quads, node.sw, body_ptr)
    calculate_net_force(bodies, quads, node.nw, body_ptr)
    calculate_net_force(bodies, quads, node.se, body_ptr)
    calculate_net_force(bodies, quads, node.ne, body_ptr)
    return 1
  end

  var d_force = gee * body_ptr.mass * node.mass / (dist * dist)
  var xn = (node.mass_x - body_ptr.mass_x) / dist
  var yn = (node.mass_y - body_ptr.mass_y) / dist
  var d_force_x = d_force * xn
  var d_force_y = d_force * yn

  body_ptr.force_x += d_force_x
  body_ptr.force_y += d_force_y
  return 1
end

task update_body_positions(bodies : region(body), quads : region(quad(wild)))
  where
  reads(bodies),
  writes(bodies),
  reads(quads)
do
  var root = dynamic_cast(ptr(quad(quads), quads), 0)
  for body in bodies do
    calculate_net_force(bodies, quads, root, body)
    body.mass_x = body.mass_x + body.speed_x * delta
    body.mass_y = body.mass_y + body.speed_y * delta
    body.speed_x = body.speed_x + body.force_x / body.mass * delta
    body.speed_y = body.speed_y + body.force_y / body.mass * delta
  end
end

task run_iteration(bodies : region(body), body_index : ispace(ptr), boundaries : region(boundary), verbose: bool)
  where
  reads(bodies),
  writes(bodies),
  reads(boundaries),
  writes(boundaries)
do
  boundaries[0] = { min_x = bodies[0].mass_x, min_y = bodies[0].mass_y, max_x = bodies[0].mass_x, max_y = bodies[0].mass_y }
  
  var bodies_partition = partition(equal, bodies, body_index)
  for i in body_index do
    update_boundaries(bodies_partition[i], boundaries)
  end

  var size_x = boundaries[0].max_x - boundaries[0].min_x
  var size_y = boundaries[0].max_y - boundaries[0].min_y
  var size = max(size_x, size_y)

  for i in body_index do
    assign_sectors(bodies_partition[i], boundaries[0].min_x, boundaries[0].min_y, size)
  end

  if verbose then
    c.printf("Calculating required size of quad tree\n")
  end
  
  var sector_index = ispace(int1d, sector_precision * sector_precision + 1)
  var bodies_by_sector = partition(bodies.sector, sector_index)
  var max_size = region(ispace(ptr, 1), int)
  
  for i in sector_index do
    var sector_x = i % sector_precision
    var sector_y: int64 = cmath.floor(i / sector_precision)
    var center_x = boundaries[0].min_x + (sector_x + 0.5) * size / sector_precision
    var center_y = boundaries[0].min_y + (sector_y + 0.5) * size / sector_precision
    size_quad(bodies_by_sector[i], center_x, center_y, size, max_size)
  end

  var quad_size = max_size[0] * sector_precision * sector_precision
  
  if verbose then
    c.printf("Quad tree size: %d\n", quad_size)
  end

  var quads = region(ispace(ptr, quad_size), quad(quads))
  fill(quads.{nw, sw, ne, se}, null(ptr(quad(quads), quads)))
  fill(quads.type, 0)

  for i in sector_index do
    var sector_x = i % sector_precision
    var sector_y: int64 = cmath.floor(i / sector_precision)
    var center_x = boundaries[0].min_x + (sector_x + 0.5) * size / sector_precision
    var center_y = boundaries[0].min_y + (sector_y + 0.5) * size / sector_precision
    size_quad(bodies_by_sector[i], center_x, center_y, size, max_size)
  end

  build_quad(bodies, quads, boundaries[0].min_x + size / 2, boundaries[0].min_y + size / 2, size)

  for i in body_index do
    update_body_positions(bodies_partition[i], quads)
  end
end

task main()
  var conf : Config
  conf.num_bodies = 128
  conf.random_seed = 213
  conf.iterations = 10
  conf.output_dir_set = false

  conf = parse_input_args(conf)
  
  if not conf.output_dir_set then
    c.printf("output dir must be specified\n")
    c.exit(-1)
  end
  
  if conf.verbose then
    c.printf("settings: output dir=%s bodies=%d seed=%d\n\n", conf.output_dir, conf.num_bodies, conf.random_seed)
  end

  var body_index = ispace(ptr, conf.num_bodies)
  var bodies = region(body_index, body)

  init_2_galaxies(bodies, conf)

  if conf.verbose then
    print_bodies_initial(bodies)
  end

  for i=0,conf.iterations do
      var boundaries = region(ispace(ptr, 1), boundary)
      fill(bodies.{force_x, force_y}, 0)
      run_iteration(bodies, body_index, boundaries, conf.verbose)

      var boundary = boundaries[0]
      if conf.verbose then
        c.printf("boundaries: min_x=%f min_y=%f max_x=%f max_y=%f\n\n", boundary.min_x, boundary.min_y, boundary.max_x, boundary.max_y)
        print_update(i, bodies)
      end

      var fp = open(i, conf.output_dir)
      c.fprintf(fp, "<svg viewBox=\"0 0 850 850\" xmlns=\"http://www.w3.org/2000/svg\">")

      var size_x = boundary.max_x - boundary.min_x
      var size_y = boundary.max_y - boundary.min_y
      var size = max(size_x, size_y)
      var scale = 800.0 / size

      for body in bodies do
        var color = "black"
        if body.color == 1 then
          color = "blue"
        elseif body.color == 2 then
          color = "orange"
        end
        c.fprintf(fp, "<circle cx=\"%f\" cy=\"%f\" r=\"10\" fill=\"%s\" />", (body.mass_x - boundary.min_x) * scale + 25,  (body.mass_y - boundary.min_y) * scale + 25, color)
      end

      c.fprintf(fp, "</svg>")
      c.fclose(fp)
  end
end
regentlib.start(main)
