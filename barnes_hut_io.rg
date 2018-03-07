import "regent"
require("barnes_hut_common")

local c = regentlib.c
local cstring = terralib.includec("string.h")
local std = terralib.includec("stdlib.h")

struct Config {
  time_steps : uint,
  input_file : rawstring,
  input_file_set : bool,
  svg_dir : rawstring,
  svg_dir_set : bool,
  csv_dir : rawstring,
  csv_dir_set : bool,
  parallelism : uint,
  N : uint,
  leaf_size : uint
  fixed_partition_size : uint,
  max_depth : uint,
}

terra parse_input_args()
  var conf : Config
  conf.time_steps = 10
  conf.input_file_set = false
  conf.csv_dir_set = false
  conf.svg_dir_set = false
  conf.leaf_size = 32
  conf.N = 4
  conf.parallelism = 16
  conf.fixed_partition_size = -1
  conf.max_depth = -1

  var args = c.legion_runtime_get_input_args()

  for i = 0, args.argc do
    if cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      conf.input_file = args.argv[i]
      conf.input_file_set = true
    elseif cstring.strcmp(args.argv[i], "-c") == 0 then
      i = i + 1
      conf.csv_dir = args.argv[i]
      conf.csv_dir_set = true
    elseif cstring.strcmp(args.argv[i], "-s") == 0 then
      i = i + 1
      conf.svg_dir = args.argv[i]
      conf.svg_dir_set = true
    elseif cstring.strcmp(args.argv[i], "-t") == 0 then
      i = i + 1
      conf.time_steps = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      conf.parallelism = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-N") == 0 then
      i = i + 1
      conf.N = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-l") == 0 then
      i = i + 1
      conf.leaf_size = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-x") == 0 then
      i = i + 1
      conf.fixed_partition_size = std.atoi(args.argv[i])
    elseif cstring.strcmp(args.argv[i], "-m") == 0 then
      i = i + 1
      conf.max_depth = std.atoi(args.argv[i])
    end
  end

  if not conf.input_file_set then
    c.printf("input file must be set")
    c.exit(-1)
  end

  c.printf("settings: time_steps=%d parallelism=%d N=%d leaf_size=%d input_file=%s", conf.time_steps, conf.parallelism, conf.N, conf.leaf_size, conf.input_file)

  if conf.csv_dir_set then
    c.printf(" csv_dir=%s", conf.csv_dir)
  end

  if conf.svg_dir_set then
    c.printf(" svg_dir=%s", conf.svg_dir)
  end

  c.printf("\n\n")

  return conf
end

terra get_number_of_bodies(conf : Config)
  var fp = c.fopen(conf.input_file, "r")
  var line : int8[1024]
  var num_lines = 0
  
  while c.fgets(line, 1024, fp) ~= nil do
     num_lines = num_lines + 1
  end

  c.fclose(fp)

  return num_lines
end

task load_bodies(bodies : region(body), conf : Config, num_bodies : int)
  where writes(bodies.{index, mass_x, mass_y, speed_x, speed_y, mass})
do
  var index = 0
  var mass_x = 0.0
  var mass_y = 0.0
  var speed_x = 0.0
  var speed_y = 0.0
  var mass = 0.0

  var fp = c.fopen(conf.input_file, "r")
  var line : int8[1024]

  for i=0,num_bodies do
    c.fgets(line, 1024, fp)
    var index = std.atoi(cstring.strtok(line, ","))
    bodies[index].index = index
    bodies[index].mass_x = std.atof(cstring.strtok([&int8](0), ","))
    bodies[index].mass_y = std.atof(cstring.strtok([&int8](0), ","))
    bodies[index].speed_x = std.atof(cstring.strtok([&int8](0), ","))
    bodies[index].speed_y = std.atof(cstring.strtok([&int8](0), ","))
    bodies[index].mass = std.atof(cstring.strtok([&int8](0), ","))
  end

  c.fclose(fp)
end

task print_bodies_csv_initial(bodies : region(body), conf : Config)
  where reads(bodies.{index, mass_x, mass_y, speed_x, speed_y, mass})
do
  var output_path : int8[1000]
  c.sprintf([&int8](output_path), "%s/0.csv", conf.csv_dir)

  var fp = c.fopen(output_path, "w")
  for body in bodies do
    c.fprintf(fp, "%d,%f,%f,%f,%f,%f\n", body.index, body.mass_x, body.mass_y, body.speed_x, body.speed_y, body.mass)
  end

  c.fclose(fp)
end

task print_bodies_csv_update(bodies : region(body), conf : Config, time_step : uint)
  where reads(bodies.{index, mass_x, mass_y, speed_x, speed_y, eliminated})
do
  var output_path : int8[1000]
  c.sprintf([&int8](output_path), "%s/%d.csv", conf.csv_dir, time_step)

  var fp = c.fopen(output_path, "w")
  for body in bodies do
    if [int](body.eliminated) == 0 then
      c.fprintf(fp, "%d,%f,%f,%f,%f\n", body.index, body.mass_x, body.mass_y, body.speed_x, body.speed_y)
    end
  end

  c.fclose(fp)
end

task print_bodies_svg(bodies : region(body), boundaries : region(boundary), conf : Config, time_step : uint)
  where
    reads(bodies.{index, mass_x, mass_y, speed_x, speed_y}),
    reads(boundaries)
do
  var output_path : int8[1000]
  c.sprintf([&int8](output_path), "%s/%d.svg", conf.svg_dir, time_step)

  var fp = c.fopen(output_path, "w")
  c.fprintf(fp, "<svg viewBox=\"0 0 850 850\" xmlns=\"http://www.w3.org/2000/svg\">")

  var boundary = boundaries[0]
  var size_x = boundary.max_x - boundary.min_x
  var size_y = boundary.max_y - boundary.min_y
  var size = max(size_x, size_y)
  var scale = 800.0 / size

  for body in bodies do
    c.fprintf(fp, "<circle cx=\"%f\" cy=\"%f\" r=\"10\" fill=\"blue\" />", (body.mass_x - boundary.min_x) * scale + 25,  (body.mass_y - boundary.min_y) * scale + 25)
  end

  c.fprintf(fp, "</svg>")
  c.fclose(fp)
end
