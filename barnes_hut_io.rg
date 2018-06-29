import "regent"
require("barnes_hut_common")

local c = regentlib.c
local cstring = terralib.includec("string.h")
local std = terralib.includec("stdlib.h")

struct Config {
  time_steps : uint,
  input_file : rawstring,
  input_file_set : bool,
  num_bodies : uint,
  svg_dir : rawstring,
  svg_dir_set : bool,
  csv_dir : rawstring,
  csv_dir_set : bool,
  parallelism : uint,
  leaf_size : uint
  fixed_partition_size : uint,
  max_depth : uint,
}

terra parse_input_args()
  var conf : Config
  conf.time_steps = 10
  conf.input_file_set = false
  conf.num_bodies = -1
  conf.csv_dir_set = false
  conf.svg_dir_set = false
  conf.leaf_size = 32
  conf.parallelism = 16
  conf.fixed_partition_size = -1
  conf.max_depth = -1

  var args = c.legion_runtime_get_input_args()

  for i = 0, args.argc do
    if cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      conf.input_file = args.argv[i]
      conf.input_file_set = true
    elseif cstring.strcmp(args.argv[i], "-n") == 0 then
      i = i + 1
      conf.num_bodies = std.atoi(args.argv[i])
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

  return conf
end
