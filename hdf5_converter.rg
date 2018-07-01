import "regent"

require("barnes_hut_common")

local c = regentlib.c
local cstring = terralib.includec("string.h")
local std = terralib.includec("stdlib.h")

local hdf5 = terralib.includec("hdf5.h")
hdf5.H5F_ACC_TRUNC = 2
hdf5.H5T_STD_I32LE = hdf5.H5T_STD_I32LE_g
hdf5.H5T_IEEE_F64LE = hdf5.H5T_IEEE_F64LE_g
hdf5.H5P_DEFAULT = 0

struct Config {
  input_file : rawstring,
  input_file_set : bool,
  output_file : rawstring,
  output_file_set : bool,
  num_bodies : uint,
}

terra parse_input_args()
  var conf : Config
  conf.input_file_set = false
  conf.output_file_set = false
  conf.num_bodies = -1

  var args = c.legion_runtime_get_input_args()

  for i = 0, args.argc do
    if cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1
      conf.input_file = args.argv[i]
      conf.input_file_set = true
    elseif cstring.strcmp(args.argv[i], "-o") == 0 then
      i = i + 1
      conf.output_file = args.argv[i]
      conf.output_file_set = true
    end
  end

  return conf
end

terra generate_hdf5_file(filename : rawstring, num_bodies: uint)
  var fid = hdf5.H5Fcreate(filename, hdf5.H5F_ACC_TRUNC, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)

  var h_dims : hdf5.hsize_t[1]
  h_dims[0] = num_bodies
  var did = hdf5.H5Screate_simple(1, h_dims, [&uint64](0))

  var ds1id = hdf5.H5Dcreate2(fid, "index", hdf5.H5T_STD_I32LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds1id)

  var ds2id = hdf5.H5Dcreate2(fid, "mass_x", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds2id)

  var ds3id = hdf5.H5Dcreate2(fid, "mass_y", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds3id)

  var ds4id = hdf5.H5Dcreate2(fid, "speed_x", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds4id)

  var ds5id = hdf5.H5Dcreate2(fid, "speed_y", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds5id)

  var ds6id = hdf5.H5Dcreate2(fid, "mass", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds6id)

  hdf5.H5Sclose(did)
  hdf5.H5Fclose(fid)
end

task main()
  var conf = parse_input_args()

  var input = region(ispace(ptr, conf.num_bodies), body)

  var fp = c.fopen(conf.input_file, "r")
  var line : int8[1024]

  for i=0,conf.num_bodies do
    c.fgets(line, 1024, fp)
    var index = std.atoi(cstring.strtok(line, ","))
    input[index].index = index
    input[index].mass_x = std.atof(cstring.strtok([&int8](0), ","))
    input[index].mass_y = std.atof(cstring.strtok([&int8](0), ","))
    input[index].speed_x = std.atof(cstring.strtok([&int8](0), ","))
    input[index].speed_y = std.atof(cstring.strtok([&int8](0), ","))
    input[index].mass = std.atof(cstring.strtok([&int8](0), ","))
  end

  c.fclose(fp)

  generate_hdf5_file(conf.output_file, conf.num_bodies)

  var output = region(ispace(ptr, conf.num_bodies), body)
  attach(hdf5, output.{mass, mass_x, mass_y, speed_x, speed_y, index}, conf.output_file, regentlib.file_read_write)
  acquire(output.{mass, mass_x, mass_y, speed_x, speed_y, index})
  copy(input.{mass, mass_x, mass_y, speed_x, speed_y, index}, output.{mass, mass_x, mass_y, speed_x, speed_y, index})
  release(output.{mass, mass_x, mass_y, speed_x, speed_y, index})
  detach(hdf5, output.{mass, mass_x, mass_y, speed_x, speed_y, index})
end

regentlib.start(main)
