local c = regentlib.c
local floor = regentlib.floor(float)
local log10 = regentlib.log10(float)

local cstring = terralib.includec("string.h")

terra open(iteration: uint, dirname: rawstring)
  var svg_path : int8[1000]
  c.sprintf(svg_path, "%s/%d.svg", dirname, iteration)
  return c.fopen(svg_path, "w")
end