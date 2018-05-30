import "regent"

fspace body {
  {mass, mass_x, mass_y, speed_x, speed_y} : double,
  {sector, eliminated} : int1d,
  index : uint,
}
