import "regent"

fspace body {
  {mass, mass_x, mass_y, speed_x, speed_y} : double,
  {sector, eliminated} : int1d,
  index : uint,
}

fspace boundary {
  {min_x, min_y, max_x, max_y} : double
}
