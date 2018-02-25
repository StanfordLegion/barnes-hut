import "regent"

fspace body {
  {mass, mass_x, mass_y, speed_x, speed_y, force_x, force_y} : double,
  sector : int1d,
  index : uint,
}
