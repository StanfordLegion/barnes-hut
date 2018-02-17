import "regent"
require("quad_tree")

local assert = regentlib.assert
local QuadTreeSizer = require("quad_tree_sizer")

local epsilon = 0.00001

terra double_equal(a: double, b: double)
  return a - b < epsilon and a - b > -epsilon
end

task test_quad_sizer()
  var root = create_placeholder()
  root.center_x = 1
  root.center_y = 1
  root.size = 2
  root.type = 2

  var body_quad = create_placeholder()
  body_quad.mass_x = 0.5
  body_quad.mass_y = 1.5
  body_quad.type = 1

  add_placeholder(root, body_quad)
  assert(count(root, false) == 2, "first insert failed")

  body_quad = create_placeholder()
  body_quad.mass_x = 0.3
  body_quad.mass_y = 1.7
  body_quad.type = 1

  add_placeholder(root, body_quad)
  assert(count(root, false) == 4, "second insert failed")

  body_quad = create_placeholder()
  body_quad.mass_x = 1.8
  body_quad.mass_y = 0.8
  body_quad.type = 1

  add_placeholder(root, body_quad)
  assert(count(root, false) == 5, "third insert failed")

  body_quad = create_placeholder()
  body_quad.mass_x = 0.2
  body_quad.mass_y = 1.8
  body_quad.type = 1

  add_placeholder(root, body_quad)
  assert(count(root, true) == 7, "fourth insert failed")
end

task test_quad_tree()
  var quads = region(ispace(int1d, 10), quad)
  fill(quads.{nw, sw, ne, se}, -1)
  fill(quads.type, 0)

  quads[0].center_x = 1
  quads[0].center_y = 1
  quads[0].size = 2
  quads[0].type = 2

  quads[1].mass_x = 0.5
  quads[1].mass_y = 1.5
  quads[1].mass = 1
  quads[1].type = 1
  quads[1].total = 1
  var index = add_node(quads, 0, 1, 1)
  assert(index == 1, "first insert failed")

  quads[index + 1].mass_x = 0.3
  quads[index + 1].mass_y = 1.7
  quads[index + 1].mass = 1
  quads[index + 1].type = 1
  quads[index + 1].total = 1
  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 3, "second insert failed")

  quads[index + 1].mass_x = 1.8
  quads[index + 1].mass_y = 0.8
  quads[index + 1].mass = 2
  quads[index + 1].type = 1
  quads[index + 1].total = 1
  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 4, "third insert failed")

  quads[index + 1].mass_x = 0.2
  quads[index + 1].mass_y = 1.8
  quads[index + 1].mass = 2
  quads[index + 1].type = 1
  quads[index + 1].total = 1
  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 6, "fourth insert failed")

  regentlib.c.printf("total %d mass %f %f %f", quads[0].total, quads[0].mass, quads[0].mass_x, quads[0].mass_y)
  assert(quads[0].total == 4, "root total failed")
  assert(quads[0].mass == 6, "root mass failed")
  assert(double_equal(quads[0].mass_x, 0.8), "root mass_x failed")
  assert(double_equal(quads[0].mass_y, 1.4), "root mass_y failed")
end

-- 32: x: -1800.000000, y: -1200.000000, speed_x: 0.000000, speed_y: 0.000000, mass: 224.000000
-- 33: x: -2017.915453, y: -1316.219045, speed_x: -5.538046, speed_y: 10.384067, mass: 1.272257
-- 34: x: -1642.261817, y: -1433.410078, speed_x: -9.863039, speed_y: -6.665425, mass: 1.033917

task test_divergence()
  var quads = region(ispace(int1d, 10), quad)
  fill(quads.{nw, sw, ne, se}, -1)
  fill(quads.type, 0)

  quads[0].center_x = -2151.000000
  quads[0].center_y = -1271.616941
  quads[0].size = 600.941984
  quads[0].type = 2

  quads[1].mass_x = -1800.000000
  quads[1].mass_y = -1200.000000
  quads[1].mass = 224.000000
  quads[1].type = 1
  quads[1].total = 1

  var root = create_placeholder()
  root.center_x = quads[0].center_x
  root.center_y = quads[0].center_y
  root.size = quads[0].size
  root.type = 2

  var body_quad = create_placeholder()
  body_quad.mass_x = quads[1].mass_x
  body_quad.mass_y = quads[1].mass_y
  body_quad.type = 1

  var index = add_node(quads, 0, 1, 1)
  assert(index == 1, "first insert failed")
  assert(quads[0].ne == 1, "first insert ne failed")

  add_placeholder(root, body_quad)
  assert(count(root, false) == 2, "first insert failed")
  assert(root.ne == body_quad, "first insert ne failed")

  quads[index + 1].mass_x = -2017.915453
  quads[index + 1].mass_y = -1316.219045
  quads[index + 1].mass = 1.272257
  quads[index + 1].type = 1
  quads[index + 1].total = 1

  body_quad = create_placeholder()
  body_quad.mass_x = quads[index + 1].mass_x
  body_quad.mass_y = quads[index + 1].mass_y
  body_quad.type = 1

  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 2, "second insert failed")
  assert(quads[0].se == 2, "second insert se failed")

  add_placeholder(root, body_quad)
  assert(count(root, false) == 3, "second insert failed")
  assert(root.se == body_quad, "second insert se failed")

  quads[index + 1].mass_x = -1642.261817
  quads[index + 1].mass_y = -1433.410078
  quads[index + 1].mass = 1.033917
  quads[index + 1].type = 1
  quads[index + 1].total = 1

  var old_body_quad = body_quad
  body_quad = create_placeholder()
  body_quad.mass_x = quads[index + 1].mass_x
  body_quad.mass_y = quads[index + 1].mass_y
  body_quad.type = 1

  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 4, "third insert failed")
  assert(quads[quads[0].se].nw == 2, "third insert se.nw failed")
  assert(quads[quads[0].se].se == 3, "third insert se.se failed")

  add_placeholder(root, body_quad)
  assert(root.se.nw == old_body_quad, "second insert se.nw failed")
  assert(root.se.se == body_quad, "second insert se.se failed")
  assert(count(root, true) == 5, "second insert failed")
end

task main()
  test_quad_sizer()
  test_quad_tree()
  test_divergence()
end
regentlib.start(main)