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
  var index = add_node(quads, 0, 1, 1)
  assert(index == 1, "first insert failed")

  quads[index + 1].mass_x = 0.3
  quads[index + 1].mass_y = 1.7
  quads[index + 1].mass = 1
  quads[index + 1].type = 1
  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 3, "second insert failed")

  quads[index + 1].mass_x = 1.8
  quads[index + 1].mass_y = 0.8
  quads[index + 1].mass = 2
  quads[index + 1].type = 1
  index = add_node(quads, 0, index + 1, index + 1)
  assert(index == 4, "third insert failed")

  assert(quads[0].total == 3, "root total failed")
  regentlib.c.printf("mass %f %f %f", quads[0].mass, quads[0].mass_x, quads[0].mass_y)
  assert(quads[0].mass == 4, "root mass failed")
  assert(double_equal(quads[0].mass_x, 1.1), "root mass_x failed")
  assert(double_equal(quads[0].mass_y, 1.2), "root mass_y failed")
end

task main()
  test_quad_sizer()
  test_quad_tree()
end
regentlib.start(main)