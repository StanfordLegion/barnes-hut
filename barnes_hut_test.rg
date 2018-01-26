import "regent"
require("barnes_hut_common")

local assert = regentlib.assert

local epsilon = 0.00001

terra double_equal(a: double, b: double)
	return a - b < epsilon and a - b > -epsilon
end

task main()
  var quads = region(ispace(ptr, 10), quad(quads))
  fill(quads.{nw, sw, ne, se}, null(ptr(quad(quads), quads)))
  var root = dynamic_cast(ptr(quad(quads), quads), 0)
  root.center_x = 1
  root.center_y = 1
  root.size = 2
  root.type = 2

  var body_quad = dynamic_cast(ptr(quad(quads), quads), 1)
  body_quad.mass_x = 0.5
  body_quad.mass_y = 1.5
  body_quad.mass = 1
  body_quad.type = 1
  var index = add_node(quads, root, body_quad, 1)
  assert(index == 1, "first insert failed")

  body_quad = dynamic_cast(ptr(quad(quads), quads), index + 1)
  body_quad.mass_x = 0.3
  body_quad.mass_y = 1.7
  body_quad.mass = 1
  body_quad.type = 1
  index = add_node(quads, root, body_quad, index + 1)
  assert(index == 3, "second insert failed")

  body_quad = dynamic_cast(ptr(quad(quads), quads), index + 1)
  body_quad.mass_x = 1.8
  body_quad.mass_y = 0.8
  body_quad.mass = 2
  body_quad.type = 1
  index = add_node(quads, root, body_quad, index + 1)
  assert(index == 4, "third insert failed")

  assert(root.total == 3, "root total failed")
  assert(root.mass == 4, "root mass failed")
  assert(double_equal(root.mass_x, 1.1), "root mass_x failed")
  assert(double_equal(root.mass_y, 1.2), "root mass_y failed")
end
regentlib.start(main)