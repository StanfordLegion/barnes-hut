import "regent"
require("barnes_hut_common")

local assert = regentlib.assert

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
  body_quad.type = 1
  var index = add_node(quads, root, body_quad, 1)
  assert(index == 1, "test failed")

  body_quad = dynamic_cast(ptr(quad(quads), quads), index + 1)
  body_quad.mass_x = 0.3
  body_quad.mass_y = 1.7
  index = add_node(quads, root, body_quad, index + 1)
  assert(index == 3, "test failed")
end
regentlib.start(main)