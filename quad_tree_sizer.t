import "regent"

struct quad_placeholder {
  mass_x: double,
  mass_y: double,
  center_x: double,
  center_y: double,
  size: double,
  type: uint,
  num_elements: uint
  ne: &quad_placeholder,
  nw: &quad_placeholder,
  se: &quad_placeholder,
  sw: &quad_placeholder
}

terra create_placeholder(): quad_placeholder
  var placeholder : quad_placeholder
  placeholder.nw = nil
  placeholder.sw = nil
  placeholder.ne = nil
  placeholder.se = nil
  placeholder.num_elements = 1
  return placeholder
end

terra add_placeholder(parent: quad_placeholder, child: quad_placeholder): quad_placeholder
  var half_size = parent.size / 2
  if child.mass_x <= parent.center_x then
    if child.mass_y <= parent.center_y then
      if parent.sw == nil then
        parent.sw = &child
        parent.num_elements = parent.num_elements + 1
      elseif parent.sw.type == 1 then
        var new_fork = create_placeholder()
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        new_fork = add_placeholder(new_fork, @parent.sw)
        new_fork = add_placeholder(new_fork, child)
        parent.sw = &new_fork
        parent.num_elements = parent.num_elements + new_fork.num_elements - 1 
      else
        var old_child = @parent.sw
        var before = old_child.num_elements
        var new_child = add_placeholder(old_child, child)
        parent.sw = &new_child
        parent.num_elements = parent.num_elements - before + new_child.num_elements
      end
    else
      if parent.nw == nil then
        parent.nw = &child
        parent.num_elements = parent.num_elements + 1
      elseif parent.nw.type == 1 then
        var new_fork = create_placeholder()
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y + half_size / 2
        new_fork.size = half_size
        new_fork = add_placeholder(new_fork, @parent.nw)
        new_fork = add_placeholder(new_fork, child)
        parent.nw = &new_fork
        parent.num_elements = parent.num_elements + new_fork.num_elements - 1 
      else
        var old_child = @parent.nw
        var before = old_child.num_elements
        var new_child = add_placeholder(old_child, child)
        parent.nw = &new_child
        parent.num_elements = parent.num_elements - before + new_child.num_elements
      end      
    end
  else
    if child.mass_y <= parent.center_y then
      if parent.se == nil then
        parent.se = &child
        parent.num_elements = parent.num_elements + 1
      elseif parent.se.type == 1 then
        var new_fork = create_placeholder()
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        new_fork = add_placeholder(new_fork, @parent.se)
        new_fork = add_placeholder(new_fork, child)
        parent.se = &new_fork
        parent.num_elements = parent.num_elements + new_fork.num_elements - 1 
      else
        var old_child = @parent.se
        var before = old_child.num_elements
        var new_child = add_placeholder(old_child, child)
        parent.se = &new_child
        parent.num_elements = parent.num_elements - before + new_child.num_elements
      end
    else
      if parent.ne == nil then
        parent.ne = &child
        parent.num_elements = parent.num_elements + 1
      elseif parent.ne.type == 1 then
        var new_fork = create_placeholder()
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_x + half_size / 2
        new_fork.size = half_size
        new_fork = add_placeholder(new_fork, @parent.ne)
        new_fork = add_placeholder(new_fork, child)
        parent.ne = &new_fork
        parent.num_elements = parent.num_elements + new_fork.num_elements - 1 
      else
        var old_child = @parent.ne
        var before = old_child.num_elements
        var new_child = add_placeholder(old_child, child)
        parent.ne = &new_child
        parent.num_elements = parent.num_elements - before + new_child.num_elements
      end      
    end
  end

  return parent
end
