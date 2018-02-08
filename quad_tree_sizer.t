import "regent"

struct quad_placeholder {
  mass_x: double,
  mass_y: double,
  center_x: double,
  center_y: double,
  size: double,
  type: uint,
  ne: &quad_placeholder,
  nw: &quad_placeholder,
  se: &quad_placeholder,
  sw: &quad_placeholder
}

terra add_placeholder(parent: quad_placeholder, child: quad_placeholder): uint
  var half_size = parent.size / 2
  var added = 0
  if child.mass_x <= parent.center_x then
    if child.mass_y <= parent.center_y then
      if parent.sw == nil then
        parent.sw = &child
      elseif parent.sw.type == 1 then
        var new_fork : quad_placeholder
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        added = added + 1
        added = added + add_placeholder(new_fork, @parent.sw)
        added = added + add_placeholder(new_fork, child)
        parent.sw = &new_fork
      else
        added = added + add_placeholder(@parent.sw, child)
      end
    else
      if parent.nw == nil then
        parent.nw = &child
      elseif parent.nw.type == 1 then
        var new_fork : quad_placeholder
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y + half_size / 2
        new_fork.size = half_size
        added = added + 1
        added = added + add_placeholder(new_fork, @parent.nw)
        added = added + add_placeholder(new_fork, child)
        parent.nw = &new_fork
      else
        added = added + add_placeholder(@parent.nw, child)
      end      
    end
  else
    if child.mass_y <= parent.center_y then
      if parent.se == nil then
        parent.se = &child
      elseif parent.se.type == 1 then
        var new_fork : quad_placeholder
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        added = added + 1
        added = added + add_placeholder(new_fork, @parent.se)
        added = added + add_placeholder(new_fork, child)
        parent.se = &new_fork
      else
        added = added + add_placeholder(@parent.se, child)
      end
    else
      if parent.ne == nil then
        parent.ne = &child
      elseif parent.ne.type == 1 then
        var new_fork : quad_placeholder
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_x + half_size / 2
        new_fork.size = half_size
        added = added + 1
        added = added + add_placeholder(new_fork, @parent.ne)
        added = added + add_placeholder(new_fork, child)
        parent.ne = &new_fork
      else
        added = added + add_placeholder(@parent.ne, child)
      end      
    end
  end

  return added
end
