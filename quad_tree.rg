import "regent"

local assert = regentlib.assert

fspace quad {
  {mass_x, mass_y, mass, center_x, center_y, size} : double,
  {total, type, index} : uint,
  {ne, nw, se, sw} : int,
}

task add_node(quads : region(ispace(int1d), quad), parent_index : int, child_index : int, last_used : int): int
where
  reads(quads),
  writes(quads)
do
  -- regentlib.c.printf("entered parent: %d, child: %d, %d\n", parent_index, child_index, last_used)
  var parent = quads[parent_index]
  var child = quads[child_index]

  var half_size = parent.size / 2
  var new_last_used = last_used
  if child.mass_x <= parent.center_x then
    if child.mass_y <= parent.center_y then
      if parent.sw == -1 then
        -- regentlib.c.printf("new sw")
        parent.sw = child_index
      elseif quads[parent.sw].type == 1 then
        var new_fork = quads[last_used + 1]
        -- regentlib.c.printf("allocating %d\n", last_used + 1)
        assert(new_fork.type == 0, "region already allocated")
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        new_last_used = add_node(quads, last_used + 1, parent.sw, last_used + 1)
        new_last_used = add_node(quads, last_used + 1, child_index, new_last_used)
        parent.sw = last_used + 1
      else
        -- regentlib.c.printf("other %d %d\n", parent.sw, quads[parent.sw].type)
        new_last_used = add_node(quads, parent.sw, child_index, last_used)
      end
    else
      if parent.nw == -1 then
        -- regentlib.c.printf("new nw")
        parent.nw = child_index
      elseif quads[parent.nw].type == 1 then
        var new_fork = quads[last_used + 1]
        -- regentlib.c.printf("allocating %d\n", last_used + 1)
        assert(new_fork.type == 0, "region already allocated")
        new_fork.type = 2
        new_fork.center_x = parent.center_x - half_size / 2
        new_fork.center_y = parent.center_y + half_size / 2
        new_fork.size = half_size
        new_last_used = add_node(quads, last_used + 1, parent.nw, last_used + 1)
        new_last_used = add_node(quads, last_used + 1, child_index, new_last_used)
        parent.nw = last_used + 1
      else
        -- regentlib.c.printf("other %d %d\n", parent.nw, quads[parent.nw].type)
        new_last_used = add_node(quads, parent.nw, child_index, last_used)
      end      
    end
  else
    if child.mass_y <= parent.center_y then
      if parent.se == -1 then
        -- regentlib.c.printf("new se")
        parent.se = child_index
      elseif quads[parent.se].type == 1 then
        var new_fork = quads[last_used + 1]
        -- regentlib.c.printf("allocating %d\n", last_used + 1)
        assert(new_fork.type == 0, "region already allocated")
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_y - half_size / 2
        new_fork.size = half_size
        new_last_used = add_node(quads, last_used + 1, parent.se, last_used + 1)
        new_last_used = add_node(quads, last_used + 1, child_index, new_last_used)
        parent.se = last_used + 1
      else
        -- regentlib.c.printf("other %d %d\n", parent.se, quads[parent.se].type)
        new_last_used = add_node(quads, parent.se, child_index, last_used)
      end
    else
      if parent.ne == -1 then
        -- regentlib.c.printf("new ne")
        parent.ne = child_index
      elseif quads[parent.ne].type == 1 then
        var new_fork = quads[last_used + 1]
        -- regentlib.c.printf("allocating %d\n", last_used + 1)
        assert(new_fork.type == 0, "region already allocated")
        new_fork.type = 2
        new_fork.center_x = parent.center_x + half_size / 2
        new_fork.center_y = parent.center_y + half_size / 2
        new_fork.size = half_size
        new_last_used = add_node(quads, last_used + 1, parent.ne, last_used + 1)
        new_last_used = add_node(quads, last_used + 1, child_index, new_last_used)
        parent.ne = last_used + 1
      else
        -- regentlib.c.printf("other %d %d\n", parent.ne, quads[parent.ne].type)
        new_last_used = add_node(quads, parent.ne, child_index, last_used)
      end      
    end
  end

  parent.mass = 0
  parent.mass_x = 0
  parent.mass_y = 0
  
  if not parent.sw == -1 then
    var mass = quads[parent.sw].mass
    parent.mass += mass
    parent.mass_x += quads[parent.sw].mass_x * mass
    parent.mass_y += quads[parent.sw].mass_y * mass
  end 

  if not parent.nw == -1 then
    var mass = quads[parent.nw].mass
    parent.mass += mass
    parent.mass_x += quads[parent.nw].mass_x * mass
    parent.mass_y += quads[parent.nw].mass_y * mass
  end

  if not parent.se == -1 then
    var mass = quads[parent.se].mass
    parent.mass += mass
    parent.mass_x += quads[parent.se].mass_x * mass
    parent.mass_y += quads[parent.se].mass_y * mass
  end

  if not parent.ne == -1 then
    var mass = quads[parent.ne].mass
    parent.mass += mass
    parent.mass_x += quads[parent.ne].mass_x * mass
    parent.mass_y += quads[parent.ne].mass_y * mass
  end

  if parent.mass > 0 then
    parent.mass_x = parent.mass_x / parent.mass
    parent.mass_y = parent.mass_y / parent.mass
  end

  parent.total += 1

  -- regentlib.c.printf("returning %d\n", new_last_used)
  return new_last_used
end
