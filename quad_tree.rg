import "regent"

local assert = regentlib.assert

fspace quad {
  {mass_x, mass_y, mass, center_x, center_y, size} : double,
  {total, leaf_count, type, index} : uint,
  {ne, nw, se, sw, next_in_leaf} : int,
}

task add_node(quads : region(ispace(int1d), quad), parent_index : int, child_index : int, leaf_size : uint, last_used : int): int
where
  reads(quads),
  writes(quads)
do
  -- regentlib.c.printf("entered parent: %d, child: %d, %d\n", parent_index, child_index, last_used)
  var half_size = quads[parent_index].size / 2
  var new_last_used = last_used
  if quads[child_index].mass_x <= quads[parent_index].center_x then
    if quads[child_index].mass_y <= quads[parent_index].center_y then
      if quads[parent_index].sw == -1 then
        -- regentlib.c.printf("new sw\n")
        quads[child_index].leaf_count = 1
        quads[parent_index].sw = child_index
      elseif quads[quads[parent_index].sw].type == 1 then
        if quads[quads[parent_index].sw].leaf_count < leaf_size then
          quads[child_index].leaf_count = quads[quads[parent_index].sw].leaf_count + 1
          quads[child_index].next_in_leaf = quads[parent_index].sw
          quads[parent_index].sw = child_index
        else
          -- regentlib.c.printf("allocating %d\n", last_used + 1)
          assert(quads[last_used + 1].type == 0, "region already allocated")
          quads[last_used + 1].type = 2
          quads[last_used + 1].center_x = quads[parent_index].center_x - half_size / 2
          quads[last_used + 1].center_y = quads[parent_index].center_y - half_size / 2
          quads[last_used + 1].size = half_size
          
          var current = quads[parent_index].sw
          while current ~= -1 do
            var next_in_leaf = quads[current].next_in_leaf
            quads[current].next_in_leaf = -1
            new_last_used = add_node(quads, last_used + 1, current, leaf_size, last_used + 1)
            current = next_in_leaf
          end

          new_last_used = add_node(quads, last_used + 1, child_index, leaf_size, new_last_used)
          
          quads[parent_index].sw = last_used + 1
        end
      else
        -- regentlib.c.printf("other %d %d\n", quads[parent_index].sw, quads[quads[parent_index].sw].type)
        new_last_used = add_node(quads, quads[parent_index].sw, child_index, leaf_size, last_used)
      end
    else
      if quads[parent_index].nw == -1 then
        -- regentlib.c.printf("new nw\n")
        quads[child_index].leaf_count = 1
        quads[parent_index].nw = child_index
      elseif quads[quads[parent_index].nw].type == 1 then
        if quads[quads[parent_index].nw].leaf_count < leaf_size then
          quads[child_index].leaf_count = quads[quads[parent_index].nw].leaf_count + 1
          quads[child_index].next_in_leaf = quads[parent_index].nw
          quads[parent_index].nw = child_index
        else
          -- regentlib.c.printf("allocating %d\n", last_used + 1)
          assert(quads[last_used + 1].type == 0, "region already allocated")
          quads[last_used + 1].type = 2
          quads[last_used + 1].center_x = quads[parent_index].center_x - half_size / 2
          quads[last_used + 1].center_y = quads[parent_index].center_y + half_size / 2
          quads[last_used + 1].size = half_size

          var current = quads[parent_index].nw
          while current ~= -1 do
            var next_in_leaf = quads[current].next_in_leaf
            quads[current].next_in_leaf = -1
            new_last_used = add_node(quads, last_used + 1, current, leaf_size, last_used + 1)
            current = next_in_leaf
          end

          new_last_used = add_node(quads, last_used + 1, child_index, leaf_size, new_last_used)
          quads[parent_index].nw = last_used + 1
        end
      else
        -- regentlib.c.printf("other %d %d\n", quads[parent_index].nw, quads[quads[parent_index].nw].type)
        new_last_used = add_node(quads, quads[parent_index].nw, child_index, leaf_size, last_used)
      end      
    end
  else
    if quads[child_index].mass_y <= quads[parent_index].center_y then
      if quads[parent_index].se == -1 then
        -- regentlib.c.printf("new se\n")
        quads[child_index].leaf_count = 1
        quads[parent_index].se = child_index
      elseif quads[quads[parent_index].se].type == 1 then
        if quads[quads[parent_index].se].leaf_count < leaf_size then
          quads[child_index].leaf_count = quads[quads[parent_index].se].leaf_count + 1
          quads[child_index].next_in_leaf = quads[parent_index].se
          quads[parent_index].se = child_index
        else
          -- regentlib.c.printf("allocating %d\n", last_used + 1)
          assert(quads[last_used + 1].type == 0, "region already allocated")
          quads[last_used + 1].type = 2
          quads[last_used + 1].center_x = quads[parent_index].center_x + half_size / 2
          quads[last_used + 1].center_y = quads[parent_index].center_y - half_size / 2
          quads[last_used + 1].size = half_size

          var current = quads[parent_index].se
          while current ~= -1 do
            var next_in_leaf = quads[current].next_in_leaf
            quads[current].next_in_leaf = -1
            new_last_used = add_node(quads, last_used + 1, current, leaf_size, last_used + 1)
            current = next_in_leaf
          end

          new_last_used = add_node(quads, last_used + 1, child_index, leaf_size, new_last_used)
          quads[parent_index].se = last_used + 1
        end
      else
        -- regentlib.c.printf("other %d %d\n", quads[parent_index].se, quads[quads[parent_index].se].type)
        new_last_used = add_node(quads, quads[parent_index].se, child_index, leaf_size, last_used)
      end
    else
      if quads[parent_index].ne == -1 then
        -- regentlib.c.printf("new ne\n")
        quads[child_index].leaf_count = 1
        quads[parent_index].ne = child_index
      elseif quads[quads[parent_index].ne].type == 1 then
        if quads[quads[parent_index].ne].leaf_count < leaf_size then
          quads[child_index].leaf_count = quads[quads[parent_index].ne].leaf_count + 1
          quads[child_index].next_in_leaf = quads[parent_index].ne
          quads[parent_index].ne = child_index
        else
          -- regentlib.c.printf("allocating %d\n", last_used + 1)
          assert(quads[last_used + 1].type == 0, "region already allocated")
          quads[last_used + 1].type = 2
          quads[last_used + 1].center_x = quads[parent_index].center_x + half_size / 2
          quads[last_used + 1].center_y = quads[parent_index].center_y + half_size / 2
          quads[last_used + 1].size = half_size

          var current = quads[parent_index].ne
          while current ~= -1 do
            var next_in_leaf = quads[current].next_in_leaf
            quads[current].next_in_leaf = -1
            new_last_used = add_node(quads, last_used + 1, current, leaf_size, last_used + 1)
            current = next_in_leaf
          end

          new_last_used = add_node(quads, last_used + 1, child_index, leaf_size, new_last_used)
          quads[parent_index].ne = last_used + 1
        end
      else
        -- regentlib.c.printf("other %d %d\n", quads[parent_index].ne, quads[quads[parent_index].ne].type)
        new_last_used = add_node(quads, quads[parent_index].ne, child_index, leaf_size, last_used)
      end      
    end
  end

  quads[parent_index].mass = 0
  quads[parent_index].mass_x = 0
  quads[parent_index].mass_y = 0
  
  if quads[parent_index].sw ~= -1 then
    -- regentlib.c.printf("Updating mass sw\n")
    var mass = quads[quads[parent_index].sw].mass
    quads[parent_index].mass += mass
    quads[parent_index].mass_x += quads[quads[parent_index].sw].mass_x * mass
    quads[parent_index].mass_y += quads[quads[parent_index].sw].mass_y * mass
  end 

  if quads[parent_index].nw ~= -1 then
    -- regentlib.c.printf("Updating mass nw\n")
    var mass = quads[quads[parent_index].nw].mass
    quads[parent_index].mass += mass
    quads[parent_index].mass_x += quads[quads[parent_index].nw].mass_x * mass
    quads[parent_index].mass_y += quads[quads[parent_index].nw].mass_y * mass
  end

  if quads[parent_index].se ~= -1 then
    -- regentlib.c.printf("Updating mass se\n")
    var mass = quads[quads[parent_index].se].mass
    quads[parent_index].mass += mass
    quads[parent_index].mass_x += quads[quads[parent_index].se].mass_x * mass
    quads[parent_index].mass_y += quads[quads[parent_index].se].mass_y * mass
  end

  if quads[parent_index].ne ~= -1 then
    -- regentlib.c.printf("Updating mass ne\n")
    var mass = quads[quads[parent_index].ne].mass
    quads[parent_index].mass += mass
    quads[parent_index].mass_x += quads[quads[parent_index].ne].mass_x * mass
    quads[parent_index].mass_y += quads[quads[parent_index].ne].mass_y * mass
  end

  if quads[parent_index].mass > 0 then
    quads[parent_index].mass_x = quads[parent_index].mass_x / quads[parent_index].mass
    quads[parent_index].mass_y = quads[parent_index].mass_y / quads[parent_index].mass
  end

  quads[parent_index].total += quads[child_index].total

  -- regentlib.c.printf("returning %d\n", new_last_used)
  return new_last_used
end
