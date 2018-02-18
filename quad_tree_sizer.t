import "regent"

local c = regentlib.c

struct quad_placeholder {
  mass_x : double,
  mass_y : double,
  center_x : double,
  center_y : double,
  size : double,
  type : uint,
  ne : &quad_placeholder,
  nw : &quad_placeholder,
  se : &quad_placeholder,
  sw : &quad_placeholder,
  next_in_leaf : &quad_placeholder,
  leaf_count : uint,
}

struct quad_chunk {
  prev : &quad_chunk,
  placeholders : &quad_placeholder
  size : uint,
  counter : uint,
}

terra create_quad_chunk(size : uint): &quad_chunk
  var chunk : &quad_chunk = [&quad_chunk](c.malloc(sizeof(quad_chunk)))
  chunk.placeholders = [&quad_placeholder](c.malloc(size * sizeof(quad_placeholder)))
  chunk.counter = 0
  chunk.size = size
  chunk.prev = nil
  return chunk
end

terra init_placeholder(chunk : &quad_chunk): &quad_placeholder
  if chunk.counter == chunk.size then
    var new_chunk = create_quad_chunk(chunk.size)

    var temp_placeholders = chunk.placeholders
    chunk.placeholders = new_chunk.placeholders
    new_chunk.placeholders = temp_placeholders

    new_chunk.prev = chunk.prev
    chunk.prev = new_chunk

    chunk.counter = 0
    new_chunk.counter = chunk.size
  end

  chunk.placeholders[chunk.counter].nw = nil
  chunk.placeholders[chunk.counter].sw = nil
  chunk.placeholders[chunk.counter].ne = nil
  chunk.placeholders[chunk.counter].se = nil
  chunk.placeholders[chunk.counter].next_in_leaf = nil

  chunk.counter = chunk.counter + 1
  
  return &chunk.placeholders[chunk.counter - 1]
end

terra count(chunk : &quad_chunk, free : bool): uint
  var total = 0

  var current = chunk
  while current ~= nil do
    total = total + current.counter
    var old_current = current
    current = current.prev

    if free then
      c.free(old_current.placeholders)
      c.free(old_current)
    end
  end

  return total
end

terra insert_in_leaf(child : &quad_placeholder, old_child : &quad_placeholder)
  child.leaf_count = old_child.leaf_count + 1
  child.next_in_leaf = old_child
end

terra insert_leaf_into_fork(new_fork : &quad_placeholder, leaf : &quad_placeholder, chunk : &quad_chunk, leaf_size : uint)
  var current = leaf
  while current ~= nil do
    var next_in_leaf = current.next_in_leaf
    current.next_in_leaf = nil
    add_placeholder(new_fork, current, chunk, leaf_size)
    current = next_in_leaf
  end
end

terra add_placeholder(parent: &quad_placeholder, child: &quad_placeholder, chunk : &quad_chunk, leaf_size : uint): uint
  var half_size = parent.size / 2
  if child.mass_x <= parent.center_x then
    if child.mass_y <= parent.center_y then
      if parent.sw == nil then
        child.leaf_count = 1
        parent.sw = child
      elseif parent.sw.type == 1 then
        if parent.sw.leaf_count < leaf_size then
          insert_in_leaf(child, parent.sw)
          parent.sw = child
        else
          var new_fork = init_placeholder(chunk)
          new_fork.type = 2
          new_fork.center_x = parent.center_x - half_size / 2
          new_fork.center_y = parent.center_y - half_size / 2
          new_fork.size = half_size

          insert_leaf_into_fork(new_fork, parent.sw, chunk, leaf_size)
          add_placeholder(new_fork, child, chunk, leaf_size)
          parent.sw = new_fork
        end
      else
        add_placeholder(parent.sw, child, chunk, leaf_size)
      end
    else
      if parent.nw == nil then
        child.leaf_count = 1
        parent.nw = child
      elseif parent.nw.type == 1 then
        if parent.nw.leaf_count < leaf_size then
          insert_in_leaf(child, parent.nw)
          parent.nw = child
        else
          var new_fork = init_placeholder(chunk)
          new_fork.type = 2
          new_fork.center_x = parent.center_x - half_size / 2
          new_fork.center_y = parent.center_y + half_size / 2
          new_fork.size = half_size

          insert_leaf_into_fork(new_fork, parent.nw, chunk, leaf_size)
          add_placeholder(new_fork, child, chunk, leaf_size)
          parent.nw = new_fork
        end
      else
        add_placeholder(parent.nw, child, chunk, leaf_size)
      end      
    end
  else
    if child.mass_y <= parent.center_y then
      if parent.se == nil then
        child.leaf_count = 1
        parent.se = child
      elseif parent.se.type == 1 then
        if parent.se.leaf_count < leaf_size then
          insert_in_leaf(child, parent.se)
          parent.se = child
        else
          var new_fork = init_placeholder(chunk)
          new_fork.type = 2
          new_fork.center_x = parent.center_x + half_size / 2
          new_fork.center_y = parent.center_y - half_size / 2
          new_fork.size = half_size

          insert_leaf_into_fork(new_fork, parent.se, chunk, leaf_size)
          add_placeholder(new_fork, child, chunk, leaf_size)
          parent.se = new_fork
        end
      else
        add_placeholder(parent.se, child, chunk, leaf_size)
      end
    else
      if parent.ne == nil then
        child.leaf_count = 1
        parent.ne = child
      elseif parent.ne.type == 1 then
        if parent.ne.leaf_count < leaf_size then
          insert_in_leaf(child, parent.ne)
          parent.ne = child
        else
          var new_fork = init_placeholder(chunk)
          new_fork.type = 2
          new_fork.center_x = parent.center_x + half_size / 2
          new_fork.center_y = parent.center_y + half_size / 2
          new_fork.size = half_size

          insert_leaf_into_fork(new_fork, parent.ne, chunk, leaf_size)
          add_placeholder(new_fork, child, chunk, leaf_size)
          parent.ne = new_fork
        end
      else
        add_placeholder(parent.ne, child, chunk, leaf_size)
      end      
    end
  end

  return 1
end
