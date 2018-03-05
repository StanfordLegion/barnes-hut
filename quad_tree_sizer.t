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

terra add_placeholder(root: &quad_placeholder, body: &quad_placeholder, chunk : &quad_chunk, leaf_size : uint): uint
  var parent_list : (&quad_placeholder)[1024]
  var child_list : (&quad_placeholder)[1024]
  parent_list[0] = root
  child_list[0] = body
  var traverse_index = 0

  while traverse_index >= 0 do
    regentlib.assert(traverse_index < 1024 - leaf_size - 1, "possible overflow")
    var parent = parent_list[traverse_index]
    var child = child_list[traverse_index]
    traverse_index = traverse_index - 1

    var half_size = parent.size / 2
    if child.mass_x <= parent.center_x then
      if child.mass_y <= parent.center_y then
        if parent.sw == nil then
          child.leaf_count = 1
          parent.sw = child
        elseif parent.sw.type == 1 then
          if parent.sw.leaf_count < leaf_size or parent.size < 0.00001 then
            child.leaf_count = parent.sw.leaf_count + 1
            child.next_in_leaf = parent.sw
            parent.sw = child
          else
            var new_fork = init_placeholder(chunk)
            new_fork.type = 2
            new_fork.center_x = parent.center_x - half_size / 2
            new_fork.center_y = parent.center_y - half_size / 2
            new_fork.size = half_size

            var current = parent.sw
            while current ~= nil do
              var next_in_leaf = current.next_in_leaf
              current.next_in_leaf = nil
              traverse_index = traverse_index + 1
              parent_list[traverse_index] = new_fork
              child_list[traverse_index] = current
              current = next_in_leaf
            end

            traverse_index = traverse_index + 1
            parent_list[traverse_index] = new_fork
            child_list[traverse_index] = child

            parent.sw = new_fork
          end
        else
          traverse_index = traverse_index + 1
          parent_list[traverse_index] = parent.sw
          child_list[traverse_index] = child
        end
      else
        if parent.nw == nil then
          child.leaf_count = 1
          parent.nw = child
        elseif parent.nw.type == 1 then
          if parent.nw.leaf_count < leaf_size or parent.size < 0.00001 then
            child.leaf_count = parent.nw.leaf_count + 1
            child.next_in_leaf = parent.nw
            parent.nw = child
          else
            var new_fork = init_placeholder(chunk)
            new_fork.type = 2
            new_fork.center_x = parent.center_x - half_size / 2
            new_fork.center_y = parent.center_y + half_size / 2
            new_fork.size = half_size

            var current = parent.nw
            while current ~= nil do
              var next_in_leaf = current.next_in_leaf
              current.next_in_leaf = nil
              traverse_index = traverse_index + 1
              parent_list[traverse_index] = new_fork
              child_list[traverse_index] = current
              current = next_in_leaf
            end

            traverse_index = traverse_index + 1
            parent_list[traverse_index] = new_fork
            child_list[traverse_index] = child

            parent.nw = new_fork
          end
        else
          traverse_index = traverse_index + 1
          parent_list[traverse_index] = parent.nw
          child_list[traverse_index] = child
        end  
      end
    else
      if child.mass_y <= parent.center_y then
        if parent.se == nil then
          child.leaf_count = 1
          parent.se = child
        elseif parent.se.type == 1 then
          if parent.se.leaf_count < leaf_size or parent.size < 0.00001 then
            child.leaf_count = parent.se.leaf_count + 1
            child.next_in_leaf = parent.se
            parent.se = child
          else
            var new_fork = init_placeholder(chunk)
            new_fork.type = 2
            new_fork.center_x = parent.center_x + half_size / 2
            new_fork.center_y = parent.center_y - half_size / 2
            new_fork.size = half_size

            var current = parent.se
            while current ~= nil do
              var next_in_leaf = current.next_in_leaf
              current.next_in_leaf = nil
              traverse_index = traverse_index + 1
              parent_list[traverse_index] = new_fork
              child_list[traverse_index] = current
              current = next_in_leaf
            end

            traverse_index = traverse_index + 1
            parent_list[traverse_index] = new_fork
            child_list[traverse_index] = child

            parent.se = new_fork
          end
        else
          traverse_index = traverse_index + 1
          parent_list[traverse_index] = parent.se
          child_list[traverse_index] = child
        end
      else
        if parent.ne == nil then
          child.leaf_count = 1
          parent.ne = child
        elseif parent.ne.type == 1 then
          if parent.ne.leaf_count < leaf_size or parent.size < 0.00001 then
            child.leaf_count = parent.ne.leaf_count + 1
            child.next_in_leaf = parent.ne
            parent.ne = child
          else
            var new_fork = init_placeholder(chunk)
            new_fork.type = 2
            new_fork.center_x = parent.center_x + half_size / 2
            new_fork.center_y = parent.center_y + half_size / 2
            new_fork.size = half_size

            var current = parent.ne
            while current ~= nil do
              var next_in_leaf = current.next_in_leaf
              current.next_in_leaf = nil
              traverse_index = traverse_index + 1
              parent_list[traverse_index] = new_fork
              child_list[traverse_index] = current
              current = next_in_leaf
            end

            traverse_index = traverse_index + 1
            parent_list[traverse_index] = new_fork
            child_list[traverse_index] = child

            parent.ne = new_fork
          end
        else
          traverse_index = traverse_index + 1
          parent_list[traverse_index] = parent.ne
          child_list[traverse_index] = child
        end
      end
    end
  end

  return 1
end
