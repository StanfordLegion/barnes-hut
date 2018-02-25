import "regent"
require("barnes_hut_common")

local assert = regentlib.assert

local cmath = terralib.includec("math.h")

fspace quad {
  {mass_x, mass_y, mass, center_x, center_y, size} : double,
  {total, leaf_count, type, index} : uint,
  {ne, nw, se, sw, next_in_leaf} : int,
}

task build_quad(bodies : region(ispace(int1d), body), quads : region(ispace(int1d), quad), min_x : double, min_y : double, size : double, sector_precision : uint, leaf_size : uint, sector : int1d, partition_size: uint)
  where
  reads(bodies.{mass_x, mass_y, mass, index}),
  reads writes(quads)
do
  var sector_x = sector % sector_precision
  var sector_y: int64 = cmath.floor(sector / sector_precision)

  var index = sector * partition_size
  var root_index = index
  assert(quads[root_index].type == 0, "root already allocated")
  quads[root_index].center_x = min_x + (sector_x + 0.5) * size / sector_precision
  quads[root_index].center_y = min_y + (sector_y + 0.5) * size / sector_precision
  quads[root_index].size = size / sector_precision
  quads[root_index].type = 2
  
  var parent_list : int1d[1024]
  var child_list : int1d[1024]
  var traverse_index = 0

  for body in bodies do
    index = index + 1
    assert(quads[index].type == 0, "body already allocated")
    quads[index].mass_x = body.mass_x
    quads[index].mass_y = body.mass_y
    quads[index].mass = body.mass
    quads[index].total = 1
    quads[index].type = 1
    quads[index].index = body.index

    parent_list[traverse_index] = root_index
    child_list[traverse_index] = index
    
    while traverse_index >= 0 do  
      var parent_index = parent_list[traverse_index]
      var child_index = child_list[traverse_index]
      traverse_index = traverse_index - 1

      var half_size = quads[parent_index].size / 2
      if quads[child_index].mass_x <= quads[parent_index].center_x then
        if quads[child_index].mass_y <= quads[parent_index].center_y then
          if quads[parent_index].sw == -1 then
            quads[child_index].leaf_count = 1
            quads[parent_index].sw = child_index
          elseif quads[quads[parent_index].sw].type == 1 then
            if quads[quads[parent_index].sw].leaf_count < leaf_size then
              quads[child_index].leaf_count = quads[quads[parent_index].sw].leaf_count + 1
              quads[child_index].next_in_leaf = quads[parent_index].sw
              quads[parent_index].sw = child_index
            else
              index += 1
              assert(quads[index].type == 0, "region already allocated")
              quads[index].type = 2
              quads[index].center_x = quads[parent_index].center_x - half_size / 2
              quads[index].center_y = quads[parent_index].center_y - half_size / 2
              quads[index].size = half_size
              quads[parent_index].sw = index

              var current = quads[parent_index].sw
              while current ~= -1 do
                var next_in_leaf = quads[current].next_in_leaf
                quads[current].next_in_leaf = -1
                traverse_index += 1
                parent_list[traverse_index] = index
                child_list[traverse_index] = current
                current = next_in_leaf
              end

              traverse_index += 1
              parent_list[traverse_index] = index
              child_list[traverse_index] = child_index              
            end
          else
            traverse_index += 1
            parent_list[traverse_index] = quads[parent_index].sw
            child_list[traverse_index] = child_index             
          end
        else
          if quads[parent_index].nw == -1 then
            quads[child_index].leaf_count = 1
            quads[parent_index].nw = child_index
          elseif quads[quads[parent_index].nw].type == 1 then
            if quads[quads[parent_index].nw].leaf_count < leaf_size then
              quads[child_index].leaf_count = quads[quads[parent_index].nw].leaf_count + 1
              quads[child_index].next_in_leaf = quads[parent_index].nw
              quads[parent_index].nw = child_index
            else
              index += 1
              assert(quads[index].type == 0, "region already allocated")
              quads[index].type = 2
              quads[index].center_x = quads[parent_index].center_x - half_size / 2
              quads[index].center_y = quads[parent_index].center_y + half_size / 2
              quads[index].size = half_size
              quads[parent_index].nw = index

              var current = quads[parent_index].nw
              while current ~= -1 do
                var next_in_leaf = quads[current].next_in_leaf
                quads[current].next_in_leaf = -1
                traverse_index += 1
                parent_list[traverse_index] = index
                child_list[traverse_index] = current
                current = next_in_leaf
              end

              traverse_index += 1
              parent_list[traverse_index] = index
              child_list[traverse_index] = child_index 
            end
          else
            traverse_index += 1
            parent_list[traverse_index] = quads[parent_index].nw
            child_list[traverse_index] = child_index 
          end      
        end
      else
        if quads[child_index].mass_y <= quads[parent_index].center_y then
          if quads[parent_index].se == -1 then
            quads[child_index].leaf_count = 1
            quads[parent_index].se = child_index
          elseif quads[quads[parent_index].se].type == 1 then
            if quads[quads[parent_index].se].leaf_count < leaf_size then
              quads[child_index].leaf_count = quads[quads[parent_index].se].leaf_count + 1
              quads[child_index].next_in_leaf = quads[parent_index].se
              quads[parent_index].se = child_index
            else
              index += 1
              assert(quads[index].type == 0, "region already allocated")
              quads[index].type = 2
              quads[index].center_x = quads[parent_index].center_x + half_size / 2
              quads[index].center_y = quads[parent_index].center_y - half_size / 2
              quads[index].size = half_size
              quads[parent_index].se = index

              var current = quads[parent_index].se
              while current ~= -1 do
                var next_in_leaf = quads[current].next_in_leaf
                quads[current].next_in_leaf = -1
                traverse_index += 1
                parent_list[traverse_index] = index
                child_list[traverse_index] = current
                current = next_in_leaf
              end

              traverse_index += 1
              parent_list[traverse_index] = index
              child_list[traverse_index] = child_index 
            end
          else
            traverse_index += 1
            parent_list[traverse_index] = quads[parent_index].se
            child_list[traverse_index] = child_index 
          end 
        else
          if quads[parent_index].ne == -1 then
            quads[child_index].leaf_count = 1
            quads[parent_index].ne = child_index
          elseif quads[quads[parent_index].ne].type == 1 then
            if quads[quads[parent_index].ne].leaf_count < leaf_size then
              quads[child_index].leaf_count = quads[quads[parent_index].ne].leaf_count + 1
              quads[child_index].next_in_leaf = quads[parent_index].ne
              quads[parent_index].ne = child_index
            else
              index += 1
              assert(quads[index].type == 0, "region already allocated")
              quads[index].type = 2
              quads[index].center_x = quads[parent_index].center_x + half_size / 2
              quads[index].center_y = quads[parent_index].center_y + half_size / 2
              quads[index].size = half_size
              quads[parent_index].ne = index

              var current = quads[parent_index].ne
              while current ~= -1 do
                var next_in_leaf = quads[current].next_in_leaf
                quads[current].next_in_leaf = -1
                traverse_index += 1
                parent_list[traverse_index] = index
                child_list[traverse_index] = current
                current = next_in_leaf
              end

              traverse_index += 1
              parent_list[traverse_index] = index
              child_list[traverse_index] = child_index 
            end
          else
            traverse_index += 1
            parent_list[traverse_index] = quads[parent_index].se
            child_list[traverse_index] = child_index
          end     
        end
      end

      var old_mass = quads[parent_index].mass
      var new_mass = quads[parent_index].mass + quads[child_index].mass
      quads[parent_index].mass_x = (quads[parent_index].mass_x * old_mass + quads[child_index].mass_x * quads[child_index].mass) / new_mass
      quads[parent_index].mass_y = (quads[parent_index].mass_y * old_mass + quads[child_index].mass_y * quads[child_index].mass) / new_mass
      quads[parent_index].mass = new_mass

      quads[parent_index].total += 1
    end
  end
end
