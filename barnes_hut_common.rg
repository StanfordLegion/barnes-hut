import "regent"

local c = regentlib.c

fspace quad(r : region(quad(wild))) {
  mass_x: double,
  mass_y: double,
  mass: double,
  center_x: double,
  center_y: double,
  size: double,
  total: uint,
  type: uint,
  ne: ptr(quad(wild), r),
  nw: ptr(quad(wild), r),
  se: ptr(quad(wild), r),
  sw: ptr(quad(wild), r)
}

task add_node(quads : region(quad(wild)), cur: ptr(quad(wild), quads), body: ptr(quad(wild), quads), last_used: uint): uint
where
  reads(quads),
  writes(quads)
do
  c.printf("Inserting cur: %d, body: %d, last_used: %d\n", cur, body, last_used)
  var half_size = cur.size / 2
  var new_last_used = last_used
  if body.mass_x <= cur.center_x then
    if body.mass_y <= cur.center_y then
      if isnull(cur.sw) then
        cur.sw = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.sw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x - half_size / 2
        new_fork.center_y = cur.center_y - half_size / 2
        new_fork.size = half_size
        cur.sw = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.sw), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.sw, body, last_used)
      end
    else
      if isnull(cur.nw) then
        cur.nw = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.nw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x - half_size / 2
        new_fork.center_y = cur.center_y + half_size / 2
        new_fork.size = half_size
        cur.nw = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.nw), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.nw, body, last_used)
      end      
    end
  else
    if body.mass_y <= cur.center_y then
      if isnull(cur.se) then
        cur.se = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.se).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x + half_size / 2
        new_fork.center_y = cur.center_y - half_size / 2
        new_fork.size = half_size
        cur.se = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.se), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.se, body, last_used)
      end
    else
      if isnull(cur.ne) then
        cur.ne = dynamic_cast(ptr(quad(quads), quads), last_used)
      elseif dynamic_cast(ptr(quad(quads), quads), cur.sw).type == 1 then
        var new_fork = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_fork.center_x = cur.center_x + half_size / 2
        new_fork.center_y = cur.center_x + half_size / 2
        new_fork.size = half_size
        cur.ne = dynamic_cast(ptr(quad(quads), quads), last_used + 1)
        new_last_used = add_node(quads, new_fork, dynamic_cast(ptr(quad(quads), quads), cur.ne), last_used + 1)
        new_last_used = add_node(quads, new_fork, body, new_last_used)
      else
        new_last_used = add_node(quads, cur.ne, body, last_used)
      end      
    end
  end  

  return new_last_used
end
