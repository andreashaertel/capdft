

Ideas for discussion by AH, 2021


?? Constructor(Density profile)



functional_vector


register_functional(Functional new_functional)
  check, wether densities are equal

register extpot
  check, wether densities are equal

set_iteration_method()




iteration_step(density, extpot)
  calculate functional derivatives -> density objects
  add density objects

  recalculate density using mixing
  for all functionals: update <-- ??? wichtig, oder egal? (calc potential)





Constructor extpot(density profile)
  calc/return/... exp(-extpot)
  calc/return/... extpot
  if infinity, throw infinity exception



density object: 
  needs +
  needs -
  needs *
  needs exp()
They only work, if the objects are equal (same properties, same functional kind) 




