// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_spherical.hpp"
#include "dft.hpp"
// _____________________________________________________________________________
Dft::Dft(System* system) {
  // initialize internal parameters
  functional_index = 0;
  // Store the system. 
  this->system = system;
  // add a functional_ideal object
  //AH:TODO
}
// _____________________________________________________________________________
Dft::~Dft() {
  // Remove all added functionals
  for (std::map<size_t,Functional*>::iterator it = functional.begin(), 
        next_it = it; it != functional.cend(); it = next_it) {
    ++next_it;
    // destroy functional object
    delete(it->second);
    // erase element from map
    functional.erase(it);
  }
}
// _____________________________________________________________________________
template<typename AnyFunctional>
size_t Dft::add_excess_functional() {
  // Add 1 to the lastly used index to obtain the new index. 
  size_t new_index = functional_index + 1;
  // All functionals are initialized by providing the physical system. 
  AnyFunctional* new_functional = new AnyFunctional(system);
  //x->show();
  //delete(x);
  functional.insert( std::pair<size_t,Functional*>(new_index, new_functional) );
  return new_index;
}
// _____________________________________________________________________________
// Definition of all specific implementations of functionals that can be add to
// a Dft object. 
template size_t Dft::add_excess_functional<FunctionalFMTSpherical>();
// _____________________________________________________________________________
bool Dft::remove_excess_functional(size_t index) {
  // Search for the functional the index is mapped to
  std::map<size_t,Functional*>::iterator it = functional.find(index);
  // If no mapping is found, return false
  if (it == functional.end())
    return false;
  // Delete respective functional and remove it from the map. 
  delete(it->second);
  functional.erase(it);
  return true;
}
// _____________________________________________________________________________
void Dft::set_chempots_from_bulk_densities(
      std::vector<double>* bulk_densities) {
  //
}
// _____________________________________________________________________________
void Dft::set_chempots_from_bulk_densities() {
  //
}
// _____________________________________________________________________________
double Dft::iterate_densities() {
  //
  return 0.0;
}
// _____________________________________________________________________________
double Dft::calculate_gc_energy() {
  //
  return 0.0;
}
// _____________________________________________________________________________
