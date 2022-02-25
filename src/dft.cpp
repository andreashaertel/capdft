// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <math.h>
#include "functional_fmt_spherical.hpp"
#include "dft.hpp"
// _____________________________________________________________________________
Dft::Dft(System* system):system(system) {
  // initialize internal parameters
  functional_index = 0;
  // Store the system. 
  //this->system = system;
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
void Dft::set_fugacities_from_bulk_densities(
      std::vector<double>* bulk_densities) {
  // Determine fugacities z_nu according to the Euler equation: 
  //  rho(r) = z_nu * e^(c^(1)(r)-beta v(r))
  // Set c1 for each species derivative to zero. 
  std::vector<double> fugacities;
  std::vector<double> direct_corr_c1;
  std::vector<double> bulk_derivatives;
  for (size_t i = 0; i < system->get_species_properties().size(); i++) {
    fugacities.push_back(0.0);
    direct_corr_c1.push_back(0.0);
    bulk_derivatives.push_back(0.0);
  }
  // Determine bulk derivatives from each functional
  for (auto it = functional.begin(); it != functional.end(); ++it) {
    it->second->calc_bulk_derivative(&bulk_derivatives);
    auto it2 = bulk_derivatives.begin();
    for (auto it1 = direct_corr_c1.begin(); it1 != direct_corr_c1.end(); ++it1) 
                      {
      // The direct correlation is defined as the negative functional derivative
      // and a prefactor 1/kT. The functional is stored in units of kT. 
      *it1 = *it1 - *it2;
      *it2 = 0.0;
      ++it2;
    }
  }
  // Calculate the fugacities
  auto dens_it = bulk_densities->begin();
  auto c1_it = direct_corr_c1.begin();
  for (auto it = fugacities.begin(); it != fugacities.end(); ++it) {
    *it = *dens_it / exp(*c1_it);
    ++dens_it;
    ++c1_it;
  }
  // Set/Update fugacities
  system->set_fugacities(&fugacities);
}
// _____________________________________________________________________________
void Dft::set_fugacities_from_bulk_densities() {
  // Store bulk densities in a vector and call 
  // set_fugacities_from_bulk_densities(std::vector<double>* bulk_densities)
  std::vector<double> bulk_densities;
  const std::vector<Properties> species_properties 
        = system->get_species_properties();
  for (auto it = species_properties.begin();
        it != species_properties.end(); ++it) {
    double value;
    if (it->get_property<double>("bulk_density", &value))
      bulk_densities.push_back(value);
  }
  this->set_fugacities_from_bulk_densities(&bulk_densities);
}
// _____________________________________________________________________________
double Dft::iterate_densities() {
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // TODO: + Picard mixing parameter, now hard coded. 
  //       + Check for numerical correct solution, maybe first calculate the 
  //         result, then check, then store, ... 
  // AH, 24.02.2022
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Follow the Euler equation to calculate new density profiles: 
  //  rho(r) = z_nu * e^(c^(1)(r; [rho_old])-beta v(r))
  // Mixing: rho_new(r) = a * rho_new + (1.0-a) * rho_old
  double mixing_factor = 0.05;
  // TODO: 
  // copy construct functional derivatives AND c1 AND new densities 
  //                                                    from density profiles. 
  // iterate over all functionals and add up derivatives. 
  // get exp(beta external potential) from the system 
  //       -> the system has to handle the external potential. 
  // Determine the new densities. 
  // Determine the deviation of both densities (a norm operator in DataField)
  // mix densities and store. 
  return 0.0;
}
// _____________________________________________________________________________
double Dft::calculate_gc_energy() {
  //
  return 0.0;
}
// _____________________________________________________________________________
