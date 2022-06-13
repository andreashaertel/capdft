// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "dft.hpp"  // NOLINT
#include <cmath>
#include "data_frame.hpp"  // NOLINT
#include "functional_fmt_spherical.hpp"  // NOLINT
// _____________________________________________________________________________
Dft::Dft(System<DataFrame<double>>* system)  // is not solely input variable
  : system(system) {
  // TODO(Andreas): add a functional_ideal object
}
// _____________________________________________________________________________
Dft::~Dft() {
}
// _____________________________________________________________________________
void Dft::add_excess_functional() {
  functionals.push_back(FunctionalFMTSpherical(*system));
}
// _____________________________________________________________________________
void Dft::remove_excess_functional(size_t index) {
  functionals.erase(functionals.begin() + index);
}
// _____________________________________________________________________________
void Dft::set_fugacities_from_bulk_densities(
      std::vector<double>* bulk_densities) {
  // Determine fugacities z_nu according to the Euler equation:
  //  rho(r) = z_nu * e^(c^(1)(r)-beta v(r))
  // Set c1 for each species derivative to zero.
  size_t species_count = system->get_species_properties().size();
  std::vector<double> fugacities(species_count, 0.);
  std::vector<double> direct_corr_1(species_count, 0.);  // c^(1)
  std::vector<double> bulk_derivatives(species_count, 0.);
  // Calculate bulk c1 for each species
  for (auto it = functionals.begin(); it != functionals.end(); ++it) {
    // Determine bulk derivatives from each functional
    it->calc_bulk_derivative(&bulk_derivatives);
    for (size_t i = 0; i < species_count; ++i) {
      // The direct correlation is defined as the negative functional derivative
      // and a prefactor 1/kT. The functional is stored in units of kT.
      direct_corr_1.at(i) -= bulk_derivatives.at(i);
    }
  }
  // Calculate the fugacities
  for (size_t i = 0; i < species_count; ++i) {
    fugacities.at(i) = bulk_densities->at(i) * exp(-direct_corr_1.at(i));
  }
  // Set/Update fugacities
  system->set_fugacities(fugacities);
}
// _____________________________________________________________________________
void Dft::set_fugacities_from_bulk_densities() {
  // Store bulk densities in a vector and call
  // set_fugacities_from_bulk_densities(std::vector<double>* bulk_densities)
  std::vector<double> bulk_densities;
  const std::vector<Properties> species_properties =
      system->get_species_properties();
  for (auto it = species_properties.begin();
        it != species_properties.end(); ++it) {
    double value;
    if (it->get_property<double>("bulk_density", &value))
      bulk_densities.push_back(value);
  }
  set_fugacities_from_bulk_densities(&bulk_densities);
}
// _____________________________________________________________________________
double Dft::iterate_densities() {
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // TODO(Andreas): + Picard mixing parameter, now hard coded.
  //                + Check for numerical correct solution, maybe first
  //                  calculate the result, then check, then store, ...
  // AH, 24.02.2022
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Follow the Euler equation to calculate new density profiles:
  //  rho(r) = z_nu * e^(c^(1)(r; [rho_old])-beta v(r))
  // Mixing: rho_new(r) = a * rho_new + (1.0-a) * rho_old
  double mixing_factor = 0.05;
  // TODO(Andreas):
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
