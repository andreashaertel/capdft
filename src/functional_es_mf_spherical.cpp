// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "src/functional_es_mf_spherical.hpp"
// _____________________________________________________________________________
FunctionalESMFSpherical::FunctionalESMFSpherical(
    const Properties& system_properties,
    const std::vector<Properties>& species_properties,
    double** density_profile)
  : density_profile(density_profile) {
  // Clear all std::vectors
  valencies.clear();
  diameters.clear();
  bulk_densities.clear();
  // Get system properties
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  system_properties.get_property("temperature", &temperature);
  system_properties.get_property("dielectric constant", &dielectric);
  system_properties.get_property("bjerrum length", &bjerrum_length);
  // TODO(Moritz): cases for bjerrum length and temperature + dielectric
  // Get species properties; excludes all species without diameter
  double diameter{0.};
  double bulk_density{0.};
  double valency{0.};
  for (auto it = species_properties.begin(); it != species_properties.end();
      ++it) {
    if (it->get_property("valency", &valency)) {  // only species with valency
      if (!it->get_property("bulk density", &bulk_density)) {
        // TODO(Moritz): throw error, because species has valency, but no dens
      } else {
        bulk_densities.push_back(bulk_density);
      }
      if (!it->get_property("diameter", &diameter)) {
        diameters.push_back(0.);
      } else {
        diameters.push_back(diameter);
      }
    }
  }
  // Count species that interact via the hard sphere potential
  species_count = diameters.size();
  // TODO(Moritz): Initialize the following objects
  // TODO(Moritz): functional_derivative
  // TODO(Moritz): charge_density_profile
  // TODO(Moritz): species_potentials
}
// _____________________________________________________________________________
FunctionalESMFSpherical::~FunctionalESMFSpherical() {
  //
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_derivative(
    DataField<double>* functional_derivative) {
  //
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  //
}
// _____________________________________________________________________________
double FunctionalESMFSpherical::calc_energy() {
  //
  return 0.;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
