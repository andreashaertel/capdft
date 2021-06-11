// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_spherical.hpp"
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical() {
  //
}
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical(
    const Properties& system_properties,
    const std::vector<Properties>& species_properties,
    double** density_profile)
  : density_profile(density_profile) {
  // Clear all std::vectors
  diameters.clear();
  bulk_densities.clear();
  // Get system properties
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  // Get species properties; excludes all species without diameter
  double diameter{0.};
  double bulk_density{0.};
  for (auto it = species_properties.begin(); it != species_properties.end();
      ++it) {
    if (it->get_property("diameter", &diameter)) {  // only species with diam.
      if (!it->get_property("bulk density", &bulk_density)) {
        // TODO(Moritz): throw error, because species has diameter, but no dens
      }
      diameters.push_back(diameter);
      bulk_densities.push_back(bulk_density);
    }
  }
  // Count species that interact via the hard sphere potential
  species_count = diameters.size();
  // TODO(Moritz): Initialize the following objects
  // TODO(Moritz): functional_derivative
}
// _____________________________________________________________________________
FunctionalFMTSpherical::~FunctionalFMTSpherical() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_derivative() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_bulk_derivative() {
  //
}
// _____________________________________________________________________________
double FunctionalFMTSpherical::calc_energy() {
  //
  return 0.;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
