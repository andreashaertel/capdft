// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_spherical.hpp"
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical() {
  //
}
// _____________________________________________________________________________
FunctionalFMTSpherical::FunctionalFMTSpherical(System* system)
  : density_profile(system->get_density_profile_pointer()) {
  // Create temporary Properties objects to transfer information
  Properties system_properties(system->get_system_properties());
  std::vector<Properties> species_properties = system->get_species_properties();
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
        exit(1);  // TODO(Moritz): throw error
      }
      diameters.push_back(diameter);
      bulk_densities.push_back(bulk_density);
    }
  }
  // Count species that interact via the hard sphere potential
  species_count = diameters.size();
  // TODO(Moritz): Initialize the following objects
  // TODO(Moritz): functional_derivative
  // TODO(Moritz): weighted densities
  // TODO(Moritz): density_profile_times_r
  // TODO(Moritz):
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
  double integral{0.};
  DataField free_energy_density(1, grid_count);
  // Calculate the weighted densities for the current density profile
  calc_weighted_densities();
  check_weighted_densities();
  // Store free energy density
  //for (int i = 0; i < grid_no; ++i) {
  //  free_energy_density.at(0, i) = calc_local_energy_density(i);
  //}
  //integral = radial_integration(free_energy_density, grid_no, dr);
  //
  return 0.;  // TODO(Moritz): correct value
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_weighted_densities() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::check_weighted_densities() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTSpherical::calc_partial_derivatives() {
  //
}
// _____________________________________________________________________________
//double FunctionalFMTSpherical::calc_local_energy_density(size_t position) {
//  // Define a lot of dummy variables
//  double Phi1, Phi2, Phi3;  // energy density terms
//  double phi2, phi3;  // factors in the energy density terms
//  double phi2_num, phi3_num;  // factors in the energy density terms
//  double n0, n1, n2, n3; // scalar weighted densities
//  double n3n3, n3n3n3, oneMn3, logOneMn3, oneMn3squared;  // auxiliary variables
//  double n2n2, n2n2n2;  // auxiliary variables
//  double trace2, trace3;  // auxiliary variables
//  double nvec1, nvec2;  // vectorial weighted densities (only z-comp. is non-0)
//  double ntensorm2first;  // tensorial weighted densities (tensor is diagonal)
//  double ntensorm2second;
//  double ntensorm2third;
//  // Assign values to the weighted densities
//  n3 = weight_nrs[0][i];  // n3
//  n2 = weight_nrs[1][i];  // n2
//  n1 = weight_nrs[2][i];  // n1
//  n0 = weight_nrs[3][i];  // n0
//  nvec2 = weight_nrs[4][i];  // nvec2
//  nvec1 = weight_nrs[5][i];  // nvec1
//  ntensorm2first = weight_nrs[6][i];  // tensor
//  ntensorm2second = weight_nrs[6][i];  // tensor
//  ntensorm2third = weight_nrs[7][i];  // tensor
//  // Calculate auxiliary variables
//  n3n3 = n3 * n3;
//  n3n3n3 = n3 * n3 * n3;
//  oneMn3 = 1. - n3;
//  logOneMn3 = log(oneMn3);
//  oneMn3squared = oneMn3 * oneMn3;
//  n2n2 = n2 * n2;
//  n2n2n2 = n2 * n2n2;
//  trace2 = ntensorm2first * ntensorm2first + ntensorm2second * ntensorm2second +
//      ntensorm2third * ntensorm2third;
//  trace3 = ntensorm2first * ntensorm2first * ntensorm2first +
//      ntensorm2second * ntensorm2second * ntensorm2second +
//      ntensorm2third * ntensorm2third * ntensorm2third;
//  // Calculate factors in the energy density terms
//	if (n3 < 1e-5) {  // avoiding logarithm of very small numbers
//    phi2 = 1. + .5 * n3 + .3 * n3n3 + .2 * n3n3n3;  // +O(n^4)
//    phi3 = 1. - .125 * n3 - .05 * n3n3 - .025 * n3n3n3;  // +O(n^4)
//  } else {
//    phi2 = (6. * n3 - 3. * n3n3 + 6. * oneMn3 * logOneMn3) / n3n3n3;
//    phi3 = (6. * n3 - 9. * n3n3 + 6. * n3n3n3 + 6. * oneMn3squared * logOneMn3) /
//        (4. * n3n3n3);
//  }
//  phi2_num = 1. + n3n3 * phi2 / 9.;
//  phi3_num = 1. - 4. * n3 * phi3 / 9.;
//  // Calculate free energy density terms
//  Phi1 = -n0 * logOneMn3;
//  Phi2 = phi2_num *  (n1 * n2 - nvec1 * nvec2) / oneMn3;
//  Phi3 = phi3_num * (n2n2n2 - 3. * n2 * nvec2 * nvec2 + 4.5 * (
//      n2n2 * ntensorm2third - n2 * nvec2 * nvec2 - trace3 + n2 * trace2)) /
//      (24. * M_PI * oneMn3squared);
//  return Phi1 + Phi2 + Phi3;
//}
// _____________________________________________________________________________
