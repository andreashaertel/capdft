// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_planar.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "data_frame.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalFMTPlanar::FunctionalFMTPlanar() {
  //
}
// _____________________________________________________________________________
FunctionalFMTPlanar::FunctionalFMTPlanar(
    const std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties,
    const std::vector<size_t>& affected_species)
  : affected_species(affected_species),
    density_profiles_pointer(density_profiles) {
  // Clear all std::vectors
  diameters.clear();
  bulk_densities.clear();
  // Get system properties
  extract_system_properties(system_properties);
  // Get species properties; excludes all species without diameter property
  extract_species_properties(species_properties);
  // Initialize all data frames
  initialize_all_data_frames();
  // Calculate weights
  calc_weights();
}
// _____________________________________________________________________________
FunctionalFMTPlanar::FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties)
  : FunctionalFMTPlanar(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalFMTPlanar::~FunctionalFMTPlanar() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::extract_system_properties(
    const Properties& system_properties) {
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dz = length / static_cast<double>(grid_count);
  dkz = 2. * M_PI / length;
  // Calculate normalization factors for Fourier transforms
  // TODO(Moritz)
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::extract_species_properties(
    const std::vector<Properties>& species_properties) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  const std::vector<Properties>& spec_prop = species_properties;
  double diameter{0.};
  double bulk_density{0.};
  if (affected_species.empty()) {
    for (auto it = spec_prop.begin(); it != spec_prop.end(); ++it) {
      if (it->get_property("diameter", &diameter)) {  // only species with diam.
        if (!it->get_property("bulk density", &bulk_density)) {
          std::cerr << "FunctionalFMTPlanar::extract_species_properties(): ";
          std::cerr << "\"Error: A species with a diameter but no density ";
          std::cerr << "was detected.\"" << std::endl;
          exit(1);
        }
        affected_species.push_back(it - spec_prop.begin());
      }
    }
  }
  // Extract properties
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (!spec_prop.at(*it).get_property("diameter", &diameter) ||
        !spec_prop.at(*it).get_property("bulk density", &bulk_density)) {
      std::cerr << "FunctionalFMTPlanar::extract_species_properties(): ";
      std::cerr << "\"Error: One species is lacking a required parameter.";
      std::cerr << std::endl;
      exit(1);
    }
    diameters.push_back(diameter);
    bulk_densities.push_back(bulk_density);
  }
  // Count species that interact via the hard sphere potential
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::initialize_all_data_frames() {
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_weights() {
  // Declare auxiliary variables. R means radius of the hard sphere.
  // kz is the z-component of the position in Fourier space.
  // RR is the square of R, kzkz is the square of kz, and so on.
  double kz{0.}, kzkz{0.}, kzkzkz{0.};
  double R{0.}, RR{0.}, RRR{0.};
  double Rkz{0.};
  //for (size_t i = 0; i != species_count; ++i) {
  //  R = diameters.at(i) / 2.;
  //  RR = R * R;
  //  RRR = RR * R;
  //  // Calculate limits kr --> 0
  //  // scalar weight 3 (theta)
  //  weights_four.at(i).at(0).at(0) = 4. * M_PI * RRR / 3.;
  //  // scalar weight 2 (delta)
  //  weights_four.at(i).at(1).at(0) = 4. * M_PI * RR;
  //  // tensorial
  //  weights_four.at(i).at(2).at(0) =
  //      weights_four.at(i).at(1).at(0) -
  //      3. * weights_four.at(i).at(0).at(0) / R;
  //  for (size_t j = 1; j != grid_count + 1; ++j) {
  //    kr = dkr * static_cast<double>(j);
  //    Rkr = R * kr;
  //    krkr = kr * kr;
  //    krkrkr = krkr * kr;
  //    // scalar weight 3 (theta)
  //    weights_four.at(i).at(0).at(j) = 4. * M_PI * (sin(Rkr) - Rkr * cos(Rkr)) /
  //        krkrkr;
  //    // scalar weight 2 (delta)
  //    weights_four.at(i).at(1).at(j) = 4. * M_PI * R * sin(Rkr) / kr;
  //    // tensorial
  //    weights_four.at(i).at(2).at(j) =
  //        weights_four.at(i).at(1).at(j) -
  //        3. * weights_four.at(i).at(0).at(j) / R;
  //  }
  //}
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != grid_count) {
      std::cerr << "FunctionalFMTPlanar::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  //// Calculate the weighted densities
  //calc_weighted_densities();
  //// From the weighted densities calculate the partial derivatives of the
  //// excess free energy
  //calc_partial_derivatives();
  //// Calculate the weighted partial derivatives
  //calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_derivative_warnings(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != grid_count) {
      std::cerr << "FunctionalFMTPlanar::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  //// Calculate the weighted densities
  //calc_weighted_densities();
  //// Check whether there are unphysical values in the  weighted densities
  //check_weighted_densities();
  //// From the weighted densities calculate the partial derivatives of the
  //// excess free energy
  //calc_partial_derivatives();
  //// Calculate the weighted partial derivatives
  //calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  //// Auxiliary variables
  //double diameter{0.}, diameter2{0.}, diameter3{0.};
  //// Calculate the fluid weighted densities and save them to the first position
  //for (auto it = scalar_weighted_dens_real.begin();
  //    it != scalar_weighted_dens_real.end(); ++it) {
  //  it->at(0) = 0.;
  //}
  //for (size_t i = 0; i != species_count; ++i) {
  //  // Define auxiliary values
  //  diameter = diameters.at(i);
  //  diameter2 = diameter * diameter;
  //  diameter3 = diameter2 * diameter;
  //  // All non-scalar weighted densities are 0 the others are obviously constant
  //  scalar_weighted_dens_real.at(0).at(0) +=
  //      bulk_densities.at(i) * diameter3 * M_PI / 6.0;
  //  scalar_weighted_dens_real.at(1).at(0) +=
  //      bulk_densities.at(i) * diameter2 * M_PI;
  //  scalar_weighted_dens_real.at(2).at(0) +=
  //      bulk_densities.at(i) * diameter / 2.0;
  //  scalar_weighted_dens_real.at(3).at(0) += bulk_densities.at(i);
  //}
  //// Evaluate the partial derivatives at the first position
  //calc_local_partial_derivatives(0);
  //std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
  //size_t spec_i{0};
  //for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
  //  spec_i = it - affected_species.begin();
  //  // Define auxiliary values
  //  diameter = diameters.at(spec_i);
  //  diameter2 = diameter * diameter;
  //  diameter3 = diameter2 * diameter;
  //  // Calculate the derivatives of the excess free energy functional w.r.t. the
  //  // weighted densities.
  //  bulk_derivative->at(*it) = scalar_partial_derivative_real.at(0).at(0) *
  //      diameter3 * M_PI / 6.0;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(1).at(0) *
  //      diameter2 * M_PI;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(2).at(0) *
  //      diameter / 2.0;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(3).at(0);
  //}
}
// _____________________________________________________________________________
double FunctionalFMTPlanar::calc_energy() {
  //double integral{0.};
  //DataFrame<1, double> free_energy_density(grid_count);
  //// Calculate the weighted densities for the current density profile
  //calc_weighted_densities();
  //check_weighted_densities();
  //// Store free energy density
  //for (size_t i = 0; i < grid_count; ++i) {
  //  free_energy_density.at(i) = calc_local_energy_density(i);
  //}
  //integral = radial_integration(free_energy_density.array(), grid_count, dr);
  //return integral;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
