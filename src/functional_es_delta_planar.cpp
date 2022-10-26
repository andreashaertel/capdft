// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file functional_es_delta_planar.cpp
 *  \brief Source file for the FunctionalESDeltaPlanar class.
 *
 *  The file contains the class definitions of the FunctionalESDeltaPlanar
 *  class.
 */
#include "functional_es_delta_planar.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include "constants.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "integration.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalESDeltaPlanar::FunctionalESDeltaPlanar() {
}
// _____________________________________________________________________________
FunctionalESDeltaPlanar::FunctionalESDeltaPlanar(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties,
    std::vector<size_t> affected_species)
  : affected_species(affected_species),
    density_profiles_pointer(density_profiles) {
  // Get system properties
  extract_system_properties(system_properties);
  // Get species properties
  extract_species_properties(species_properties);
  // Initialize all data frames and update charge densities
  calc_extended_system_bounds();
  initialize_all_data_frames();
  calc_charge_densities();
  initialize_poisson_solver();
  initialize_weights();
}
// _____________________________________________________________________________
FunctionalESDeltaPlanar::FunctionalESDeltaPlanar(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties)
  : FunctionalESDeltaPlanar(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalESDeltaPlanar::~FunctionalESDeltaPlanar() {
  delete poisson_solver;
  delete total_poisson_solver;
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::extract_system_properties(
    const Properties& system_properties) {
  // Extract system properties directly
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dz = length / static_cast<double>(grid_count);
  // Calculate electrical properties of the system
  extract_electrical_properties(system_properties);
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::extract_electrical_properties(
    const Properties& system_properties) {
  // This variable saves a number representing, which of the el. variables are
  // available.
  int cases{0};
  // Get the variables
  try {
    if (system_properties.get_property("temperature", &temperature)) {
      cases += 1;
    }
  } catch (...) {}
  try {
    if (system_properties.get_property("bjerrum length", &bjerrum)) {
      cases += 2;
    }
  } catch (...) {}
  try {
    if (system_properties.get_property("dielectric constant", &dielectric)) {
      cases += 4;
    }
  } catch (...) {}
  // Auxiliary variable e^2 / (4 * \epsilon_0 * \pi * k_B)
  double e2_by_4pi_boltz = 1e-3 * ELECTRON_CHARGE * ELECTRON_CHARGE /
      (4. * M_PI * VACUUM_PERMITIVITY * BOLTZMANN);
  // Calculate missing variables
  switch (cases) {
    case 3:  // dielectric constant missing
      dielectric = e2_by_4pi_boltz / (bjerrum * 1e-9 * temperature);
      break;
    case 5:  // bjerrum length missing
      bjerrum = 1e9 * e2_by_4pi_boltz / (dielectric * temperature);
      break;
    case 6:  // temperature missing
      temperature = e2_by_4pi_boltz / (dielectric * bjerrum * 1e-9);;
      break;
    case 7:  // nothing missing
      std::cerr << "FunctionalESDeltaPlanar::extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "FunctionalESDeltaPlanar::extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::extract_species_properties(
    const std::vector<Properties>& species_properties) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  double diameter{0.};
  double bulk_density{0.};
  double valency{0.};
  if (affected_species.empty()) {
    for (auto& properties : species_properties) {
      try {
        if (properties.get_property("valency", &valency) &&
            properties.get_property("diameter", &diameter)) {
          affected_species.push_back(&properties - &species_properties[0]);
        }
      } catch(...) {}
    }
  }
  // Get species specific properties like diameter and valency
  for (auto& species : affected_species) {
    species_properties.at(species).get_property("diameter", &diameter);
    species_properties.at(species).get_property("bulk density", &bulk_density);
    species_properties.at(species).get_property("valency", &valency);
    diameters.push_back(diameter);
    bulk_densities.push_back(bulk_density);
    valencies.push_back(valency);
  }
  // Count species that interact via this functional
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_extended_system_bounds() {
  if (diameters.size() == 0) {
    std::cerr << "FunctionalESDeltaPlanar::calc_extended_system_bounds():";
    std::cerr << "\"ERROR: Diameters have not been initialized yet.\"";
    std::cerr << std::endl;
    exit(1);
  }
  double max_diameter{0.};
  max_diameter = *std::max_element(diameters.begin(), diameters.end());
  extended_system_offset = static_cast<size_t>(max_diameter / dz + .5);
  extended_grid_count = grid_count + 2 * extended_system_offset;
  grid_count_four = extended_grid_count / 2 + 1;
  dkz = 2. * M_PI / (static_cast<double>(extended_grid_count) * dz);
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::initialize_all_data_frames() {
  if (extended_grid_count < grid_count) {
    std::cerr << "FunctionalESDeltaPlanar::initialize_all_data_frames():";
    std::cerr << "\"ERROR: Extended system has not been initialized yet.\"";
    std::cerr << std::endl;
    exit(1);
  }
  total_charge_density_profile.~DataFrame<1, double>();
  new(&total_charge_density_profile) DataFrame<1, double>(grid_count);
  total_poisson_rhs.~DataFrame<1, double>();
  new(&total_poisson_rhs) DataFrame<1, double>(grid_count);
  total_potential.~DataFrame<1, double>();
  new(&total_potential) DataFrame<1, double>(grid_count);
  charge_density_profiles =
      std::vector(species_count, DataFrame<1, double>(extended_grid_count));
  poisson_rhs =
      std::vector(species_count, DataFrame<1, double>(extended_grid_count));
  potentials = 
      std::vector(species_count, DataFrame<1, double>(extended_grid_count));
  weighted_densities_real = 
      std::vector(species_count, DataFrame<1, double>(extended_grid_count));
  weights_four = std::vector(species_count,
      std::vector(species_count,
      DataFrame<1, fftw_complex>(grid_count_four)));
  weighted_densities_four = std::vector(species_count,
      std::vector(species_count,
      DataFrame<1, fftw_complex>(grid_count_four)));
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::initialize_poisson_solver() {
  poisson_solver = new PlanarPoissonSolver(
      extended_grid_count, dz, DIRICHLET_DIRICHLET);
  total_poisson_solver = new PlanarPoissonSolver(
      grid_count, dz, DIRICHLET_DIRICHLET);
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::initialize_weights() {
  double kx{0.}, ky{0.}, kz{0.}, kabs{0.};
  double average_diameter{0.};
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 0; j < species_count; ++j) {
      average_diameter = .5 * (diameters.at(i) + diameters.at(j));
      for (size_t k = 0; k < grid_count_four; ++k) {
        kz = dkz * static_cast<double>(k);
        kabs = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
        if (k == 0) {
          weights_four.at(i).at(j).at(k)[0] = 1.;
          weights_four.at(i).at(j).at(k)[0] = 0.;
        } else {
          weights_four.at(i).at(j).at(k)[0] =
              sin(kabs * average_diameter) / (kabs * average_diameter);
          weights_four.at(i).at(j).at(k)[1] = 0.;
        }
      }
    }
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_charge_densities() {
  size_t index{0};  // index of the extended system
  total_charge_density_profile.zero();
  for (size_t i = 0; i < species_count; ++i) {
    total_charge_density_profile += valencies.at(i) *
          density_profiles_pointer->at(affected_species.at(i));
    charge_density_profiles.at(i).zero();
    for (size_t j = 0; j < grid_count; ++j) {
      index = extended_system_offset + j;
      charge_density_profiles.at(i).at(index) = valencies.at(i) *
          density_profiles_pointer->at(affected_species.at(i)).at(j);
    }
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  size_t i{0};
  // From the charge densities calculate the weighted densities
  // The charge densities are calculated after the FFTW plans are established.
  calc_weighted_densities();
  calc_poisson_rhs();
  // From the charge densities calculate the electrostatic potential
  calc_potential();
  // Calculate the derivative from the electrostatic potential
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    i = it - affected_species.begin();
    for (size_t j = 0; j < grid_count; ++j) {
      functional_derivative->at(*it).at(j) = valencies.at(i) *
          potentials.at(i).at(extended_system_offset+j);
    }
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // The bulk values of this functional's derivative are always zero
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);  // TODO(Moritz): correct
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_weighted_densities() {
  // Create fftw plans
  std::vector<fftw_plan> plans_forward;
  std::vector<fftw_plan> plans_backward;
  for (size_t i = 0; i < species_count; ++i) {
    plans_forward.push_back(
        fftw_plan_dft_r2c_1d(
            extended_grid_count, charge_density_profiles.at(i).array(),
            weighted_densities_four.at(i).at(0).array(), FFTW_PATIENT));
    plans_backward.push_back(
        fftw_plan_dft_c2r_1d(
            extended_grid_count, weighted_densities_four.at(0).at(i).array(),
            weighted_densities_real.at(i).array(), FFTW_PATIENT));
  }
  calc_charge_densities();
  // Execute Fourier transform
  for (auto& plan : plans_forward) {
    fftw_execute(plan);
  }
  // Normalize
  for(auto& weighted_density_four : weighted_densities_four) {
    weighted_density_four.at(0) *= dz;
  }
  // Extend data
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 1; j < species_count; ++j) {  // j=1
      weighted_densities_four.at(i).at(j) = weighted_densities_four.at(i).at(0);
    }
  }
  // Convolution of charge densities and weights
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 0; j < species_count; ++j) {
      weighted_densities_four.at(i).at(j) *= weights_four.at(i).at(j);
    }
  }
  // Summation over all species of the weight functions before back transforming
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 1; j < species_count; ++j) {  // j=1
      weighted_densities_four.at(0).at(i) +=
          weighted_densities_four.at(j).at(i);
    }
  }
  // Back transform into real space
  for (auto& plan : plans_backward) {
    fftw_execute(plan);
  }
  // Normalization
  for (auto& weighted_density : weighted_densities_real) {
    weighted_density *= dkz / (2. * M_PI);
  }
  for (auto& plan : plans_forward) { fftw_destroy_plan(plan); }
  for (auto& plan : plans_backward) { fftw_destroy_plan(plan); }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_poisson_rhs() {
  total_poisson_rhs = -4. * M_PI * bjerrum * total_charge_density_profile;
  for (size_t i = 0; i < species_count; ++i) {
    poisson_rhs.at(i) = -4. * M_PI * bjerrum * weighted_densities_real.at(i);
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaPlanar::calc_potential() {
  // Both boundary conditions equal 0 (Dirichlet). Note, that there is no
  // external potential involved inside the functional.
  double left_boundary{0.}, left_boundary_extended{0.};
  double right_boundary{0.}, right_boundary_extended{0.};
  // Solve the regular Poisson equation numerically
  total_poisson_solver->solve(left_boundary, right_boundary,
      total_poisson_rhs.array(), total_potential.array());
  // Calculate the boundary conditions via extrapolation
  // Keep in mind, that the boundary condition is set to be half a bin outside
  // the system.
  left_boundary_extended = left_boundary +
      static_cast<double>(extended_system_offset) *
      2. * (left_boundary - total_potential.at(0));
  right_boundary_extended = right_boundary +
      static_cast<double>(extended_system_offset) *
      2. * (right_boundary - total_potential.at(grid_count - 1));
  // Now solve the Poisson equation for the weighted densities separately
  // TODO(Moritz): set potential to zero at regular boundaries
  for (size_t i = 0; i < species_count; ++i) {
    poisson_solver->solve(left_boundary_extended, right_boundary_extended,
        poisson_rhs.at(i).array(), potentials.at(i).array());
  }
}
// _____________________________________________________________________________
double FunctionalESDeltaPlanar::calc_energy() {
  DataFrame<1, double> energy_density(grid_count);
  double integral{0.};
  calc_charge_densities();
  //calc_potential();
  //energy_density = .5 * potential * charge_density_profile;
  //integral = integration_1d_closed(energy_density, dz);
  return integral;
}
// _____________________________________________________________________________
