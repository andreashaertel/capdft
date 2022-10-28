// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file functional_es_delta_spherical.cpp
 *  \brief Source file for the FunctionalESDeltaSpherical class.
 *
 *  The file contains the definitions of the FunctionalESDeltaSpherical class.
 */
#include "functional_es_delta_spherical.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "constants.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "integration.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalESDeltaSpherical::FunctionalESDeltaSpherical() {
}
// _____________________________________________________________________________
FunctionalESDeltaSpherical::FunctionalESDeltaSpherical(
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
  // Initialize all data frames and update charge densities etc.
  initialize_all_data_frames();
  calc_charge_densities();
  initialize_poisson_solver();
  initialize_weights();
}
// _____________________________________________________________________________
FunctionalESDeltaSpherical::FunctionalESDeltaSpherical(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties)
  : FunctionalESDeltaSpherical(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalESDeltaSpherical::~FunctionalESDeltaSpherical() {
  delete poisson_solver;
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::extract_system_properties(
    const Properties& system_properties) {
  // Extract system properties directly
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dr = length / static_cast<double>(grid_count);
  dkr = (2. * M_PI) / (2. * static_cast<double>(grid_count+1) * dr);
  // Calculate electrical properties of the system
  extract_electrical_properties(system_properties);
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::extract_electrical_properties(
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
      std::cerr << "FunctionalESDeltaSpherical::extract_electrical_properties(";
      std::cerr << "): \"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "FunctionalESDeltaSpherical::extract_electrical_properties(";
      std::cerr << "): \"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::extract_species_properties(
    const std::vector<Properties>& species_properties) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  double valency{0.};
  double bulk_density{0.};
  double diameter{0.};
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
  // Get valencies and diameters
  for (auto& species : affected_species) {
    species_properties.at(species).get_property("valency", &valency);
    species_properties.at(species).get_property("bulk density", &bulk_density);
    species_properties.at(species).get_property("diameter", &diameter);
    valencies.push_back(valency);
    bulk_densities.push_back(bulk_density);
    diameters.push_back(diameter);
  }
  // Count species that interact via the electrostatic forces
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::initialize_all_data_frames() {
  charge_density_profiles = std::vector<DataFrame<1, double>>(species_count,
      DataFrame<1, double>(grid_count));
  potentials = std::vector<DataFrame<1, double>>(species_count,
      DataFrame<1, double>(grid_count));
  poisson_rhs = std::vector<DataFrame<1, double>>(species_count,
      DataFrame<1, double>(grid_count));
  weights_delta = std::vector<std::vector<DataFrame<1, double>>>(species_count,
      std::vector<DataFrame<1, double>>(species_count,
          DataFrame<1, double>(grid_count)));
  weighted_densities = std::vector<std::vector<DataFrame<1, double>>>(
      species_count, std::vector<DataFrame<1, double>>(species_count,
          DataFrame<1, double>(grid_count)));
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::initialize_poisson_solver() {
  poisson_solver = new RadialPoissonSolver(grid_count, dr, .5 * dr);
  poisson_solver->set_radial_laplacian(NEUMANN_DIRICHLET);
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::initialize_weights() {
  double kr = 0.;
  double average_diameter = 0.;
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 0; j < species_count; ++j) {
      average_diameter = .5 * (diameters.at(i) + diameters.at(j));
      for (size_t k = 0; k < grid_count; ++k) {
        kr = dkr * static_cast<double>(k+1);
        weights_delta.at(i).at(j).at(k) = sin(kr * average_diameter) /
            (kr * average_diameter);
      }
    }
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_charge_densities() {
  for (size_t i = 0; i < species_count; ++i) {
    charge_density_profiles.at(i) = valencies.at(i) *
        density_profiles_pointer->at(affected_species.at(i));
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_poisson_rhs() {
  for (size_t i = 0; i < species_count; ++i) {
    poisson_rhs.at(i) = -4. * M_PI * bjerrum * weighted_densities.at(0).at(i);
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_potentials() {
  double outer_boundary{0.};
  double inner_boundary{0.};
  std::vector<double> total_charges{calc_charge_weight_dens()};
  for (size_t i = 0; i < species_count; ++i) {
    inner_boundary = 0.;
    outer_boundary = bjerrum *  total_charges.at(i) /
        (dr * static_cast<double>(grid_count+1));
    // Solve the radial Poisson equation numerically
    poisson_solver->solve(inner_boundary, outer_boundary,
        poisson_rhs.at(i).array(), potentials.at(i).array());
  }
}
// _____________________________________________________________________________
std::vector<double> FunctionalESDeltaSpherical::calc_charge_weight_dens() {
  double integral{0.};
  double r{0.};
  std::vector<double> integrals(0);
  for (auto& weighted_density : weighted_densities.at(0)) {
    integral = 0.;
    // Trapezoidal integral rule, spherical integral
    r = static_cast<double>(0 + 1) * dr;
    integral += .5 * r * r * weighted_density.at(0) * dr;
    for (size_t i = 1; i < grid_count - 1; ++i) {
      r = static_cast<double>(i + 1) * dr;
      integral += r * r * weighted_density.at(i) * dr;
    }
    integral += .5 * r * r * dr * weighted_density.at(grid_count - 1);
    // Spherical symmetry
    integral *= 4. * M_PI;
    integrals.push_back(integral);
  }
  return integrals;
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_weighted_densities() {
  double r = 0.;
  // Create fftw plans
  std::vector<fftw_plan> plans_forward;
  std::vector<fftw_plan> plans_backward;
  for (size_t i = 0; i < species_count; ++i) {
    plans_forward.push_back(
        fftw_plan_r2r_1d(
            grid_count, weighted_densities.at(i).at(0).array(),
            weighted_densities.at(i).at(0).array(), FFTW_RODFT00,
            FFTW_PATIENT));
    plans_backward.push_back(
        fftw_plan_r2r_1d(
            grid_count, weighted_densities.at(0).at(i).array(),
            weighted_densities.at(0).at(i).array(), FFTW_RODFT00,
            FFTW_PATIENT));
  }
  // Prepare charge densities for Fourier transform
  for (size_t i = 0; i < grid_count; ++i) {
    r = dr * static_cast<double>(i+1);
    for (size_t j = 0; j < species_count; ++j) {
      weighted_densities.at(j).at(0).at(i) = r *
          charge_density_profiles.at(j).at(i);
    }
  }
  // Execute Fourier transform
  for (auto& plan : plans_forward) {
    fftw_execute(plan);
  }
  // Prefactors (4\pi: spherical sym.; 1/2: fftw internal factor)
  // we dont divide by the radial wave number kr, because when we do the
  // backward transform we need to multiply with it again.
  for (size_t j = 0; j < species_count; ++j) {
    weighted_densities.at(j).at(0) *= 2. * M_PI * dr;
  }
  // Extend data
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 1; j < species_count; ++j) {  // j=1
      weighted_densities.at(i).at(j) = weighted_densities.at(i).at(0);
    }
  }
  // Convolution of charge densities and weights
  // Here we would multiply with kr again.
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 0; j < species_count; ++j) {
      weighted_densities.at(i).at(j) *= weights_delta.at(i).at(j);
    }
  }
  // Summation over all species of the weight functions before back transforming
  for (size_t i = 0; i < species_count; ++i) {
    for (size_t j = 1; j < species_count; ++j) {  // j=1
      weighted_densities.at(0).at(i) += weighted_densities.at(j).at(i);
    }
  }
  // Back transform into real space
  for (auto& plan : plans_backward) {
    fftw_execute(plan);
  }
  // Normalization
  for (size_t i = 0; i < grid_count; ++i) {
    r = dr * static_cast<double>(i+1);
    for (size_t j = 0; j < species_count; ++j) {
      weighted_densities.at(0).at(j).at(i) *= 4. * M_PI / r;
      weighted_densities.at(0).at(j).at(i) *= dkr / 2.;
      weighted_densities.at(0).at(j).at(i) /= pow(2. * M_PI, 3);
    }
  }
  for (auto& plan : plans_forward) { fftw_destroy_plan(plan); }
  for (auto& plan : plans_backward) { fftw_destroy_plan(plan); }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  size_t i{0};
  // Calculate the total charge density
  calc_charge_densities();
  // From the charge densities calculate the weighted densities
  calc_weighted_densities();
  calc_poisson_rhs();
  // From the charge densities calculate the electrostatic potential
  calc_potentials();
  // Calculate the derivative from the electrostatic potential
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    i = it - affected_species.begin();
    functional_derivative->at(*it) = valencies.at(i) * potentials.at(i);
  }
}
// _____________________________________________________________________________
void FunctionalESDeltaSpherical::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  size_t index{0};
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
  for (size_t i = 0; i < species_count; ++i) {
    index = affected_species.at(i);
    for (size_t j = 0; j < species_count; ++j) {
      bulk_derivative->at(index) += valencies.at(j) * bulk_densities.at(j) *
          pow(diameters.at(i) + diameters.at(j), 2);
    }
    bulk_derivative->at(index) *= -M_PI * valencies.at(i) * bjerrum / 6.;
  }
}
// _____________________________________________________________________________
double FunctionalESDeltaSpherical::calc_energy() {
  DataFrame<1, double> energy_density(grid_count);
  double integral{0.};
  calc_charge_densities();
  calc_weighted_densities();
  calc_potentials();
  for (size_t i = 0; i < species_count; ++i) {
    energy_density += .5 * potentials.at(i) * charge_density_profiles.at(i);
  }
  integral = integration_1d_radial_open_closed(energy_density, dr);
  return integral;
}
// _____________________________________________________________________________
