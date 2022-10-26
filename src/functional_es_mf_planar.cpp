// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_es_mf_planar.hpp"  // NOLINT
#include <cmath>
#include <algorithm>
#include "constants.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "integration.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalESMFPlanar::FunctionalESMFPlanar() {
}
// _____________________________________________________________________________
FunctionalESMFPlanar::FunctionalESMFPlanar(
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
  initialize_all_data_frames();
  calc_charge_densities();
  initialize_poisson_solver();
}
// _____________________________________________________________________________
FunctionalESMFPlanar::FunctionalESMFPlanar(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties)
  : FunctionalESMFPlanar(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalESMFPlanar::~FunctionalESMFPlanar() {
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::extract_system_properties(
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
void FunctionalESMFPlanar::extract_electrical_properties(
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
      std::cerr << "FunctionalESMFPlanar::extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "FunctionalESMFPlanar::extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::extract_species_properties(
    const std::vector<Properties>& species_properties) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  double valency{0.};
  if (affected_species.empty()) {
    for (auto& properties : species_properties) {
      try {
        if (properties.get_property("valency", &valency)) {
          affected_species.push_back(&properties - &species_properties[0]);
        }
      } catch(...) {}
    }
  }
  // Get valencies
  for (auto& species : affected_species) {
    species_properties.at(species).get_property("valency", &valency);
    valencies.push_back(valency);
  }
  // Count species that interact via the electrostatic forces
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::initialize_all_data_frames() {
  // The initialization of the variables (header) called the empty constructor.
  // The following construction initializes the DataFrames correctly.
  charge_density_profile.~DataFrame<1, double>();
  new(&charge_density_profile) DataFrame<1, double>(grid_count);
  poisson_rhs.~DataFrame<1, double>();
  new(&poisson_rhs) DataFrame<1, double>(grid_count);
  potential.~DataFrame<1, double>();
  new(&potential) DataFrame<1, double>(grid_count);
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::initialize_poisson_solver() {
  poisson_solver = new PlanarPoissonSolver(grid_count, dz, DIRICHLET_DIRICHLET);
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::calc_charge_densities() {
  charge_density_profile.zero();
  for (size_t i = 0; i < species_count; ++i) {
    charge_density_profile += valencies.at(i) *
        density_profiles_pointer->at(affected_species.at(i));
  }
  poisson_rhs = -4. * M_PI * bjerrum * charge_density_profile;
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::calc_potential() {
  // Both boundary conditions equal 0 (Dirichlet). Note, that there is no
  // external potential involved inside the functional.
  double left_boundary{0.};
  double right_boundary{0.};
  // Solve the Poisson equation numerically
  poisson_solver->solve(
      left_boundary, right_boundary, poisson_rhs.array(), potential.array());
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  size_t i{0};
  // Calculate the total charge density
  calc_charge_densities();
  // From the charge densities calculate the electrostatic potential
  calc_potential();
  // Calculate the functional derivative which is the electrostatic potential
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    i = it - affected_species.begin();
    functional_derivative->at(*it) = valencies.at(i) * potential;
  }
}
// _____________________________________________________________________________
void FunctionalESMFPlanar::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // The bulk values of this functional's derivative are always zero
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
}
// _____________________________________________________________________________
double FunctionalESMFPlanar::calc_energy() {
  DataFrame<1, double> energy_density(grid_count);
  double integral{0.};
  calc_charge_densities();
  calc_potential();
  energy_density = .5 * potential * charge_density_profile;
  integral = integration_1d_closed(energy_density, dz);
  return integral;
}
// _____________________________________________________________________________
