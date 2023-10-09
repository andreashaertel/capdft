// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_es_mf_cartesian.hpp"  // NOLINT
#include <cmath>
#include <algorithm>
#include "constants.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "integration.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalESMFCartesian::FunctionalESMFCartesian() {
}
// _____________________________________________________________________________
FunctionalESMFCartesian::FunctionalESMFCartesian(
    std::vector<DataFrame<3, double>>* density_profiles,
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
FunctionalESMFCartesian::FunctionalESMFCartesian(
    std::vector<DataFrame<3, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties)
  : FunctionalESMFCartesian(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalESMFCartesian::~FunctionalESMFCartesian() {
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Internal index (*it is external index)
  size_t j{0};
  // Calculate the total charge density
  calc_charge_densities();
  // From the charge densities calculate the electrostatic potential
  calc_potential();
  // Calculate the functional derivative which is the electrostatic potential
  for (size_t i = 0; i < voxel_count; ++i) {
    for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
      j = it - affected_species.begin();
      functional_derivative->at(*it).at(i) = valencies.at(j) * potential.at(i);
    }
  }
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // The bulk values of this functional's derivative are always zero
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
}
// _____________________________________________________________________________
double FunctionalESMFCartesian::calc_energy() {
  DataFrame<1, double> energy_density(voxel_count);
  double integral{0.};
  calc_charge_densities();
  calc_potential();
  for (size_t i = 0; i < voxel_count; ++i) {
    energy_density.at(i) = .5 * potential.at(i) * charge_density_profile.at(i);
  }
  // TODO(Moritz): 3D integration
  //integral = integration_1d_closed(energy_density, dz);
  return integral;
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::extract_system_properties(
    const Properties& system_properties) {
  bool dummy{false};  // std::vector<bool> has a bug, thus we need a workaround
  // Extract properties
  lengths.resize(3);
  grid_counts.resize(3);
  periodic_boundaries.resize(3);
  system_properties.get_property("length x", &lengths.at(0));
  system_properties.get_property("length y", &lengths.at(1));
  system_properties.get_property("length z", &lengths.at(2));
  system_properties.get_property("grid count x", &grid_counts.at(0));
  system_properties.get_property("grid count y", &grid_counts.at(1));
  system_properties.get_property("grid count z", &grid_counts.at(2));
  system_properties.get_property("PBC x", &dummy);
  periodic_boundaries.at(0) = dummy;  // workaround
  system_properties.get_property("PBC y", &dummy);
  periodic_boundaries.at(1) = dummy;  // workaround
  system_properties.get_property("PBC z", &dummy);
  periodic_boundaries.at(2) = dummy;  // workaround
  voxel_count = grid_counts.at(0) * grid_counts.at(1) * grid_counts.at(2);
  // Calculate bin sizes
  for (size_t i = 0; i < lengths.size(); ++i) {
    bin_sizes.push_back(lengths.at(i) / static_cast<double>(grid_counts.at(i)));
  }
  // Calculate electrical properties of the system
  extract_electrical_properties(system_properties);
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::extract_electrical_properties(
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
      std::cerr << "FunctionalESMFCartesian::extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "FunctionalESMFCartesian::extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::extract_species_properties(
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
void FunctionalESMFCartesian::initialize_all_data_frames() {
  // The initialization of the variables (header) called the empty constructor.
  // The following construction initializes the std::vectors correctly.
  charge_density_profile.clear();
  charge_density_profile.resize(voxel_count, 0.);
  poisson_rhs.clear();
  poisson_rhs.resize(voxel_count, 0.);
  potential.clear();
  potential.resize(voxel_count, 0.);
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::initialize_poisson_solver() {
  poisson_solver = new CartesianPoissonSolver(
      grid_counts, bin_sizes, periodic_boundaries);
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::calc_charge_densities() {
  // Reset charge density array
  std::fill(charge_density_profile.begin(), charge_density_profile.end(), 0.);
  // Calculate charge density and the rhs of the Poisson equation.
  // The 3D DataFrame is converted into a long 1D array.
  for (size_t i = 0; i < voxel_count; ++i) {
    for (size_t j = 0; j < species_count; ++j) {
      charge_density_profile.at(i) += valencies.at(j) *
        density_profiles_pointer->at(affected_species.at(j)).element(i);
    }
    poisson_rhs.at(i) = -4. * M_PI * bjerrum * charge_density_profile.at(i);
  }
}
// _____________________________________________________________________________
void FunctionalESMFCartesian::calc_potential() {
  // All six sides have vanishing boundary conditions (Dirichlet),
  // since there is no external potential involved inside the functional.
  std::vector<std::vector<double>> boundary_values{
    {0., 0.}, {0., 0.}, {0., 0.}};
  // Solve the Poisson equation numerically
  poisson_solver->solve(poisson_rhs, boundary_values, potential);
}
// _____________________________________________________________________________
