// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_es_mf_spherical.hpp"  // NOLINT
#include <cmath>
// Define some natural constants
#define ELECTRON_CHARGE 1.602176634  // *1e-19
#define BOLTZMANN 1.38064852  // *1e-23
#define VACUUM_PERMITIVITY 8.8541878128  // *1e-12
#define AVOGADRO 6.02214076  // *1e23
// _____________________________________________________________________________
FunctionalESMFSpherical::FunctionalESMFSpherical() {
}
// _____________________________________________________________________________
FunctionalESMFSpherical::FunctionalESMFSpherical(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties,
    std::vector<size_t> affected_species)
  : affected_species(affected_species),
    density_profiles_pointer(density_profiles) {
  // Get system properties
  extract_system_properties(system_properties);
  // Get species properties; excludes all species without diameter property
  extract_species_properties(species_properties);
  // Initialize all data frames and update charge densities
  initialize_all_data_frames();
  calc_charge_densities();
  initialize_poisson_solver();
}
// _____________________________________________________________________________
FunctionalESMFSpherical::FunctionalESMFSpherical(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties)
  : FunctionalESMFSpherical(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalESMFSpherical::~FunctionalESMFSpherical() {
  //
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::extract_system_properties(
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
void FunctionalESMFSpherical::extract_electrical_properties(
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
      std::cerr << "FunctionalESMFSpherical::extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "FunctionalESMFSpherical::extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::extract_species_properties(
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
void FunctionalESMFSpherical::initialize_all_data_frames() {
  charge_density_profile.~DataFrame<1, double>();
  new(&charge_density_profile) DataFrame<1, double>(grid_count);
  poisson_rhs.~DataFrame<1, double>();
  new(&poisson_rhs) DataFrame<1, double>(grid_count);
  potential.~DataFrame<1, double>();
  new(&potential) DataFrame<1, double>(grid_count);
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::initialize_poisson_solver() {
  poisson_solver = new RadialPoissonSolver(grid_count, dr, .5 * dr);
  poisson_solver->set_radial_laplacian(NEUMANN_DIRICHLET);
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_charge_densities() {
  charge_density_profile.zero();
  for (size_t i = 0; i < species_count; ++i) {
    charge_density_profile += valencies.at(i) *
        density_profiles_pointer->at(affected_species.at(i));
  }
  poisson_rhs = -4. * M_PI * bjerrum * charge_density_profile;
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_potential() {
  double outer_boundary{0.};
  double inner_boundary{0.};
  // Calculate the boundary conditions of the Poisson equation.
  // The inner one equals 0 (Neumann) due to the radial symmetry and the outer
  // one equals net charge divided by radial position (Dirichlet) due to Gauss'
  // theorem. Note, that there is no external charge at the center.
  inner_boundary = 0.;
  outer_boundary = bjerrum * calc_net_charge() /
      (dr * static_cast<double>(grid_count+1));
  // Solve the radial Poisson equation numerically
  poisson_solver->solve(
      inner_boundary, outer_boundary, poisson_rhs.array(), potential.array());
}
// _____________________________________________________________________________
double FunctionalESMFSpherical::calc_net_charge() {
  double integral{0.};
  double r{0.};
  // Trapezoidal integral rule, spherical integral
  r = static_cast<double>(0 + 1) * dr;
  integral += .5 * r * r * charge_density_profile.at(0) * dr;
  for (size_t i = 1; i < grid_count - 1; ++i) {
    r = static_cast<double>(i + 1) * dr;
    integral += r * r * charge_density_profile.at(i) * dr;
  }
  integral += .5 * r * r * dr * charge_density_profile.at(grid_count - 1);
  // Spherical symmetry
  integral *= 4. * M_PI;
  return integral;
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  size_t i{0};
  // Calculate the total charge density
  calc_charge_densities();
  // From the charge densities calculate the electrostatic potential
  calc_potential();
  // Calculate the derivative from the electrostatic potential
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    i = it - affected_species.begin();
    functional_derivative->at(*it) = valencies.at(i) * potential;
  }
}
// _____________________________________________________________________________
void FunctionalESMFSpherical::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // The bulk values of this functional's derivative are always zero
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
}
// _____________________________________________________________________________
double FunctionalESMFSpherical::calc_energy() {
  DataFrame<1, double> integrand(grid_count);
  double integral{0.};
  double r{0.};
  calc_charge_densities();
  calc_potential();
  integrand = potential * charge_density_profile;
  // Trapezoidal integral rule, spherical integral
  r = static_cast<double>(0 + 1) * dr;
  integral += .5 * r * r * integrand.at(0) * dr;
  for (size_t i = 1; i < grid_count - 1; ++i) {
    r = static_cast<double>(i + 1) * dr;
    integral += r * r * integrand.at(i) * dr;
  }
  integral += .5 * r * r * dr * integrand.at(grid_count - 1);
  // Spherical symmetry
  integral *= 4. * M_PI;
  // 1/2 occuring in the MF functional
  integral /= 2.;
  return integral;
}
// _____________________________________________________________________________
