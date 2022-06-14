// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "iterator.hpp"  // NOLINT
#include <cmath>
#include "data_frame.hpp"  // NOLINT
// _____________________________________________________________________________
Iterator::Iterator(
    std::vector<DataFrame<1, double>>* density_profiles,
    const std::vector<DataFrame<1, double>>& exp_external_potentials,
    const std::vector<Properties>& species_properties)
  : density_profiles(density_profiles),
    exp_external_potentials(&exp_external_potentials),
    species_properties(&species_properties) {
}
// _____________________________________________________________________________
Iterator::~Iterator() {
}
// _____________________________________________________________________________
void Iterator::add_excess_functional(Functional& functional) {
  size_t species_count{species_properties->size()};
  size_t index{0};
  excess_functionals.push_back(&functional);
  // For every functional and every species a DataFrame is created
  bulk_derivatives.push_back(std::vector<double>(0));
  functional_derivatives.push_back(std::vector<DataFrame<1, double>>(0));
  index = functional_derivatives.size() - 1;
  for (size_t i = 0; i < species_count; ++i) {
    bulk_derivatives.at(index).push_back(0.);
    functional_derivatives.at(index).push_back(DataFrame<1, double>(
        density_profiles->at(0)));  // copy constructor
  }
}
// _____________________________________________________________________________
void Iterator::remove_excess_functional(size_t index) {
  excess_functionals.erase(excess_functionals.begin() + index);
  functional_derivatives.erase(functional_derivatives.begin() + index);
  bulk_derivatives.erase(bulk_derivatives.begin() + index);
}
// _____________________________________________________________________________
double Iterator::run() {
  DataFrame<1, double> proposed_density;
  DataFrame<1, double> deviations;
  std::vector<double> maximum_deviations(0);
  double maximum_deviation{std::numeric_limits<double>::max()};
  double target_deviation{1.0e-6};
  double mixing{.05};
  double bulk_density{0.};
  size_t steps{0};
  size_t species_count{species_properties->size()};
  // Calculate bulk derivatives
  for (size_t i = 0; i < excess_functionals.size(); ++i) {
    excess_functionals.at(i)->calc_bulk_derivative(&bulk_derivatives.at(i));
  }
  // Iterator loop
  while (maximum_deviation > target_deviation) {
    ++steps;
    // Calculate the functional derivatives.
    for (size_t i = 0; i < excess_functionals.size(); ++i) {
      excess_functionals.at(i)->calc_derivative(&functional_derivatives.at(i));
    }
    // Picard iteration scheme
    for (size_t i = 0; i < species_count; ++i) {
      species_properties->at(i).get_property("bulk density", &bulk_density);
      // Calculate the right hand side of the update formula of cDFT.
      // No iteration over the grid is needed, because of the overloaded
      // operators of DFSpherical.
      proposed_density = bulk_density * exp_external_potentials->at(i);
      for (size_t j = 0; j < functional_derivatives.size(); ++j) {
        proposed_density *= exp(bulk_derivatives.at(j).at(i)) *
            exp(-1. * functional_derivatives.at(j).at(i));
      }
      // Mix the new solution with the old one
      density_profiles->at(i) = density_profiles->at(i) * (1. - mixing) +
          mixing * proposed_density;
      deviations = abs(density_profiles->at(i) - proposed_density);
      maximum_deviations.push_back(max(deviations));
    }
    maximum_deviation = *std::max_element(
        maximum_deviations.begin(), maximum_deviations.end());
    maximum_deviations.clear();
    // Print the deviation and the number of iteration steps taken
    std::cout << "Picard iteration step " << steps;
    std::cout << ". Deviation: " << maximum_deviation << std::endl;
  }
  return maximum_deviation;
}
// _____________________________________________________________________________
double Iterator::calculate_gc_energy() {
  //
  return 0.0;
}
// _____________________________________________________________________________
