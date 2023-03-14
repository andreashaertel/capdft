// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functionals_some_tools/src/main.cpp
 *  \brief Main file of the example of the spherical functionals.
 *  
 *  This main file contains examples to show how the FunctionalFMTSpherical
 *  class works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_spherical.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  // Set the desired system properties
  /* We start by defining the necessary geometric and physical properties.
   * The general properties of a system are put into one Properties container.
   * The properties of every species are put in separate Properties container.
   * The containers containing the species properties are encapsulated by an
   * std::vector.
   *
   * We will define three particle species.
   * Two of them will interact via this functional. The other species has no
   * diameter and will be ignored.
   *
   * In general programs should obtain their parameters via the command line
   * or a parameter file (e.g. via the ParameterHandler).
   * However, for this example we hard-coded the numbers to reduce complexity.
   */
// _____________________________________________________________________________
  size_t grid_count = 10001;
  double system_length = 19.821782178217823;  // external pot. between two bins
  double bjerrum_length = 1.;  // not really needed
  // Create objects of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length", system_length);
  system_properties.add_property<double>("bjerrum length", bjerrum_length);
  system_properties.add_property<size_t>("grid count", grid_count);
  // First species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .1);
  properties.add_property<double>("valency", -1.);  // not really needed
  species_properties.push_back(properties);
  properties.clear();
  // Second species
  properties.add_property<double>("bulk density", .2);
  properties.add_property<double>("valency", +1.);  // not really needed
  species_properties.push_back(properties);
  properties.clear();
  // Third species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .3);
  species_properties.push_back(properties);
  properties.clear();
// _____________________________________________________________________________
  /* The functional objects are created. As input they require the Properties of
   * the system and the std::vector<Properties> from the species. Moreover
   * the pointer of the density profiles needs to be passed on, as well as a
   * std::vector<size_t> with the species number,
   * that are affected by the functional.
   */
// _____________________________________________________________________________
  // Create density profiles
  std::vector<DataFrame<1, double>> density_profiles(0);
  for (size_t i = 0; i < species_properties.size(); ++i) {
    density_profiles.push_back(DataFrame<1, double>(grid_count));
  }
  // Create an FMT Functional object. Either specify the species
  std::vector<size_t> affected_species{0, 2};  // selected species 0 and 2
  FunctionalFMTSpherical my_fmt_functional(&density_profiles,
      species_properties, system_properties, affected_species);
  // or let the constructor decide which species can interact via this
  // functional
  // FunctionalFMTSpherical my_fmt_functional(&density_profiles,
  //     species_properties, system_properties);
// _____________________________________________________________________________
  // Picard iteration setup
  /* Now that we processed all the necessary parameters, it is time to set up
   * everything that is required prior to the Picard iteration scheme.
   *
   * We first have to specify an external potential, since otherwise
   * our profiles would be completely flat (bulk case).
   *
   * We decided to mimic the shape of a hard spherical capacitor
   * (without any charges). This external potential can also be used to
   * calculate pair correlation functions via the Percus trick, if the inner
   * part and the outer part are far away from each other.
   * From a numerical point of view it is also more stable to limit the system
   * in such a way.
   *
   * We also have to make an initial guess that makes sense, i.e. one where the
   * density profiles vanish beyond the hard potential walls.
   * Moreover the functional derivative values for the bulk case are required.
   */
// _____________________________________________________________________________
  std::vector<double> maximum_deviations(0);
  double maximum_deviation{std::numeric_limits<double>::max()};
  double target_deviation{1.0e-6};
  double r{0.};
  double dr{system_length / static_cast<double>(grid_count)};
  double diameter{0.};
  double bulk_density{0.};
  double mixing{2.0e-1};
  size_t step{0};
  DataFrame<1, double> proposed_density(grid_count);
  DataFrame<1, double> deviations(grid_count);
  std::vector<double> bulk_derivatives(species_properties.size());
  std::vector<DataFrame<1, double>> fmt_derivatives;
  std::vector<DataFrame<1, double>> exp_ext_potential;
  for (size_t i = 0; i != species_properties.size(); ++i) {
    fmt_derivatives.push_back(DataFrame<1, double>(grid_count));
    exp_ext_potential.push_back(DataFrame<1, double>(grid_count));
  }
  // Set external (hard) potential
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    species_properties.at(*it).get_property("diameter", &diameter);
    for (size_t j = 0; j != grid_count; ++j) {
      r = dr * static_cast<double>(j + 1);
      if (r < diameter) {
        exp_ext_potential.at(*it).at(j) = 0.;
      } else if (system_length - r < 1.5 * diameter) {
        exp_ext_potential.at(*it).at(j) = 0.;
      } else {
        exp_ext_potential.at(*it).at(j) = 1.;
      }
    }
  }
  // Initial guess for the density profiles
  for (size_t i = 0; i < species_properties.size(); ++i) {
    species_properties.at(i).get_property("bulk density", &bulk_density);
    density_profiles.at(i).set_all_elements_to(bulk_density);
    density_profiles.at(i) *= exp_ext_potential.at(i);
  }
  // Get the bulk derivatives
  my_fmt_functional.calc_bulk_derivative(&bulk_derivatives);
// _____________________________________________________________________________
  // Picard iterations
  /* It is time to calculate a density profile by iteratively calculating the
   * functional derivatives.
   * 
   * Since we use DFSpherical as data containers, which have a lot of overloaded
   * operators, we can simply write down the well known Picard iteration formula
   * without iterating over the grid point indices.
   *
   */
// _____________________________________________________________________________
  while (maximum_deviation > target_deviation) {
    ++step;
    // Calculate the functional derivative; Suppress warnings about unphysical
    // values
    my_fmt_functional.calc_derivative(&fmt_derivatives);
    for (auto it = affected_species.begin(); it != affected_species.end();
        ++it) {
      species_properties.at(*it).get_property("bulk density", &bulk_density);
      // Calculate the right hand side of the update formula of cDFT.
      // No iteration over the grid is needed, because of the overloaded
      // operators of DFSpherical.
      proposed_density = bulk_density * exp_ext_potential.at(*it) *
          exp(bulk_derivatives.at(*it)) * exp(-1. * fmt_derivatives.at(*it));
      // Mix the new solution with the old one
      density_profiles.at(*it) = density_profiles.at(*it) * (1. - mixing) +
          mixing * proposed_density;
      // Calculate how much the new density deviates from the old one at most
      deviations = abs(density_profiles.at(*it) - proposed_density);
      maximum_deviations.push_back(max(deviations));
    }
    // Determine the largest deviation ovaer all species
    maximum_deviation = *std::max_element(
        maximum_deviations.begin(), maximum_deviations.end());
    maximum_deviations.clear();
    // Print the deviation and the number of iteration steps taken
    std::cout << "Picard iteration step " << step;
    std::cout << ". Deviation: " << maximum_deviation << std::endl;
  }
// _____________________________________________________________________________
  /* All done!
   * Now we produce some output and view it in gnuplot.
   * We also supplied this example with a pdf that shows the plot in case you do
   * not use gnuplot.
   */
// _____________________________________________________________________________
  // Write density profile to file
  std::fstream out_stream;
  out_stream.open("spherical_profile.dat", std::ios::out);
  for (size_t i = 0; i < grid_count; ++i) {
    r = dr * static_cast<double>(i+1);
    out_stream << r << " ";
    for (size_t j = 0; j < species_properties.size(); ++j) {
      out_stream << density_profiles.at(j).at(i) << " ";
    }
    out_stream << std::endl;
  }
  out_stream.close();
  // Obtain grand potential of the system
  double energy;
  energy = my_fmt_functional.calc_energy();
  std::cout << "Excess free energy of FMT functional: ";
  std::cout << energy << std::endl;
  return 0;
}
