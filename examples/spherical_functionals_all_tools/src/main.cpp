// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
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
#include "../../../src/convergence_criterion.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_spherical.hpp"
#include "../../../src/iterator.hpp"
#include "../../../src/properties.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  // Set the desired system properties
  /* We start by defining the necessary geometric and physical properties.
   * The general properties of a system are put into one Properties conainer.
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
  size_t grid_count = static_cast<size_t>(1e4+.5) + 1;  // equals = 10,001
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
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functional "my_fmt_functional" to the
   * functional list and use the run() function to carry out the Picard
   * iterations. In this case we use two spherical hard walls as external
   * potential.
   */
// _____________________________________________________________________________
  // Create external potential DataFrames
  double diameter{0.}, bulk_density{0.}, r{0.};
  double dr{system_length / static_cast<double>(grid_count)};
  std::vector<DataFrame<1, double>> exp_ext_potential;
  for (size_t i = 0; i != species_properties.size(); ++i) {
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
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential,
      species_properties);
  my_iterator.add_excess_functional(my_fmt_functional);
  my_iterator.run();
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
