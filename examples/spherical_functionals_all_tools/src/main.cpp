// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
 *  \brief Main file of the example of the spherical functionals.
 *  
 *  This main file contains examples to show how Functional classes work.
 *
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/convergence_criterion.hpp"
#include "../../../src/convergence_criterion_max_dev.hpp"
#include "../../../src/convergence_criterion_steps.hpp"
#include "../../../src/convergence_criterion_nan.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_spherical.hpp"
#include "../../../src/functional_es_mf_spherical.hpp"
#include "../../../src/functional_es_delta_spherical.hpp"
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
   *
   * Note, that the system_length is chosen such that the external potential
   * step lies exactly between two bins.
   */
// _____________________________________________________________________________
  size_t grid_count = static_cast<size_t>(1e4+.5) + 1;  // equals = 10,001
  double ext_potential_charge = 1.;
  double system_length = 19.821782178217823;  // in nm
  double bjerrum_length = 1.;  // in nm
  double temperature = 300.;  // in K
  // Create objects of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length", system_length);
  system_properties.add_property<double>("bjerrum length", bjerrum_length);
  system_properties.add_property<double>("temperature", temperature);
  system_properties.add_property<size_t>("grid count", grid_count);
  // First species
  properties.add_property<double>("diameter", .3);
  properties.add_property<double>("bulk density", 3.);
  properties.add_property<double>("valency", -1.);
  species_properties.push_back(properties);
  properties.clear();
  // Second species
  properties.add_property<double>("bulk density", 1.5);
  species_properties.push_back(properties);
  properties.clear();
  // Third species
  properties.add_property<double>("diameter", .3);
  properties.add_property<double>("bulk density", 3.);
  properties.add_property<double>("valency", +1.);
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
  // Create an FMT Functional object. For this we specify the
  // species which are interacting via this functional in affected_species.
  std::vector<size_t> affected_species_fmt{0, 2};  // selected species 0 and 2
  FunctionalFMTSpherical my_fmt_functional(&density_profiles,
      species_properties, system_properties, affected_species_fmt);
  // Create an ES functional object.
  std::vector<size_t> affected_species_es{0, 2};  // selected species 0 and 2
  //FunctionalESMFSpherical my_es_functional(&density_profiles,
  //    species_properties, system_properties, affected_species_es);
  FunctionalESDeltaSpherical my_es_functional(&density_profiles,
      species_properties, system_properties, affected_species_es);
// _____________________________________________________________________________
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functionals "my_fmt_functional" and
   * "my_es_functional" to the functional list and use the run_picard() function
   * to carry out the Picard iterations. In this case we use two spherical hard
   * walls as external potential.
   *
   * If you want to save loads of time use the Andersen mixing algorithm
   * run_andersen() instead of the Picard iterations. They usually are faster
   * by a factor of 20.
   *
   * The convergence criteria are also added to the Iterator. They are added via
   * a template argument that specifies the type of criterion and a function
   * argument that specifies the threshold of the specified criterion.
   *
   * For example: 
   * my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-4);
   * This will terminate the iterations after the largest difference
   * (ConvergenceCriterionMaxDev) of the old density profile and the new one is
   * smaller than 1.0e-4.
   */
// _____________________________________________________________________________
  // Create external potential DataFrames
  double diameter{0.}, valency{0.}, bulk_density{0.}, r{0.};
  double dr{system_length / static_cast<double>(grid_count)};
  std::vector<DataFrame<1, double>> exp_ext_potential(
      species_properties.size(), DataFrame<1, double>(grid_count));
  for (auto& potential : exp_ext_potential) {
    potential.set_all_elements_to(1.);
  }
  // Set external hard potential
  for (auto& species : affected_species_fmt) {
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t j = 0; j != grid_count; ++j) {
      r = dr * static_cast<double>(j + 1);
      if (r < diameter) {
        exp_ext_potential.at(species).at(j) = 0.;
      } else if (system_length - r < 1.5 * diameter) {
        exp_ext_potential.at(species).at(j) = 0.;
      } else {
        exp_ext_potential.at(species).at(j) = 1.;
      }
    }
  }
  // Initial guess for the density profiles
  for (size_t i = 0; i < species_properties.size(); ++i) {
    species_properties.at(i).get_property("bulk density", &bulk_density);
    density_profiles.at(i).set_all_elements_to(bulk_density);
    density_profiles.at(i) *= exp_ext_potential.at(i);
  }
  // Set external electrostatic potential
  for (auto& species : affected_species_es) {
    species_properties.at(species).get_property("valency", &valency);
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t j = 0; j != grid_count; ++j) {
      r = dr * static_cast<double>(j + 1);
      if (r < diameter) { continue; }  // prevents divergences in the walls
      exp_ext_potential.at(species).at(j) *=
          exp(-bjerrum_length * ext_potential_charge * valency / r);
    }
  }
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential,
      species_properties);
  my_iterator.add_excess_functional(&my_fmt_functional);
  my_iterator.add_excess_functional(&my_es_functional);
  my_iterator.clear_convergence_criteria();
  my_iterator.add_convergence_criterion<ConvergenceCriterionSteps>(2e3);
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-4);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  my_iterator.run_picard(1.5e-4);
  //my_iterator.run_andersen(1.5e-4, 10);
// _____________________________________________________________________________
  /* All done!
   * Now we produce some output and view it in gnuplot.
   * We also supplied this example with a pdf that shows the plot in case you do
   * not use gnuplot.
   *
   * You can also try to reproduce the density profiles of
   * [M. Bültmann and A. Härtel 2022 J. Phys.: Condens. Matter 34 235101].
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
  energy = my_es_functional.calc_energy();
  std::cout << "Excess free energy of mean-field electrostatic functional: ";
  std::cout << energy << std::endl;
  return 0;
}
