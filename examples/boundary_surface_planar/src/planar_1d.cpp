// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/boundary_surface_planar/src/planar_1d.cpp
 *  \brief Main file of the example of the planar functionals in 1d.
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/constants.hpp"
#include "../../../src/convergence_criterion.hpp"
#include "../../../src/convergence_criterion_max_dev.hpp"
#include "../../../src/convergence_criterion_steps.hpp"
#include "../../../src/convergence_criterion_nan.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_planar.hpp"
#include "../../../src/functional_es_mf_planar.hpp"
#include "../../../src/functional_es_delta_planar.hpp"
#include "../../../src/iterator.hpp"
#include "../../../src/properties.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  // Set the desired system properties
  /* The necessary geometric and physical properties are extracted from the
   * command line arguments (or a parameter file specified in the command line)
   * via the ParameterHandler. From that, a System object is constructed that
   * handles all those properties from then on.
   *
   * The System class calculates some additional properties (bin sizes, 
   * dielectric constant), sorts all species properties into separate Properties
   * objects, and determines which species are affected by the electrostatic 
   * and / or the FMT functional depending on the species' properties.
   */
// _____________________________________________________________________________
  ParameterHandler params(argc, args);
  params.process_parameters();
  System<1> system(&params);
  // Get some properties from System
  std::vector<size_t> affected_species_fmt = system.affected_species_fmt;
  std::vector<size_t> affected_species_es = system.affected_species_es;
  std::vector<Properties> species_properties = system.species_properties;
  size_t grid_count = system.grid_counts.at(0);
  double system_length = system.system_lengths.at(0);
  double potential;
  system.get_property("potential", &potential);
  std::cout << "potential: " << potential << " [V]";
  potential *= 1e4 * ELECTRON_CHARGE / (BOLTZMANN * system.temperature);
  std::cout << " = " << potential << " [reduced units]\n";
// _____________________________________________________________________________
  /* The functional objects are created. As input they require the Properties of
   * the system and the std::vector<Properties> from the species. Moreover
   * the pointer of the density profiles needs to be passed on, as well as a
   * std::vector<size_t> with the species numbers, that are affected by the
   * functional.
   */
// _____________________________________________________________________________
  // Create density profiles
  std::vector<DataFrame<1, double>> density_profiles(0);
  for (size_t i = 0; i < species_properties.size(); ++i) {
    density_profiles.push_back(DataFrame<1, double>(grid_count));
  }
  // Create an FMT Functional object. For this we specify the
  // species which are interacting via this functional in affected_species.
  FunctionalFMTPlanar my_fmt_functional(&density_profiles,
      species_properties, system, affected_species_fmt);
  // Create an ES functional object.
  FunctionalESMFPlanar my_es_functional(&density_profiles,
      species_properties, system, affected_species_es);
  //FunctionalESDeltaPlanar my_es_functional(&density_profiles,
  //    species_properties, system, affected_species_es);
// _____________________________________________________________________________
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functionals "my_fmt_functional" and
   * "my_es_functional" to the functional list and use the run_picard() function
   * to carry out the Picard iterations. In this case we use two planar hard
   * walls at distance system_length as external potential.
   *
   * If you want to save loads of time use the Anderson mixing algorithm
   * run_anderson() instead of the Picard iterations. They usually are faster
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
  double diameter{0.}, valency{0.}, bulk_density{0.}, z{0.};
  double dz{system_length / static_cast<double>(grid_count)};
  std::vector<DataFrame<1, double>> exp_ext_potential(
      species_properties.size(), DataFrame<1, double>(grid_count));
  for (auto& potential : exp_ext_potential) {
    potential.set_all_elements_to(1.);
  }
  // Set external hard potential
  for (auto& species : affected_species_fmt) {
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t j = 0; j != grid_count; ++j) {
      z = dz * (static_cast<double>(j) + 0.5);
      if (z < (diameter / 2.)) {
        exp_ext_potential.at(species).at(j) = 0.;
      } else if ((system_length - z) < (diameter / 2.)) {
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
  // Set external electrostatic potential (linear between two charged plates,
  // with maximum / minimum value +/-potential)
  for (auto& species : affected_species_es) {
    species_properties.at(species).get_property("valency", &valency);
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t j = 0; j != grid_count; ++j) {
      z = dz * (static_cast<double>(j) + 0.5);
      exp_ext_potential.at(species).at(j) *=
          exp(- valency * potential * (1. - 2. * z / system_length));
    }
  }
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential,
      species_properties);
  my_iterator.add_excess_functional(&my_fmt_functional);
  my_iterator.add_excess_functional(&my_es_functional);
  my_iterator.clear_convergence_criteria();
  my_iterator.add_convergence_criterion<ConvergenceCriterionSteps>(3e3);
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-5);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  //my_iterator.run_picard(1e-7);
  my_iterator.run_anderson(1e-4, 15);
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
  out_stream.open("planar_profile.dat", std::ios::out);
  for (size_t i = 0; i < grid_count; ++i) {
    z = dz * static_cast<double>(i);
    out_stream << z << " ";
    for (size_t j = 0; j < species_properties.size(); ++j) {
      out_stream << density_profiles.at(j).at(i) << " ";
    }
    out_stream << std::endl;
  }
  out_stream.close();
  // Obtain grand potential of the system
  double energy;
  energy = my_fmt_functional.calc_energy();
  std::cout << "Excess free energy per square nanometer of FMT functional: ";
  std::cout << energy << std::endl;
  energy = my_es_functional.calc_energy();
  std::cout << "Excess free energy of mean-field electrostatic functional: ";
  std::cout << energy << std::endl;
  return 0;
}
