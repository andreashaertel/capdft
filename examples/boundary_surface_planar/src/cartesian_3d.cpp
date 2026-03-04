// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/boundary_surface_planar/src/cartesian_3d.cpp
 *  \brief Main file of the example of the planar functionals in 3d.
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
#include "../../../src/functional_fmt_cartesian.hpp"
#include "../../../src/functional_es_mf_cartesian.hpp"
#include "../../../src/cartesian_poisson_solver.hpp"
#include "../../../src/iterator.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/stl_algorithms.hpp"
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
  System<3> system(&params);
  // Get some properties from System
  std::vector<size_t> affected_species_fmt = system.affected_species_fmt;
  std::vector<size_t> affected_species_es = system.affected_species_es;
  std::vector<Properties> species_properties = system.species_properties;
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> system_lengths = system.system_lengths;
  std::vector<bool> periodic_boundaries = system.periodic_boundaries;
  double potential;
  system.get_property("potential", &potential);
  potential *= 1e4 * ELECTRON_CHARGE / (BOLTZMANN * system.temperature);
// _____________________________________________________________________________
  /* The functional objects are created. As input they require the Properties of
   * the system and the std::vector<Properties> from the species. Moreover
   * the pointer of the three-dimensional density profiles needs to be passed
   * on, as well as a std::vector<size_t> with the species numbers, that are
   * affected by the functional.
   * 
   * Remark: Even though the DataFrames (e.g. density_profiles) are
   * three-dimensional, the underlying data arrays are still one-dimensional.
   * The only difference are the 3-dim. helper functions like
   * at(size_t,size_t,size_t).
   */
// _____________________________________________________________________________
  // Create density profiles
  std::vector<DataFrame<3, double>> density_profiles(0);
  for (size_t i = 0; i < species_properties.size(); ++i) {
    density_profiles.push_back(DataFrame<3, double>(grid_counts));
  }
  // Create an FMT Functional object. For this we specify the
  // species which are interacting via this functional in affected_species.
  FunctionalFMTCartesian my_fmt_functional(&density_profiles, system);
  // Create an ES functional object.
  CartesianPoissonSolver poisson_solver(grid_counts, system.bin_sizes,
		  periodic_boundaries); // default boundary values: zero
  FunctionalESMFCartesian my_es_functional(
      &density_profiles, system, &poisson_solver);
// _____________________________________________________________________________
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functionals "my_fmt_functional" and
   * "my_es_functional" to the functional list and use the run_picard() function
   * to carry out the Picard iterations. In this case we use two planar hard
   * walls at distance system_length as external potential.
   *
   * If you want to save loads of time use the Anderson mixing algorithm
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
  double diameter{0.}, valency{0.}, bulk_density{0.}, z{0.};
  double value;
  double dz = system.bin_sizes.at(2);
  std::vector<DataFrame<3, double>> exp_ext_potential(
      species_properties.size(), DataFrame<3, double>(grid_counts));
  for (auto& potential : exp_ext_potential) {
    potential.set_all_elements_to(1.);
  }
  // Set external hard potential
  for (auto& species : affected_species_fmt) {
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t k = 0; k < grid_counts.at(2); ++k) {
      z = dz * (static_cast<double>(k) + 0.5);
      if (z < (diameter / 2.)) {
        value = 0.;
      } else if ((system_lengths.at(2) - z) < (diameter / 2.)) {
        value = 0.;
      } else {
        value = 1.;
      }
      for (size_t i = 0; i < grid_counts.at(0); ++i) {
        for (size_t j = 0; j < grid_counts.at(1); ++j) {
	  exp_ext_potential.at(species).at(i,j,k) = value;
	}
      }
    }
  }
  // Initial guess for the density profiles
  for (size_t s = 0; s < species_properties.size(); ++s) {
    species_properties.at(s).get_property("bulk density", &bulk_density);
    density_profiles.at(s).set_all_elements_to(bulk_density);
    density_profiles.at(s) *= exp_ext_potential.at(s);
  }
  // Set external electrostatic potential (linear between two charged plates,
  // with maximum / minimum value +/-potential)
  for (auto& species : affected_species_es) {
    species_properties.at(species).get_property("valency", &valency);
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t k = 0; k < grid_counts.at(2); ++k) {
      z = dz * (static_cast<double>(k) + 0.5);
      value = exp(- valency * potential * (1. - 2. * z / system_lengths.at(2)));
      for (size_t i = 0; i < grid_counts.at(0); i++) {
	for (size_t j = 0; j < grid_counts.at(1); j++) {
	  exp_ext_potential.at(species).at(i,j,k) *= value;
	}
      }
    }
  }
  std::fstream file;
  file.open("exp_ext_potential.dat", std::ios::out);
  file << "# [x] [y] [z] [exponentiated external potential profiles]\n";
  system.print_data(exp_ext_potential, file);
  file.close();
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential,
      species_properties);
  my_iterator.add_excess_functional(&my_fmt_functional);
  my_iterator.add_excess_functional(&my_es_functional);
  my_iterator.clear_convergence_criteria();
  my_iterator.add_convergence_criterion<ConvergenceCriterionSteps>(2e3);
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-6);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  double mixing = params.get_double("mixing", 0.1);
  //my_iterator.run_picard(mixing);
  my_iterator.run_anderson(mixing, 20);
// _____________________________________________________________________________
  /* All done!
   * Now we produce some output and view it in gnuplot.
   */
// _____________________________________________________________________________
  // Write density profiles to file
  double x{0.}, y{0.};
  double dx = system.bin_sizes.at(0);
  double dy = system.bin_sizes.at(1);
  std::fstream out_stream;
  out_stream.open("3d_profiles.dat", std::ios::out);
  out_stream << "# [x] [y] [z] [density profiles]" << std::endl;
  system.print_data(density_profiles, out_stream);
  out_stream.close();
  return 0;
}
