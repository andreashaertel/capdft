// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/structured_functionals_some_tools/src/main.cpp
 *  \brief This is an example for the usage of the BoundarySurface class.
 *
 *  This example is based on ../../cartesian_functionals_all_tools/src/main.cpp
 *  and adapted to demonstrate how those sinusoidal surfaces can be handled via
 *  the BoundarySurface class.
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/boundary_surface_sine.hpp"
#include "../../../src/boundaries.hpp"
#include "../../../src/external_potential.hpp"
#include "../../../src/cartesian_poisson_solver_any.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/convergence_criterion.hpp"
#include "../../../src/convergence_criterion_max_dev.hpp"
#include "../../../src/convergence_criterion_steps.hpp"
#include "../../../src/convergence_criterion_nan.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_cartesian.hpp"
#include "../../../src/functional_es_mf_cartesian.hpp"
#include "../../../src/iterator.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/system.hpp"
#include "../../../src/stl_algorithms.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  // Set the desired system properties
  /* We start by defining the necessary geometric and physical properties.
   * All properties are put into one Properties container, which is then used to
   * initialize a System object.
   *
   * The System class calculates some additional properties (bin sizes, 
   * dielectric constant), sorts all species properties into separate Properties
   * objects, and determines which species are affected by the electrostatic 
   * and / or the FMT functional depending on the species' properties.
   *
   * We will define three particle species.
   * Two of them will interact via this functional. The other species has no
   * diameter or valency and will be ignored.
   *
   * In general programs should obtain their parameters via the command line
   * or a parameter file (e.g. via our ParameterHandler).
   * However, for this example we hard-coded the numbers to reduce complexity.
   *
   * Note, that the system_length is chosen such that the external potential
   * step lies exactly between two bins.
   */
// _____________________________________________________________________________
  std::vector<size_t> grid_counts{16, 16, 16}; //128, 128, 128};  // (x, y, z)
  std::vector<double> system_lengths{2., 2., 2.};  // in nm (x, y, z)
  std::vector<bool> periodic_boundaries{true, true, false};  // (x, y, z)
  double bjerrum_length = 1.;  // in nm
  double temperature = 300.;  // in K
  double voltage = 0.1;  // in Volt
  // Convert voltage to reduced units
  double potential =
      1e4 * voltage * ELECTRON_CHARGE / (BOLTZMANN * temperature);
  // Create objects of Properties class
  Properties system_properties;
  // Put the system properties into a Properties container
  for (size_t dir = 0; dir < 3; dir++) {
    system_properties.add_property<double>("length", system_lengths.at(0), dir);
    system_properties.add_property<size_t>("grid_count", grid_counts.at(0), dir);
    system_properties.add_property<bool>("PBC", periodic_boundaries.at(0), dir);
  }
  system_properties.add_property<double>("bjerrum length", bjerrum_length);
  system_properties.add_property<double>("temperature", temperature);
  // First species
  system_properties.add_property<double>("diameter", .3, 0);
  system_properties.add_property<double>("bulk density", 3., 0);
  system_properties.add_property<double>("valency", -1., 0);
  // Second species
  system_properties.add_property<double>("bulk density", 1.5, 1);
  // Third species
  system_properties.add_property<double>("diameter", .3, 2);
  system_properties.add_property<double>("bulk density", 3., 2);
  system_properties.add_property<double>("valency", +1., 2);
  // Surface geometry
  for (size_t side : {0, 1}) {
    system_properties.add_property<double>("amplitude", 0.5, side);
    system_properties.add_property<size_t>("maxima_count", 1, side);
    system_properties.add_property<bool>("side_inversed", (side == 1), side);
  }
  // Initialize System object
  System<3> system(system_properties);
  // Get some properties from System
  std::vector<size_t> affected_species_fmt = system.affected_species_fmt;
  std::vector<size_t> affected_species_es = system.affected_species_es;
  std::vector<Properties> species_properties = system.species_properties;
// _____________________________________________________________________________
  /* Construction of the BoundarySurface objects
   * The surface geometry is specified in the constructor. The boundary values
   * of the potential are defined via the set_boundary_values method. Using
   * this, the external potential can be calculated. Afterwards, the boundary
   * values are reset to zero as required for the calculation of the excess
   * potential by the functional_es_mf_cartesian.
   */
// _____________________________________________________________________________
  // Initialize boundaries
  std::cout << "Define system boundaries\n";
  Boundaries<3> surfaces;
  BoundarySurfaceSine surface_left(system, 0);
  BoundarySurfaceSine surface_right(system, 1);
  surface_left.set_boundary_value(potential);
  surface_right.set_boundary_value(-potential);
  surfaces.add_surface(&surface_left);
  surfaces.add_surface(&surface_right);
  // Calculate hard sphere external potential
  std::cout << "Calculate external potential\n";
  std::vector<DataFrame<3, double>> exp_ext_potential_hs(
      species_properties.size(), DataFrame<3, double>(grid_counts));
  surfaces.exp_external_potential_hs(system, &exp_ext_potential_hs);
  // Calculate external electrostatic potential
  DataFrame<3, double> ext_potential_es(grid_counts);
  CartesianPoissonSolverAny poisson_solver(system, &surfaces);
  DataFrame<3, double> rhs(grid_counts);
  rhs.set_all_elements_to(0.);
  poisson_solver.solve(rhs, ext_potential_es);
  // Calculate total external potential
  std::vector<DataFrame<3, double>> exp_ext_potential_total = 
	  exp_ext_potential_hs;
  double valency = 0.;
  for (size_t species : affected_species_es) {
    species_properties.at(species).get_property("valency", &valency);
    exp_ext_potential_total.at(species) *= exp(- valency * ext_potential_es);
  }
  // Reset boundary values
  surfaces.set_all_boundary_values(0.);
// _____________________________________________________________________________
  /* The functional objects are created. They extract all (system and species) 
   * properties they need from the System object. The electrostatics functional
   * receives an CartesianPoissonSolverAny object to perform its calculations
   * of the electrostatic potential profile. Moreover the pointer of the
   * three-dimensional density profiles needs to be passed on.
   * 
   * Remark: Even though the DataFrames (e.g. density_profiles) are
   * three-dimensional, the underlying data arrays are still one-dimensional.
   * The only difference are the 3-dim. helper functions like
   * at(size_t,size_t,size_t).
   */
// _____________________________________________________________________________
  std::cout << "Initialize functionals\n";
  // Create density profiles
  std::vector<DataFrame<3, double>> density_profiles(0);
  for (size_t i = 0; i < species_properties.size(); ++i) {
    density_profiles.push_back(DataFrame<3, double>(grid_counts));
  }
  // Create an FMT Functional object. For this we specify the
  // species which are interacting via this functional in affected_species.
  FunctionalFMTCartesian my_fmt_functional(&density_profiles, system);
  // Create an ES functional object.
  CartesianPoissonSolverAny poisson_solver_zero(system, &surfaces);
  FunctionalESMFCartesian my_es_functional(&density_profiles, system,
		  &poisson_solver_zero);
// _____________________________________________________________________________
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functionals "my_fmt_functional" and
   * "my_es_functional" to the functional list and use the run_picard() function
   * to carry out the Picard iterations.
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
  std::cout << "Calculate density profiles\n";
  double bulk_density{0.};
  // Initial guess for the density profiles
  for (size_t i = 0; i < species_properties.size(); ++i) {
    species_properties.at(i).get_property("bulk density", &bulk_density);
    density_profiles.at(i).set_all_elements_to(bulk_density);
    density_profiles.at(i) *= exp_ext_potential_hs.at(i);
  }
  std::fstream out;
  out.open("initial_densities.dat", std::ios::out);
  out << "# [x] [y] [z] [density profiles]\n";
  system.print_data(density_profiles, out);
  out.close();
  out.open("extpot.dat", std::ios::out);
  out << "# [x] [y] [z] [exp. extpot profiles]\n";
  system.print_data(exp_ext_potential_total, out);
  out.close();
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential_total,
      species_properties);
  my_iterator.add_excess_functional(&my_fmt_functional);
  my_iterator.add_excess_functional(&my_es_functional);
  my_iterator.clear_convergence_criteria();
  my_iterator.add_convergence_criterion<ConvergenceCriterionSteps>(2e3);
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-6);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  //my_iterator.run_picard(1e-1);
  my_iterator.run_anderson(1e-1, 20);
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
  // Write density profiles to file
  out.open("3d_profiles.dat", std::ios::out);
  out << "# [x] [y] [z] [density profiles]\n";
  system.print_data(density_profiles, out);
  out.close();
  return 0;
}
