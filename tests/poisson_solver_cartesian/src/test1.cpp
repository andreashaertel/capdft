// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/poisson_solver_cartesian/src/main.cpp
  * \brief Test the CartesianPoissonSolverAny for the simple case of planar
  * boundaries and a homogeneous charge distribution; compare with analytic
  * expectation
  *
  * An example parameter file is provided (../params1.txt). You can also run
  * this test with different parameter settings, but there are some restrictions
  * (e.g. this code assumes that the two planar boundary surfaces are orthogonal
  * to the z axis).
  */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundaries.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/cartesian_poisson_solver_any.hpp"
#include "../../../src/system.hpp"
#include "../../../src/discretization_handler.hpp"
#include <vector>
#include <string> // for parameter labels
#include <cmath> // std::exp
#include <math.h> // M_PI
#include <iostream>
#include <fstream> // write results to file
// _____________________________________________________________________________
// Main function
int main(int argc, char** argv) {
// _____________________________________________________________________________
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " -f <parameter file>\n";
    std::cout << "(or use flag -p <parameter>=<value>, cf. ParameterHandler)\n";
    exit(1);
  }
  std::string log_text;
// _____________________________________________________________________________
  /* Define System  
   * The required parameters are extracted from the input via ParameterHandler.
   */
// _____________________________________________________________________________
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
  // Define System
  System<Cartesian> system(&params);
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> lengths = system.system_lengths;
  std::vector<double> bin_sizes = system.bin_sizes;
// _____________________________________________________________________________
  /* Define the right hand side of the Poisson equation     
   * Here, we choose a homogeneous density distribution for the ions, equal to
   * the bulk density given in the input parameters. For simplicity, we only 
   * consider one ion species.
   */
// _____________________________________________________________________________
  // Species properties
  size_t species = 0;
  double density;
  double valency;
  // Density profiles
  system.species_properties.at(species).get_property("bulk density", &density);
  system.species_properties.at(species).get_property("valency", &valency);
  double charge_density = valency * density;
  // Poisson right hand side: homogeneous charge density 
  double rhs_value = - 4. * M_PI * system.bjerrum * charge_density;
  DataFrame<3, double> rhs(grid_counts);
  rhs.set_all_elements_to(rhs_value);
// _____________________________________________________________________________
  /* Define boundary conditions  
   * We set the boundary positions and boundary values via two
   * BoundarySurfacePlanar objects.
   */
// _____________________________________________________________________________
  Boundaries<Cartesian> boundaries;
  BoundarySurfacePlanar boundaries_left(system, 0);
  BoundarySurfacePlanar boundaries_right(system, 1);
  boundaries.add_surface(&boundaries_left);
  boundaries.add_surface(&boundaries_right);
// _____________________________________________________________________________
  /* Solve Poisson equation with CartesianPoissonSolverAny
   */
// _____________________________________________________________________________
  DataFrame<3, double> solution(grid_counts);
  CartesianPoissonSolverAny poisson(system, &boundaries);
  poisson.solve(rhs, solution);
// _____________________________________________________________________________
  /* Calculate expected potential profile analytically    
   * We assume that the boundary surfaces are positioned at z=0 and 
   * z=system_lengths.at(2). Then the electrostatic potential should
   * follow a quadratic shape Phi(z) = c2 z² + c1 z + c0, where the coefficients
   * c2,c1,c0 are determined by the boundary values and the charge density.
   */
// _____________________________________________________________________________
  double potential0 = boundaries_left.get_boundary_value();
  double potential1 = boundaries_right.get_boundary_value();
  double L_z = lengths.at(2) - bin_sizes.at(2);
  double z;
  // polynomial coefficients:
  double c2 = 0.5 * rhs_value;
  double c1 = (potential1 - potential0) / L_z - c2 * L_z;
  double c0 = potential0;
  // Calculate potential profile
  DataFrame<3, double> solution_expected(grid_counts);
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
	z = static_cast<double>(k) * bin_sizes.at(2);
	solution_expected.at(i,j,k) = c2 * z * z + c1 * z + c0; 
      }
    }
  }
// _____________________________________________________________________________
  /* Produce output
   * We write the calculated potential profiles into a data file.
   */
// _____________________________________________________________________________
  double max_dev = max(abs(solution - solution_expected));
  double rel_dev = max_dev / fabs(potential1 - potential0);
  std::cout << "Maximum deviation: " << max_dev << " = " << rel_dev * 100.;
  std::cout << "% of the applied voltage\n";
  std::fstream log;
  log.open("test1.log", std::ios::out);
  log << "Maximum deviation: " << max_dev << " = " << rel_dev * 100.;
  log << "% of the applied voltage\n";
  std::fstream ES_dat;
  std::string filename{"potential_profiles1.dat"};
  ES_dat.open(filename, std::ios::out);
  ES_dat << "# [x] [y] [z] [calculated] [expected]" << std::endl;
  // write grid dimensions to header
  ES_dat << "# grid_counts=";
  for (size_t g : grid_counts) {
    ES_dat << g << ",";
  }
  ES_dat << std::endl;
  // write data
  double x, y;
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    x = bin_sizes.at(0) * static_cast<double>(i);
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      y = bin_sizes.at(1) * static_cast<double>(j);
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
        z = bin_sizes.at(2) * static_cast<double>(k);
        ES_dat << x << " " << y << " " << z << " ";
        ES_dat << solution.at(i, j, k) << " " << solution_expected.at(i, j, k) << std::endl;
      }
    }
  }
  ES_dat.close();
// _____________________________________________________________________________
  std::cout << "Wrote calculated and expected potential profiles to "
	       "potential_profiles1.dat.\n";
  log << "Wrote calculated and expected potential profiles to "
	       "potential_profiles1.dat.\n";
  log.close();
}
