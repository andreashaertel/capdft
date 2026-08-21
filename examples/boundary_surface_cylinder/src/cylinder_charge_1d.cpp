// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/boundary_surface_cylinder/src/cylinder_charge_1d.cpp
 *  \brief Calculate the charge distribution on the cylinder surface
 *  along a few paths
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_cylinder.hpp"
#include "../../../src/surface_charge_distribution.hpp"
#include "../../../src/system.hpp"
#include "../../../src/integration.hpp"
#include "../../../src/stl_algorithms.hpp"
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
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " -f <parameter_file> ";
    std::cout << "-p potential_profile=<filename1> ";
    std::cout << "-p density_profiles=<filename2>\n";
    std::cout << "	where <filename1,2> are the data files of the total "
	 	 "electrostatic potential and the density profiles.\n";
    std::cout << "	If not "
	         "specified, the default filenames 'total_ES.dat' and "
		 "'3d_profiles.dat' are used.\n";
    exit(1);
  }
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
  ParameterHandler params(argc, argv);
  params.process_parameters();
  System<Cartesian> system(&params);
  // System properties:
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> bin_sizes = system.bin_sizes;
  std::vector<size_t> xy_counts{grid_counts.at(0), grid_counts.at(1)};
  double resolution = bin_sizes.at(0);
  std::vector<double> system_lengths = system.system_lengths;
  double potential = params.get_double("potential");
  std::cout << "potential = " << potential << " V "; // in Volt
  potential *= ELECTRON_CHARGE / (BOLTZMANN * system.temperature) * 1e4;
  std::cout << "= " << potential << " kT/e\n"; // in reduced units
// _____________________________________________________________________________
  /* Construction of the BoundarySurface objects
   * The surface geometry is specified in the constructor. The boundary values
   * of the potential are defined via the set_boundary_values method.
   */
// _____________________________________________________________________________
  BoundarySurfaceCylinder cylinder(system, 0);
  cylinder.set_boundary_value(potential);
  // extract cylinder radius and height for later usage:
  double half_height = params.get_double("0_height") * 0.5;
  double z_mid = params.get_double("0_z_mid",
		  		   system.system_lengths.at(2) * 0.5);
// _____________________________________________________________________________
  /* Extract the electrostatic potential profile from a data file
   * The filename is specified in the command line parameters.
   */
// _____________________________________________________________________________
  DataFrame<3,double> total_ES(grid_counts);
  std::string potential_file =
	  params.get_string("potential_profile", "total_ES.dat");
  size_t data_column = 3;
  system.load_data(potential_file, data_column, &total_ES);
// _____________________________________________________________________________
  /* Calculate charge distribution on the cylinder surface
   *
   * We first determine the total electrostatic potential profile in the system.
   * From that, the charge densities on the surface can be calculated via the
   * function charge_density (or alternatively charge_distribution). In order
   * to visualize the results in a plot, we split the surface into three regions
   * (cylinder mantle and the two cylinder caps) and store the results in three
   * different data files. Note that the data files also include points that do
   * not lie on the surface, where we set the value to zero.
   */
// _____________________________________________________________________________
