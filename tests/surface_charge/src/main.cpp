// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/surface_charge/src/main.cpp
 * \brief This file contains code to test the function charge_distribution.
 * 
 * Calculate the electrode charge distribution in simple problems with known
 * solutions (planar electrodes, linear or constant potential...).
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/surface_charge_distribution.hpp"
#include "../../../src/system.hpp"
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
  /* Set system properties and define boundary surface
   */
// _____________________________________________________________________________
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
  System<3> system(&params);
  // System properties:
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> bin_sizes = system.bin_sizes;
  std::vector<size_t> xy_counts{grid_counts.at(0), grid_counts.at(1)};
  double resolution = bin_sizes.at(0);
  std::vector<double> system_lengths = system.system_lengths;
  // potential on the electrodes
  double potential = params.get_double("potential");
  double temperature;
  system.get_property<double>("temperature", &temperature);
  std::cout << "potential = " << potential << " V "; // in Volt
  potential *= ELECTRON_CHARGE / (BOLTZMANN * temperature) * 1e4;
  std::cout << "= " << potential << " kT/e\n"; // in reduced units
  // Define boundary surface (planar wall at z=0)
  BoundarySurfacePlanar<3> surface(system, 0);
  surface.set_boundary_value(potential);
// _____________________________________________________________________________
  /* Define the electrostatic potential profile
   * We set a linear profile, decreasing along the z direction from value
   * 'potential' to zero.
   */
// _____________________________________________________________________________
  std::cout << "define potential profile\n";
  DataFrame<3,double> total_ES(grid_counts);
  double potential_step =  // such that value at position gc-1 is -potential:
	  -2. * potential / static_cast<double>(grid_counts.at(2) - 1);
  std::cout << "potential step " << potential_step << std::endl;
  double value = potential - potential_step;
  // first potential step is added in for loop, yielding value(z=0) = potential
  for (size_t k = 0; k < grid_counts.at(2); ++k) {
    value += potential_step;
    for (size_t i = 0; i < grid_counts.at(0); ++i) {
      for (size_t j = 0; j < grid_counts.at(1); ++j) {
	total_ES.at(i,j,k) = value;
      }
    }
  }
// _____________________________________________________________________________
  /* Calculate the charge distribution (and total charge) on the surface.
   * We also calculate the expected value for comparison.
   */
// _____________________________________________________________________________
  std::cout << "calculate charge distr\n";
  // using the normal method (calculated from the potential gradient):
  std::vector<std::pair<std::vector<double>, double>> charge_normal(0);
  double total_charge_normal;
  charge_distribution(surface, system, total_ES,
		  &charge_normal, &total_charge_normal, resolution);
  // alternative method:
//  std::vector<std::pair<std::vector<double>, double>> charge_poisson(0);
//  double total_charge_poisson;
//  charge_distribution_poisson(surface, system, total_ES,
//		  &charge_poisson, &total_charge_poisson, resolution);
  // expected result (total charge):
  double diff_expected =
	  2. * potential * system_lengths.at(0) * system_lengths.at(1) /
	  (system.bjerrum * 4. * M_PI * system_lengths.at(2));
// _____________________________________________________________________________
  /* Write the results to an output file.
   */
// _____________________________________________________________________________
  std::cout << "print results\n";
  std::string filename = "electrode_charge.dat";
  std::fstream file;
  file.open(filename, std::ios::out);
  file << "# [x] [y] [z] [surface charge density]\n";
  file << "# total charge: ";
  file << "normal: " << total_charge_normal << " ";
//  file << "Poisson: " << total_charge_poisson << " ";
  file << "- expected: +-" << diff_expected;
  file << "\n";
  file << "# grid_counts=";
  for (size_t gc : grid_counts) {
    file << gc << ",";
  }
  file << std::endl;
  for (size_t i = 0; i  < charge_normal.size(); i++) {
    for (double pos : charge_normal.at(i).first) {
      file << pos << " ";
    }
    file << charge_normal.at(i).second << std::endl; //" ";
  }
  file.close();
  // If something is wrong: check potential gradients
//  double grad = -2. * potential / system_lengths.at(2);
//  std::cout << "expected potential gradient: " << grad << std::endl;
//  std::cout << "expected potential at bin size " << bin_sizes.at(2) << ": ";
//  std::cout << potential + grad * bin_sizes.at(2) << std::endl;
//  file.open("potential.dat", std::ios::out);
//  file << "# [x] [y] [z] [electrostatic potential]\n";
//  system.print_data(total_ES, file);
//  file.close();
  // all done
  std::cout << "see results in " << filename << std::endl;
  std::cout << "relative deviation from expected results: ";
  std::cout << (total_charge_normal - diff_expected) / diff_expected << std::endl;
  return 0;
}
