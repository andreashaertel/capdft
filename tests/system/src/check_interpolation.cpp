// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/system/src/check_interpolation.cpp
 * \brief Test the interpolation methods of the System class
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/system.hpp"
#include <vector>
#include <string> // for parameter labels
#include <cmath> // std::exp
#include <math.h> // M_PI
#include <iostream>
#include <fstream> // write results to file
// _____________________________________________________________________________
// Test function to be interpolated (we only use x and z as parameters, but one
// could also include y)
double function(double scale, double x, double z) {
  return scale * (z + sin(x * 6));
}
// _____________________________________________________________________________
// Main function
int main(int argc, char** argv) {
// _____________________________________________________________________________
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
// _____________________________________________________________________________
  /* We define two grids - the second one has a finer resolution. We want to
   * calculate the values at the new grid points via interpolation between the
   * original ones.
   */
// _____________________________________________________________________________
  // Original grid:
  System<3> system(&params);
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> bin_sizes = system.bin_sizes;
  // New grid:
  std::vector<size_t> new_grid_counts(3);
  std::vector<double> new_bin_sizes(3);
  for (size_t i = 0; i < 3; i++) {
    new_grid_counts.at(i) = 2 * grid_counts.at(i);
    new_bin_sizes.at(i) = 0.5 * bin_sizes.at(i);
  }
// _____________________________________________________________________________
  /* Define the values on the original grid based on the previously defined test
   * function. For instance, these values could represent an electrostatic
   * potential profile.
   */
// _____________________________________________________________________________
  double potential = params.get_double("potential");
  std::cout << "define data field\n";
  DataFrame<3,double> total_ES(grid_counts);
  double x, z, value;
  for (size_t k = 0; k < grid_counts.at(2); ++k) {
    z = static_cast<double>(k) * bin_sizes.at(2);
    for (size_t i = 0; i < grid_counts.at(0); ++i) {
      x = static_cast<double>(i) * bin_sizes.at(0);
      value = function(potential, x, z); 
      for (size_t j = 0; j < grid_counts.at(1); ++j) {
	total_ES.at(i,j,k) = value;
      }
    }
  }
// _____________________________________________________________________________
  /* Calculate the values on the new grid via interpolation.
   */
// _____________________________________________________________________________
  std::cout << "interpolate\n";
  std::vector<double> position(3);
  std::fstream file;
  file.open("interpolation.dat", std::ios::out);
  file << "# [x] [y] [z] [calculated] [expected]\n";
  file << "# grid_counts=";
  for (size_t gc : new_grid_counts) {
    file << gc - 3 << ",";
  }
  file << std::endl;
  // optionally, also test other system methods here:
//  std::fstream file2;
//  file2.open("laplace.dat", std::ios::out);
//  file2 << "# [x] [y] [z] [calculated]\n";
//  file2 << "# grid_counts=";
//  for (size_t gc: new_grid_counts) {
//    file2 << gc - 3 << ",";
//  }
//  file2 << std::endl;
  double expected;
  double deviation;
  double max_deviation = 0.;
  for (size_t i = 1; i < new_grid_counts.at(0) - 2; ++i) {
    position.at(0) = static_cast<double>(i) * new_bin_sizes.at(0);
    for (size_t j = 1; j < new_grid_counts.at(1) - 2; ++j) {
      position.at(1) = static_cast<double>(j) * new_bin_sizes.at(1);
      for (size_t k = 1; k < new_grid_counts.at(2) - 2; ++k) {
        position.at(2) = static_cast<double>(k) * new_bin_sizes.at(2);
	// print data
        for (double pos : position) {
	  file << pos << " ";
//	  file2 << pos << " ";
	}
	value = system.interpolate(total_ES, position);
        expected = function(potential, position.at(0), position.at(2));
	deviation = std::fabs(value - expected);
	if (deviation > max_deviation) {
	  max_deviation = deviation;
	}
	file << value << " " << expected << std::endl;
//	value = system.laplace(total_ES, position, value);
//	file2 << value << std::endl;
      }
    }
  }
  file.close();
//  file2.close();
// _____________________________________________________________________________
  // all done
  std::cout << "maximum error of interpolated values: ";
  std::cout << max_deviation << std::endl;
  return 0;
}
