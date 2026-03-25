// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/system/src/interpolate_data.cpp
 * \brief Interpolate data from a given data file to a new grid
 *
 * Required input parameters are listed in ../params_interpolation.txt.
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
// Main function
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " -f <parameter file>\n";
    std::cout << "	(or use -p <parameter_name>=<value>)\n";
    exit(1);
  }
// _____________________________________________________________________________
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
// _____________________________________________________________________________
  /* We define two grids. We want to
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
    system.get_property("new_grid_count", &new_grid_counts.at(i), i);
    new_bin_sizes.at(i) = system.system_lengths.at(i) /
	    		  static_cast<double>(new_grid_counts.at(i));
  }
// _____________________________________________________________________________
  /* Extract original data from data file.
   */
// _____________________________________________________________________________
  std::string input_file = params.get_string("input");
  std::vector<size_t> col_nums(0);
  { // workaround because parameter handler has no vector<size_t>:
    std::vector<int> temp = params.get_vector_int("col_nums");
    for (int num : temp) {
      col_nums.push_back(static_cast<size_t>(num));
    }
  }
  std::vector<DataFrame<3,double>> data(col_nums.size(),
	  				DataFrame<3,double>(grid_counts));
  system.load_data(input_file, col_nums, &data);
// _____________________________________________________________________________
  /* Calculate the values on the new grid via interpolation.
   */
// _____________________________________________________________________________
  std::cout << "interpolate\n";
  std::string output_file = params.get_string("output");
  std::fstream file;
  file.open(output_file, std::ios::out);
  file << "# data interpolated from " << input_file << std::endl;
  file << "# [x] [y] [z] [interpolated data]\n";
  file << "# grid_counts=";
  for (size_t gc : new_grid_counts) {
    file << gc << ",";
  }
  file << std::endl;
  std::vector<double> position(3);
  double value;
  for (size_t i = 0; i < new_grid_counts.at(0); ++i) {
    position.at(0) = static_cast<double>(i) * new_bin_sizes.at(0);
    for (size_t j = 0; j < new_grid_counts.at(1); ++j) {
      position.at(1) = static_cast<double>(j) * new_bin_sizes.at(1);
      for (size_t k = 0; k < new_grid_counts.at(2); ++k) {
        position.at(2) = static_cast<double>(k) * new_bin_sizes.at(2);
	// print data
        for (double pos : position) {
	  file << pos << " ";
	}
	for (size_t col = 0; col < col_nums.size(); col++) {
	  value = system.interpolate(data.at(col), position);
	  file << value << " ";
	}
	file << std::endl;
      }
    }
  }
  file.close();
  return 0;
}
