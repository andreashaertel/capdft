// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/boundaries/src/main.cpp
 * \brief This file contains code that can be used to visualize and test the
 * boundary surface objects
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundaries.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_sine.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/boundary_surface_cylinder.hpp"
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
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " -f <parameter file>\n";
    std::cout << "(or use flag -p <parameter>=<value>, cf. ParameterHandler)\n";
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
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
  System<3> system(&params);
  // Get some properties from System
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> lengths = system.system_lengths;
  std::vector<double> bin_sizes = system.bin_sizes;
// _____________________________________________________________________________
  /* Construction of the BoundarySurface object
   * All necessary parameters are specified in the system object which is given
   * to the constructor. In order to test the various member functions, we
   * calculate the external hard sphere potential (exp_external_potential_hs),
   * the minimal and directed distances to the boundary (minimal_distances and
   * directed_distance) at different positions in the system and write the
   * results into some data files that can then be visualized e.g. via gnuplot.
   * For the specific case of the BoundarySurfaceCylinder, we also check the
   * cylinder coordinates that this class is using.
   */
// _____________________________________________________________________________
  BoundarySurfaceCylinder boundary(system, 0);
  // BoundarySurfacePlanar boundary(system, 0);
  // Write data to files
  std::fstream file_out;
  DataFrame<3, double> data_frame(grid_counts);
  std::vector<DataFrame<3, double>> data = {data_frame, data_frame};
  // discretize volume
  DataFrame<3, bool> data_frame_bool(grid_counts);
  boundary.discretize_volume(&data_frame_bool);
  file_out.open("discretize_volume.dat", std::ios::out);
  file_out << "# [x] [y] [z] [is within boundaries?]" << std::endl;
  system.print_data(data_frame_bool, file_out);
  file_out.close();
  // exp_external_potential_hs
  boundary.exp_external_potential_hs(system, &data);
  file_out.open("HS_potential.dat", std::ios::out);
  file_out << "# [x] [y] [z] [exp. external hard sphere potential]" << std::endl;
  system.print_data(data, file_out);
  file_out.close();
  double resolution = bin_sizes.at(0);
//  std::vector<std::vector<double>> surface_points;
//  boundary.discretize_surface(&surface_points, resolution);
  boundary.minimal_distances(&data_frame, resolution);
  file_out.open("minimal_distances.dat", std::ios::out);
  file_out << "# [x] [y] [z] [minimal distance to boundary]" << std::endl;
  system.print_data(data_frame, file_out);
  file_out.close();
  // distance_directed
  data.clear();
  std::vector<double> position(3);
  for (size_t dir = 0; dir < 3; dir++) {
    for (bool forward : {true, false}) {
      for (size_t i = 0; i < grid_counts.at(0); i++) {
	position.at(0) = bin_sizes.at(0) * static_cast<double>(i);
        for (size_t j = 0; j < grid_counts.at(1); j++) {
	  position.at(1) = bin_sizes.at(1) * static_cast<double>(j);
          for (size_t k = 0; k < grid_counts.at(2); k++) {
	    position.at(2) = bin_sizes.at(2) * static_cast<double>(k);
            data_frame.at(i,j,k) = boundary.distance_directed(position, dir, forward);
	  }
	}
      }
      data.push_back(data_frame);
    }
  }
  file_out.open("distance_directed.dat", std::ios::out);
  file_out << "# [x] [y] [z] [distances to boundary in different directions]" << std::endl;
  system.print_data(data, file_out);
  file_out.close();
  // cylinder coordinates
  data.resize(3);
  std::vector<double> cylinder_coordinates(3);
  for (size_t i = 0; i < grid_counts.at(0); i++) {
    position.at(0) = bin_sizes.at(0) * static_cast<double>(i);
    for (size_t j = 0; j < grid_counts.at(1); j++) {
      position.at(1) = bin_sizes.at(1) * static_cast<double>(j);
      for (size_t k = 0; k < grid_counts.at(2); k++) {
        position.at(2) = bin_sizes.at(2) * static_cast<double>(k);
	boundary.cylinder_coordinates(position, &cylinder_coordinates);
	for (size_t dir = 0; dir < 3; dir++) {
          data.at(dir).at(i,j,k) = cylinder_coordinates.at(dir);
	}
      }
    }
  }
  file_out.open("cylinder_coordinates.dat", std::ios::out);
  file_out << "# [x] [y] [z] [r] [phi] [z']" << std::endl;
  system.print_data(data, file_out);
  file_out.close();
  // all done
  return 0.;
}
