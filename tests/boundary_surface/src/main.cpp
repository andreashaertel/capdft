// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/boundary_surface/src/main.cpp
 * \brief This file contains code that can be used to visualize and test the
 * boundary surface objects
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_sine.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/boundary_surface_cylinder.hpp"
#include "../../../src/boundary_surface_sphere.hpp"
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
    std::cout << "	(or use -p <parameter_name>=<value>)\n";
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
  std::cout << "Initialize system...\n";
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
   * calculate the external potential acting on hard sphere particles 
   * (exp_external_potential_hs) and
   * the minimal and directed distances to the boundary (minimal_distances and
   * directed_distance) at different positions in the system. We write the
   * results into some data files that can then be visualized e.g. via gnuplot.
   * For the specific case of the BoundarySurfaceCylinder, we also
   * check the cylinder coordinates that this class is using.
   */
// _____________________________________________________________________________
  // BoundarySurfaceCylinder boundary(system, 0);
  BoundarySurfaceSphere boundary(system, 0);
  // BoundarySurfaceSine boundary(system, 0);
  double resolution;
  // Write data to files
  std::fstream file_out;
  DataFrame<3, double> data_frame(grid_counts);
  std::vector<DataFrame<3, double>> data = {data_frame, data_frame};
  // discretize volume
  std::cout << "> discretize_volume...\n";
  DataFrame<3, bool> data_frame_bool(grid_counts);
  boundary.discretize_volume(&data_frame_bool);
  file_out.open("discretize_volume.dat", std::ios::out);
  file_out << "# [x] [y] [z] [is within boundaries?]" << std::endl;
  system.print_data(data_frame_bool, file_out);
  file_out.close();
  // exp_external_potential_hs
  std::cout << "> HS_potential...\n";
  boundary.exp_external_potential_hs(system, &data);
  file_out.open("HS_potential.dat", std::ios::out);
  file_out << "# [x] [y] [z] [exp. external hard sphere potential]" << std::endl;
  system.print_data(data, file_out);
  file_out.close();
  resolution = bin_sizes.at(0);
  // minimal_distances
  std::cout << "> minimal_distances...\n";
  boundary.minimal_distances(&data_frame, resolution);
  file_out.open("minimal_distances.dat", std::ios::out);
  file_out << "# [x] [y] [z] [minimal distance to boundary]" << std::endl;
  system.print_data(data_frame, file_out);
  file_out.close();
  // distance_directed
  std::cout << "> distance_directed...\n";
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
////  // cylinder coordinates
//  std::cout << "> cylinder_coordinates...\n";
////  data.resize(3);
////  std::vector<double> cylinder_coordinates(3);
////  for (size_t i = 0; i < grid_counts.at(0); i++) {
////    position.at(0) = bin_sizes.at(0) * static_cast<double>(i);
////    for (size_t j = 0; j < grid_counts.at(1); j++) {
////      position.at(1) = bin_sizes.at(1) * static_cast<double>(j);
////      for (size_t k = 0; k < grid_counts.at(2); k++) {
////        position.at(2) = bin_sizes.at(2) * static_cast<double>(k);
////	boundary.cylinder_coordinates(position, &cylinder_coordinates);
////	for (size_t dir = 0; dir < 3; dir++) {
////          data.at(dir).at(i,j,k) = cylinder_coordinates.at(dir);
////	}
////      }
////    }
////  }
////  file_out.open("cylinder_coordinates.dat", std::ios::out);
////  file_out << "# [x] [y] [z] [r] [phi] [z']" << std::endl;
////  system.print_data(data, file_out);
////  file_out.close();
//  // surface discretization - can't be tested when
//  // BoundarySurface::discretize is a private member function!
//  // Make it public before executing the following code.
//  std::cout << "> discretize_surface...\n";
//  std::vector<std::vector<std::vector<double>>> surface;
//  resolution = params.get_double("resolution", bin_sizes.at(0));
//  boundary.discretize_surface(&surface, resolution);
//  file_out.open("discretize_surface.dat", std::ios::out);
//  file_out << "# [x] [y] [z] [surface normal vector]" << std::endl;
//  for (std::vector<std::vector<double>>& point : surface) {
//    for (double pos : point.at(0)) { // positions
//      file_out << pos << " ";
//    }
//    for (double dir : point.at(1)) { // normal vector
//      file_out << dir << " "; 
//    }
//    file_out << std::endl;
//  }
//  file_out.close();
  // all done
  return 0.;
}
