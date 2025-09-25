// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file tests/cartesian_poisson_solver/src/main.cpp
  * \brief Calculate electric potential between two structured electrodes.
  * 
  * input parameters via file/command line (cf. ParameterHandler),
  * output to file "total_ES.dat": potential profiles [kT/e]
  */
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundaries.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_sine.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/boundary_surface_cylinder.hpp"
#include "../../../src/cartesian_poisson_solver_any.hpp"
#include "../../../src/cartesian_poisson_solver.hpp"
#include "../../../src/system.hpp"
#include <vector>
#include <string> // for parameter labels
#include <cmath> // std::exp
#include <math.h> // M_PI
#include <iostream>
#include <fstream> // write results to file

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " -f <parameter file>\n";
    std::cout << "(or use flag -p <parameter>=<value>, cf. ParameterHandler)\n";
    exit(1);
  }
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
  // Define System
  System<3> system(&params);
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> lengths = system.system_lengths;
  std::vector<double> bin_sizes = system.bin_sizes;
  // potential on the electrodes
  double potential = params.get_double("potential");
  std::cout << "potential = " << potential << " V "; // in Volt
  potential *= ELECTRON_CHARGE / (BOLTZMANN * system.temperature) * 1e4;
  std::cout << "= " << potential << " kT/e\n"; // in reduced units
  // Species properties:
  std::string type_densities = params.get_string("type_densities");
  std::vector<DataFrame<3,double>> densities(0);
  // Define ion distribution (Poisson rhs)
  DataFrame<3, double> rhs(grid_counts);
  rhs.set_all_elements_to(0.);
  if (type_densities == "none") {} // leave all species empty
  else {
    std::vector<Properties> spec;
    // bulk densities
    double bulk_dens;
    double valency;
    // density profiles
    DataFrame<3,double> spec_density(grid_counts);
    for (size_t s : system.affected_species_es) {
      system.species_properties.at(s).get_property("bulk density", &bulk_dens);
      // Define densities according to input parameter "type_densities"
      if (type_densities == "bulk") { spec_density.set_all_elements_to(bulk_dens); }
      else if (type_densities == "irreg") {
	double x,y,z;
	double x0,y0,z0;
	double dist_squared;
	double width;
        system.get_property("ion_pos_x", &x0);
        system.get_property("ion_pos_y", &y0);
        system.get_property("ion_pos_z", &z0);
        system.get_property("ion_width", &width);
        for (size_t i = 0; i < grid_counts.at(0); i++) {
	  x = static_cast<double>(i) * bin_sizes.at(0);
  	  for (size_t j = 0; j < grid_counts.at(1); j++) {
	    y = static_cast<double>(j) * bin_sizes.at(1);
            for (size_t k = 0; k < grid_counts.at(2); k++) {
	      z = static_cast<double>(k) * bin_sizes.at(2);
  	      // - just some completely arbitrary, unphysical distribution
  	      //spec_density.at(i,j,k) = bulk_dens * (1.e-2 * static_cast<double>(s * k) +
  	      //    (static_cast<double>(s)-1.) * 1.e-4 * static_cast<double>(i*i));
	      // put a blob of ions at position (x0,y0,z0)
              dist_squared = pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2);    
  	      spec_density.at(i,j,k) = bulk_dens * std::exp(-dist_squared/width);
  	    }
  	  }
        }
      }
      else {
        std::cout << "invalid type_densities '" << type_densities << "'\n";
        exit(1);
      }
      system.species_properties.at(s).get_property("valency", &valency);
      rhs += spec_density * valency;
    }
  }
  rhs *= - 4. * M_PI * system.bjerrum;
  // Define Boundaries
  // Set boundary positions
  std::vector<size_t> xy_counts{grid_counts.at(0), grid_counts.at(1)};
  std::vector<DataFrame<2,double>>
	  boundary_positions(2, DataFrame<2,double>(xy_counts));
  // for CartesianPoissonSolverAny:
  Boundaries<3> boundaries;
  std::vector<double> boundary_values = {potential, -potential};
  // for CartesianPoissonSolver:
  std::vector<std::vector<double>> boundary_values2(3);
  boundary_values2.at(2) = boundary_values;
  //BoundarySurfaceSine boundaries_left(system, 0);
  //BoundarySurfaceSine boundaries_right(system, 1);
  BoundarySurfaceCylinder boundaries_left(system, 0);
  BoundarySurfacePlanar boundaries_right(system, 1);
  boundaries_left.set_boundary_value(potential);
  boundaries_right.set_boundary_value(-potential);
  boundaries.add_surface(&boundaries_left);
  boundaries.add_surface(&boundaries_right);
//  } else if (type_ES == "planar") {
//    BoundarySurfacePlanar boundaries_left(properties, 0);
//    BoundarySurfacePlanar boundaries_right(properties, 1);
//    boundaries.add_surface(boundaries_left);
//    boundaries.add_surface(boundaries_right);
//  } else if (type_ES == "zigzag") {
//    // waveform params
//    properties.add_property<double>("amplitude", params.get_double("amplitude", lengths.at(2) / 4));
//    properties.add_property<size_t>("maxima_count", params.get_int("maxima_count", 1));
//    boundaries::zigzag(syst, properties, &boundary_positions);
//  } else {
//    std::cout << "invalid type '" << type_ES << "'\n";
//    exit(1);
//  }
  // Solve Poisson equation with CartesianPoissonSolverAny and - if boundaries
  // are planar - for comparison also with the simpler CartesianPoissonSolver.
  DataFrame<3, double> solution1(grid_counts);
  DataFrame<3, double> solution2(grid_counts);
// -   std::cout << "main: boundary val " << boundaries.surfaces.front()->get_boundary_value() << std::endl;
// -   std::cout << "main: boundary val " << boundaries.surfaces.front()->boundary_value << std::endl;
// -   std::cout << "main: electrost. " << boundaries.surfaces.front()->electrostatics << std::endl;
// -   std::cout << "main: Boundaries address " << &boundaries << std::endl;
// -   std::cout << "main: Surface 1 address " << boundaries.surfaces.front() << std::endl;
// -   std::cout << "main: boundary val address " << &(boundaries.surfaces.front()->boundary_value) << std::endl;
// -   std::cout << "main: electrost. address " << &(boundaries.surfaces.front()->electrostatics) << std::endl;
  std::cout << "CartesianPoissonSolverAny:\n";
  CartesianPoissonSolverAny poisson1(system, &boundaries);
  poisson1.solve(rhs, solution1);
//  if (type_ES == "planar") {
//    std::cout << "CartesianPoissonSolver:\n";
//    CartesianPoissonSolver poisson2(grid_counts, bin_sizes, system.periodic_boundaries);
//    poisson2.solve(rhs, boundary_values2, solution2);
//  }
  // Write results to files
  std::cout << "calculations done, write to file\n";
  std::fstream ES_dat;
  std::string filename{"total_ES.dat"};
  ES_dat.open(filename, std::ios::out);
  ES_dat << "# [x] [y] [z] [electrostatic potential]" << std::endl;
  // write grid dimensions to header
  ES_dat << "# grid_counts=";
  for (size_t g : grid_counts) {
    ES_dat << g << ",";
  }
  ES_dat << std::endl;
  // write data
  double x, y, z;
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    x = bin_sizes.at(0) * static_cast<double>(i);
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      y = bin_sizes.at(1) * static_cast<double>(j);
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
        z = bin_sizes.at(2) * static_cast<double>(k);
//	size_t index = k + j * grid_counts.at(2) + i * grid_counts.at(2) * grid_counts.at(1);
        ES_dat << x << " " << y << " " << z << " ";
	// write both solutions, second is zero if only first Solver was used
        ES_dat << solution1.at(i, j, k) << " " << solution2.at(i, j, k) << std::endl;
      }
    }
  }
  ES_dat.close();
  // all done
  std::cout << "see results in " << filename << std::endl;
}
