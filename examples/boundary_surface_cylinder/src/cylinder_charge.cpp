// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/boundary_surface_cylinder/src/cylinder_charge.cpp
 *  \brief Calculate the charge distribution on the cylinder surface
 *  from the electrostatic potential profile.
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
  System<3> system(&params);
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
// _____________________________________________________________________________
  /* Extract the electrostatic potential profile from a data file
   * The filename is specified in the command line parameters.
   */
// _____________________________________________________________________________
  DataFrame<3,double> total_ES(grid_counts);
  std::string potential_file = params.get_string("potential_profile", "total_ES.dat");
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
  // initialize
  std::vector<double> position(3);
  std::vector<double> normal(3);
  double radius = params.get_double("0_radius");
  double radius_squared = radius * radius;
  double half_height = params.get_double("0_height") * 0.5;
  double charge;
  double x_mid = system.system_lengths.at(0) * 0.5;
  double y_mid = system.system_lengths.at(1) * 0.5;
  // mantle
  normal.at(2) = 0.;
  double phi;
  double dphi = system.bin_sizes.at(0) / radius;
  double z;
  std::fstream out;
  out.open("surface_charge_mantle.dat", std::ios::out);
  out << "# [phi] [z] [surface charge density]\n";
  z = system.system_lengths.at(2) * 0.5 - half_height;
  while (z <= system.system_lengths.at(2) * 0.5 + half_height) {
    position.at(2) = z;
    phi = 0.;
    while (phi < 2. * M_PI) {
      normal.at(0) = cos(phi);
      normal.at(1) = sin(phi);
      position.at(0) = x_mid + normal.at(0) * radius;
      position.at(1) = y_mid + normal.at(1) * radius;
      charge = charge_density(cylinder, system, total_ES, position, normal);
      out << phi << " " << z << " " << charge << std::endl;
      phi += dphi;
    }
    z += system.bin_sizes.at(2);
  }
  out.close();
  // upper cylinder cap
  normal.at(0) = 0.;
  normal.at(1) = 0.;
  normal.at(2) = 1.;
  double x, y;
  position.at(2) = system.system_lengths.at(2) / 2. + half_height;
  out.open("surface_charge_top.dat", std::ios::out);
  out << "# [x] [y] [surface charge density]\n";
  x = 0.;
  while (x < system.system_lengths.at(0)) {
    position.at(0) = x;
    y = 0.;
    while (y < system.system_lengths.at(1)) {
      position.at(1) = y;
      if (pow(x - x_mid, 2) + pow(y - y_mid, 2) > radius_squared) {
	charge = 0.;
      } else {
        charge = charge_density(cylinder, system, total_ES, position, normal);
      }
      out << x << " " << y << " " << charge << std::endl;
      y += system.bin_sizes.at(1);
    }
    x += system.bin_sizes.at(0);
  }
  out.close();
  // lower cylinder cap
  normal.at(0) = 0.;
  normal.at(1) = 0.;
  normal.at(2) = -1.;
  position.at(2) = system.system_lengths.at(2) / 2. - half_height;
  out.open("surface_charge_bottom.dat", std::ios::out);
  out << "# [x] [y] [surface charge density]\n";
  double x_diff, y_diff;
  x = 0.;
  while (x < system.system_lengths.at(0)) {
    position.at(0) = x;
    y = 0.;
    while (y < system.system_lengths.at(1)) {
      position.at(1) = y;
      x_diff = x - x_mid;
      y_diff = y - y_mid;
      if (x_diff * x_diff + y_diff * y_diff > radius_squared) {
	charge = 0.;
      } else {
        charge = charge_density(cylinder, system, total_ES, position, normal);
      }
      out << x << " " << y << " " << charge << std::endl;
      y += system.bin_sizes.at(1);
    }
    x += system.bin_sizes.at(0);
  }
  out.close();
  // Could also use the function charge_distribution, but then you might have to
  // sort the values by positions before you can plot the results:
  std::vector<std::pair<std::vector<double>, double>> charge_distr;
  double total_charge;
  charge_distribution(cylinder, system, total_ES,
		      &charge_distr, &total_charge, resolution);
  out.open("surface_charge_total.dat", std::ios::out);
  out << "# [x] [y] [z] [surface charge density]\n";
  out << "# total charge: " << total_charge << std::endl;
  for (std::pair<std::vector<double>, double>& surface_point : charge_distr) {
    for (size_t dir = 0; dir < 3; dir++) {
      out << surface_point.first.at(dir) << " ";
    }
    out << surface_point.second << std::endl;
  }
  out.close();
  std::cout << "total charge: " << total_charge << std::endl;
// _____________________________________________________________________________
  /* Alternative method to calculate total surface charge:
   * Use total charge neutrality of the system and integrate over all ion
   * charges.
  */
//____________________________________________________________________________
  // get species properties
  size_t species_count = system.species_properties.size();
  std::vector<DataFrame<3,double>> densities(species_count,
		  DataFrame<3,double>(grid_counts));
  std::vector<size_t> data_columns(species_count);
  std::vector<double> valencies(species_count);
  for (size_t s = 0; s < species_count; s++) {
    data_columns.at(s) = s + 3;
    system.species_properties.at(s).get_property("valency", &valencies.at(s));
  }
  // get density profiles
  std::string densities_file = params.get_string("density_profiles", "3d_profiles.dat");
  system.load_data(densities_file, data_columns, &densities);
  // integrate charge densities
  double ion_count;
  double total_ion_charge = 0.;
  for (size_t s = 0; s < species_count; s++) {
    ion_count = integration_3d_closed(densities.at(s), system.bin_sizes);
    total_ion_charge += ion_count * valencies.at(s);
  }
  std::cout << "total ion charge: " << total_ion_charge << std::endl;
// _____________________________________________________________________________
  /* Check consistency of the two methods' results
   */
// _____________________________________________________________________________
  std::cout << "total charge in system: ";
  std::cout << total_ion_charge + total_charge;
  std::cout << std::endl;
  return 0;
}
