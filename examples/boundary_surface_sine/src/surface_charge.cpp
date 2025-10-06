// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/boundary_surface_sine/src/surface_charge.cpp
 *  \brief Calculate the charge distribution on a system's boundary surfaces
 *  from the electrostatic potential profile.
 */
// _____________________________________________________________________________
// Includes
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_sine.hpp"
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
  BoundarySurfaceSine surface_left(system, 0);
  BoundarySurfaceSine surface_right(system, 1);
  surface_left.set_boundary_value(potential);
  surface_right.set_boundary_value(-potential);
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
  /* From the potential profile and the boundary objects, we can calculate the
   * charge distribution on the boundary surfaces.
   */
// _____________________________________________________________________________
  std::vector<std::pair<std::vector<double>, double>> charge_left(0);
  std::vector<std::pair<std::vector<double>, double>> charge_right(0);
  double total_charge_left, area_left;
  double total_charge_right, area_right;
  std::cout << "calculate charge distribution\n";
  charge_distribution(surface_left, system, total_ES,
		  &charge_left, &total_charge_left, resolution, &area_left);
  charge_distribution(surface_right, system, total_ES,
		  &charge_right, &total_charge_right, resolution, &area_right);
// _____________________________________________________________________________
  /* Write results to file
   * Since the surface points yielded by BoundarySurfaceSine::discretize_surface
   * which are used in the charge_distribution function are well-ordered on a
   * square grid, we can write the results to a data file in exactly the order
   * in which they are returned. For other BoundarySurface subclasses, this
   * might not be the case (e.g. BoundarySurfaceCylinder) and we would have to
   * be careful with the ordering if we want to use the resulting data file for
   * further analysis or plot it.
   */
// _____________________________________________________________________________
  std::cout << "print results\n";
  // Integrated over the whole surface
  std::fstream file_total;
  file_total.open("electrode_charge_total.dat", std::ios::out);
  file_total << "total charge [e]: ";
  file_total << total_charge_left << " (left), "
	     << total_charge_right << " (right)\n";
  file_total << "total area [nm²]: ";
  file_total << area_left << " (left), "
	     << area_right << " (right)\n";
  file_total << "mean charge per surface area [e/nm²]: ";
  file_total << total_charge_left / area_left << " (left), "
	     << total_charge_right / area_right << " (right)\n";
  // Spatial distribution on the surface
  std::fstream file;
  // left surface
  file.open("electrode_charge_left.dat", std::ios::out);
  file << "# Surface charge calculated from '" << potential_file << "'\n";
  file << "# [x] [y] [z] [surface charge density [e]]\n";
  file << "# total charge: " << total_charge_left << std::endl;
  file << "# grid_counts=" << grid_counts.at(0) << "," << grid_counts.at(1) << std::endl;
  for (size_t i = 0; i  < charge_left.size(); i++) {
    for (double pos : charge_left.at(i).first) {
      file << pos << " ";
    }
    file << charge_left.at(i).second << std::endl;
  }
  file.close();
  // right surface
  file.open("electrode_charge_right.dat", std::ios::out);
  file << "# Surface charge calculated from '" << potential_file << "'\n";
  file << "# [x] [y] [z] [surface charge density [e]]\n";
  file << "# total charge: " << total_charge_right << std::endl;
  file << "# grid_counts=" << grid_counts.at(0) << "," << grid_counts.at(1) << std::endl;
  for (size_t i = 0; i  < charge_right.size(); i++) {
    for (double pos : charge_right.at(i).first) {
      file << pos << " ";
    }
    file << charge_right.at(i).second << std::endl;
  }
  file.close();
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
  file_total << "total ion charge: " << total_ion_charge << std::endl;
// _____________________________________________________________________________
  /* Check consistency of the two methods' results
   */
// _____________________________________________________________________________
  file_total << "total charge in system: ";
  file_total << total_ion_charge + total_charge_left + total_charge_right;
  file_total << std::endl;
  file_total.close();
  return 0;
}
