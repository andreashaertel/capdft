// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
 *  \brief Main file of the example of the spherical functionals.
 *  
 *  This main file contains examples to show how the FunctionalESMFSpherical
 *  and FunctionalFMTSpherical classes work.
 *
 */
// _____________________________________________________________________________
// Includes
#include <iostream>
#include <vector>
#include "../../../src/data_field.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/system.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_spherical.hpp"
#include "../../../src/functional_es_mf_spherical.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
  // Set the desired system properties
  double system_length = 10.;
  double bjerrum_length = 1.;
  size_t grid_count = 1001;
  // Create objects of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length", system_length);
  system_properties.add_property<double>("bjerrum length", bjerrum_length);
  system_properties.add_property<size_t>("grid count", grid_count);
  // First species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .1);
  properties.add_property<double>("valency", -1.);
  species_properties.push_back(properties);
  properties.clear();
  // Second species
  properties.add_property<double>("bulk density", .2);
  properties.add_property<double>("valency", +1.);
  species_properties.push_back(properties);
  properties.clear();
  // Third species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .3);
  species_properties.push_back(properties);
  properties.clear();
  // All the supplied data is now brought together in the System class.
  // The system class creates a density profile (DataField) from the given
  // information and does not allow modifictation of the properties.
  System my_system(system_properties, species_properties);
  // Create FMT Functional object
  FunctionalFMTSpherical my_fmt_functional(&my_system);
  double energy;
  energy = my_fmt_functional.calc_energy();
  std::cout << "Excess free energy of FMT functional: ";
  std::cout << energy << std::endl;
  // Create mean-field electrostatic Functional object
  //FunctionalESMFSpherical my_es_functional(
  //  system_properties, species_properties, density_profile);
  return 0;
}
