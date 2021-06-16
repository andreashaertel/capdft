// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/properties/src/main.cpp
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
#include "../../../src/properties.hpp"
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
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .1);
  species_properties.push_back(properties);
  properties.clear();
  // Third species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .1);
  properties.add_property<double>("valency", +1.);
  species_properties.push_back(properties);
  properties.clear();
  // Create a density array
  double** density_profiles = new double*[species_properties.size()];
  for (size_t spec_i = 0; spec_i < species_properties.size(); ++spec_i) {
    density_profiles[spec_i] = new double[grid_count];
  }
  // Create FMT Functional object
  FunctionalFMTSpherical functional2(
    system_properties, species_properties, density_profiles);
  // Create mean-field electrostatic Functional object
  FunctionalESMFSpherical functional3(
    system_properties, species_properties, density_profiles);
  return 0;
}
