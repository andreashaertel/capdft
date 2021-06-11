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
  // Create object of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Define system properties
  system_properties.add_property<double>("length", 10.);
  system_properties.add_property<double>("bjerrum length", 1.);
  system_properties.add_property<int>("grid_count", 1001);
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
  return 0;
}
