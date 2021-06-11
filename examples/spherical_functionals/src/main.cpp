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
#include "../../../src/properties.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
  // Create object of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // First species
  properties.add_property<int>("number", 1);
  properties.add_property<double>("size", 2.71);
  properties.add_property("is bird", true);
  properties.add_property<std::string>("color", "black");
  species_properties.push_back(properties);
  properties.clear();
  return 0;
}
