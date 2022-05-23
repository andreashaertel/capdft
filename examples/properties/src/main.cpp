// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/properties/src/main.cpp
 *  \brief Main file of the example of the Properties class.
 *  
 *  This main file contains examples to show how the Properties class
 *  works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <iostream>
#include <vector>
#include "../../../src/properties.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
  // Create object of Properties class
  Properties properties;
  std::vector<Properties> species_properties;
  // First species
  properties.add_property<int>("number", 1);
  properties.add_property<double>("size", 2.71);
  properties.add_property("is bird", true);  // without type (bool)
  properties.add_property<std::string>("color", "black");
  species_properties.push_back(properties);
  properties.clear();
  // Second species (does not have color).
  properties.add_property<int>("number", 2);
  properties.add_property<double>("size", 3.14);
  properties.add_property<bool>("is bird", false);
  species_properties.push_back(properties);
  properties.clear();
  // Third species (does not have color and no size).
  properties.add_property<int>("number", 3);
  properties.add_property<bool>("is bird", true);
  species_properties.push_back(properties);
  properties.~Properties();
  // Check for number
  int species_number{0};
  for (auto it = species_properties.begin();
      it != species_properties.end(); ++it) {
    it->get_property<int>("number", &species_number);
    std::cout << "This species has the number " << species_number << ".";
    std::cout << std::endl;
  }
  // Check for size
  std::cout << "=== Property: size ===" << std::endl;
  double species_size;
  for (auto it = species_properties.begin();
      it != species_properties.end(); ++it) {
    it->get_property<int>("number", &species_number);
    std::cout << "Species " << species_number;
    if (it->get_property<double>("size", &species_size)) {
      std::cout << " has size " << species_size;
    } else {
      std::cout << " has no size";
    }
    std::cout << "." <<std::endl;
  }
  // Check whether it is a bird
  std::cout << "=== Property: is bird ===" << std::endl;
  bool species_bird;
  for (auto it = species_properties.begin();
      it != species_properties.end(); ++it) {
    it->get_property("number", &species_number);  // works without type (int)
    std::cout << "Species " << species_number;
    if (it->get_property<bool>("is bird", &species_bird)) {
      if (species_bird) {
        std::cout << " is a bird";
      } else {
        std::cout << " is not a bird";
      }
    } else {
      std::cout << " has not clarified whether it identifies as bird";
    }
    std::cout << "." << std::endl;
  }
  // Check for color
  std::cout << "=== Property: color ===" << std::endl;
  std::string species_color;
  for (auto it = species_properties.begin();
      it != species_properties.end(); ++it) {
    it->get_property<int>("number", &species_number);
    std::cout << "Species " << species_number;
    if (it->get_property<std::string>("color", &species_color)) {
      std::cerr << " is " << species_color;
    } else {
      std::cerr << " has no color";
    }
    std::cout << "." << std::endl;
  }
  return 0;
}
