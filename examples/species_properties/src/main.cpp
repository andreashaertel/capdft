/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
/** \file examples/species_properties/src/main.cpp
 *  \brief Main file of the example of the SpeciesProperties class.
 *  
 *  This main file contains examples to show how the SpeciesProperties class
 *  works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <iostream>
#include "../../../src/species_properties.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
  // Create object of SpeciesProperties class
  SpeciesProperties species_properties;
  std::vector<SpeciesProperties> species_properties_vector;
  // First species
  species_properties.add_property("number", 1);
  species_properties.add_property("size", 2.71);
  species_properties.add_property("is bird", true);
  species_properties.add_property("color", std::string("black"));
  species_properties_vector.push_back(species_properties);
  species_properties.clear();
  // Second species (does not have color).
  species_properties.add_property("number", 2);
  species_properties.add_property("size", 3.14);
  species_properties.add_property("is bird", false);
  species_properties_vector.push_back(species_properties);
  species_properties.clear();
  // Third species (does not have color and no size).
  species_properties.add_property("number", 3);
  species_properties.add_property("is bird", true);
  species_properties_vector.push_back(species_properties);
  species_properties.~SpeciesProperties();
  // Check for number
  int species_number;
  for (auto it = species_properties_vector.begin();
      it != species_properties_vector.end(); ++it) {
    it->get_property("number", &species_number);
    std::cout << "This species has the number " << species_number << ".";
    std::cout << std::endl;
  }
  // Check for size
  std::cout << "=== Property: size ===" << std::endl;
  double species_size;
  for (auto it = species_properties_vector.begin();
      it != species_properties_vector.end(); ++it) {
    it->get_property("number", &species_number);
    std::cout << "Species " << species_number;
    if (it->get_property("size", &species_size)) {
      std::cout << " has size " << species_size;
    } else {
      std::cout << " has no size";
    }
    std::cout << "." <<std::endl;
  }
  // Check whether it is a bird
  std::cout << "=== Property: is bird ===" << std::endl;
  bool species_bird;
  for (auto it = species_properties_vector.begin();
      it != species_properties_vector.end(); ++it) {
    it->get_property("number", &species_number);
    std::cout << "Species " << species_number;
    if (it->get_property("is bird", &species_bird)) {
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
  for (auto it = species_properties_vector.begin();
      it != species_properties_vector.end(); ++it) {
    it->get_property("number", &species_number);
    std::cout << "Species " << species_number;
    if (it->get_property("color", &species_color)) {
      std::cerr << " is " << species_color;
    } else {
      std::cerr << " has no color";
    }
    std::cout << "." << std::endl;
  }
  return 0;
}
