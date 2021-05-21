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
  SpeciesProperties my_properties;
  // Add two species dependent properties
  my_properties.add_property("number", std::vector<int>{1, 2});
  my_properties.add_property("size", std::vector<double>{10.2, 3.4});
  my_properties.add_property("is bird", std::vector<bool>{true, false});
  my_properties.add_property("color", std::vector<std::string>{"red", "black"});
  // Retrieve the properties
  std::vector<int> number;
  std::vector<double> size;
  std::vector<bool> is_bird;
  std::vector<std::string> color;
  my_properties.get_property("number", &number);
  my_properties.get_property("size", &size);
  my_properties.get_property("is bird", &is_bird);
  my_properties.get_property("color", &color);
  // Print the properties
  std::cout << "Number: ";
  for (auto it = number.begin(); it != number.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
  std::cout << "Size: ";
  for (auto it = size.begin(); it != size.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
  std::cout << "Is it a bird?: ";
  for (auto it = is_bird.begin(); it != is_bird.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
  std::cout << "Color: ";
  for (auto it = color.begin(); it != color.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
  return 0;
}
