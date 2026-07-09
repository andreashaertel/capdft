// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file system.cpp
 *  \brief Source file for the System class.
 *
 *  The file contains the definitions of the System class.
 */
#include "../../parameter_handler/src/parameter_handler.hpp"
#include "properties.hpp"
#include "system.hpp"
#include "constants.hpp"
#include "data_frame.hpp"
#include "discretization_handler.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

// Explicit instantiation for different geometries
template class System<Cartesian>;
template class System<Planar>;
//template class System<Spherical>;

// _____________________________________________________________________________
template <Geometry geometry>
System<geometry>::System(Properties& properties) : 
	DiscretizationHandler<geometry>(properties) {
  extract_electrical_properties();
  extract_species_properties();
}
// _____________________________________________________________________________
template <Geometry geometry>
System<geometry>::System(ParameterHandler* parameters) : 
	DiscretizationHandler<geometry>(parameters) {
  extract_electrical_properties();
  extract_species_properties();
}
// _____________________________________________________________________________
template <Geometry geometry>
System<geometry>::System(Properties& properties, 
		std::vector<size_t>& species_indices) : 
	DiscretizationHandler<geometry>(properties) {
  extract_electrical_properties();
  extract_species_properties(species_indices);
}
// _____________________________________________________________________________
template <Geometry geometry>
System<geometry>::System(ParameterHandler* parameters, 
		std::vector<size_t>& species_indices) : 
	DiscretizationHandler<geometry>(parameters) {
  extract_electrical_properties();
  extract_species_properties(species_indices);
}
// _____________________________________________________________________________
// Copied over from functional_es_mf::extract_electrical_properties() because
// this is also needed outside of the functional's method.
template <Geometry geometry>
void System<geometry>::extract_electrical_properties() {
  // This variable saves a number representing, which of the el. variables are
  // available.
  int cases{0};
  // Get the variables
  try {
    if (this->get_property("temperature", &temperature)) {
      cases += 1;
    }
  } catch (...) {}
  try {
    if (this->get_property("bjerrum length", &bjerrum)) {
      cases += 2;
    }
  } catch (...) {}
  try {
    if (this->get_property("dielectric constant", &dielectric)) {
      cases += 4;
    }
  } catch (...) {}
  // Auxiliary variable e^2 / (4 * \epsilon_0 * \pi * k_B)
  double e2_by_4pi_boltz = 1e-3 * ELECTRON_CHARGE * ELECTRON_CHARGE /
      (4. * M_PI * VACUUM_PERMITIVITY * BOLTZMANN);
  // Calculate missing variables
  switch (cases) {
    case 3:  // dielectric constant missing
      dielectric = e2_by_4pi_boltz / (bjerrum * 1e-9 * temperature);
      break;
    case 5:  // bjerrum length missing
      bjerrum = 1e9 * e2_by_4pi_boltz / (dielectric * temperature);
      break;
    case 6:  // temperature missing
      temperature = e2_by_4pi_boltz / (dielectric * bjerrum * 1e-9);;
      break;
    case 7:  // nothing missing
      std::cerr << "extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
template <Geometry geometry>
void System<geometry>::extract_species_properties(
		std::vector<size_t>& species_indices) {
  // initialize
  size_t number = species_indices.size();
  species_properties.resize(number);
  affected_species_fmt.clear();
  affected_species_es.clear();
  double bulk_density;
  double diameter;
  double valency;
  // extract properties
  for (size_t i = 0; i < number; i++) {
    species_properties.at(i).clear();
    this->get_property("bulk density", &bulk_density, species_indices.at(i));
    species_properties.at(i).add_property("bulk density", bulk_density);
    try {
      this->get_property("diameter", &diameter, species_indices.at(i));
      species_properties.at(i).add_property("diameter", diameter);
      affected_species_fmt.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    try {
      this->get_property("valency", &valency, species_indices.at(i));
      species_properties.at(i).add_property("valency", valency);
      affected_species_es.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
  }
}
// _____________________________________________________________________________
template <Geometry geometry>
void System<geometry>::extract_species_properties() {
  // initialize
  species_properties.clear();
  affected_species_fmt.clear();
  affected_species_es.clear();
  double bulk_density;
  double diameter;
  double valency;
  // extract properties (indices i = 0,1,... as long as bulk density is defined)
  size_t i = 0;
  Properties properties;
  while (true) {
    properties.clear();
    try {
      this->get_property("bulk density", &bulk_density, i);
      properties.add_property("bulk density", bulk_density);
    } catch (const Properties::MissingPropertyException*) {
      break;
    }
    try {
      this->get_property("diameter", &diameter, i);
      properties.add_property("diameter", diameter);
      affected_species_fmt.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    try {
      this->get_property("valency", &valency, i);
      properties.add_property("valency", valency);
      affected_species_es.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    species_properties.push_back(properties);
    i++;
  }
}

