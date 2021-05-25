/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
#include "species_properties.hpp"
// _____________________________________________________________________________
SpeciesProperties::SpeciesProperties() {
  //
}
// _____________________________________________________________________________
SpeciesProperties::SpeciesProperties(size_t species_count) {
  //
}
// _____________________________________________________________________________
SpeciesProperties::~SpeciesProperties() {
  //
}
// _____________________________________________________________________________
void SpeciesProperties::add_property(
    std::string property_name, Variant property_values
) {
  // Add name-value pair
  species_properties[property_name] = property_values;
}
// _____________________________________________________________________________
bool SpeciesProperties::contains_property(const std::string& property_name) {
  auto search = species_properties.find(property_name);
  if (search != species_properties.end()) {
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, int* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<int>(species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, double* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<double>(species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, bool* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<bool>(species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, std::string* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<std::string>(species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, std::vector<int>* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<std::vector<int>>(
        species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, std::vector<double>* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<std::vector<double>>(
        species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, std::vector<bool>* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<std::vector<bool>>(
        species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
bool SpeciesProperties::get_property(
    const std::string& property_name, std::vector<std::string>* property_value
) {
  if (contains_property(property_name)) {
    *property_value = std::get<std::vector<std::string>>(
        species_properties[property_name]);
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
void SpeciesProperties::clear() {
  species_properties.clear();
}
// _____________________________________________________________________________
