/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
#include "species_properties.hpp"
// _____________________________________________________________________________
SpeciesProperties::SpeciesProperties()
    : species_count(0) {
  //
}
// _____________________________________________________________________________
SpeciesProperties::SpeciesProperties(size_t species_count)
    : species_count(species_count) {
  //
}
// _____________________________________________________________________________
SpeciesProperties::~SpeciesProperties() {
  //
}
// _____________________________________________________________________________
void SpeciesProperties::add_property(
    std::string property_name, VectorVariant property_values
) {
  size_t vector_count = get_vector_size(property_values);
  // Check for special cases
  if (species_count == 0) {  // species_count has not been determined yet
    species_count = vector_count;
  } else if (species_count != vector_count) {  // invalid # species
    std::cerr << species_count << std::endl;
    // TODO(Moritz): throw error
    exit(1);
  }
  // Add name-value pair
  species_properties[property_name] = property_values;
}
// _____________________________________________________________________________
void SpeciesProperties::get_property(
    const std::string& property_name, std::vector<int>* poperty_values
) {
  *poperty_values =
      std::get<std::vector<int>>(species_properties[property_name]);
}
// _____________________________________________________________________________
void SpeciesProperties::get_property(
    const std::string& property_name, std::vector<double>* poperty_values
) {
  *poperty_values =
      std::get<std::vector<double>>(species_properties[property_name]);
}
// _____________________________________________________________________________
void SpeciesProperties::get_property(
    const std::string& property_name, std::vector<bool>* poperty_values
) {
  *poperty_values =
      std::get<std::vector<bool>>(species_properties[property_name]);
}
// _____________________________________________________________________________
void SpeciesProperties::get_property(
    const std::string& property_name, std::vector<std::string>* poperty_values
) {
  *poperty_values =
      std::get<std::vector<std::string>>(species_properties[property_name]);
}
// _____________________________________________________________________________
size_t SpeciesProperties::get_vector_size(const VectorVariant& v) {
  try {
    return std::get<std::vector<int>>(v).size();
  } catch (...) {}
  try {
    return std::get<std::vector<double>>(v).size();
  } catch (...) {}

  try {
    return std::get<std::vector<bool>>(v).size();
  } catch (...) {}
  try {
    return std::get<std::vector<std::string>>(v).size();
  } catch (...) {
    // TODO(Moritz): throw
    exit(1);
  }
}
// _____________________________________________________________________________
