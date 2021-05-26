/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
#include "species_properties.hpp"
// _____________________________________________________________________________
SpeciesProperties::SpeciesProperties() {
  clear();
}
// _____________________________________________________________________________
SpeciesProperties::~SpeciesProperties() {
  //
}
// _____________________________________________________________________________
void SpeciesProperties::clear() {
  species_properties.clear();
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
