// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "system.hpp"
// _____________________________________________________________________________
//System::System() {
//  //
//}
// _____________________________________________________________________________
System::System(
    Properties system_properties,
    std::vector<Properties> species_properties)
  : system_properties(system_properties),
    species_properties(species_properties) {
  // Get the number of grid points that were specified in system_properties
  size_t grid_count{0};
  system_properties.get_property("grid count", &grid_count);
  // From the system properties create the density Profile
  density_profile = new DataField<double>(
      species_properties.size(), grid_count);
  bulk();
}
// _____________________________________________________________________________
System::~System() {
  density_profile->~DataField();
}
// _____________________________________________________________________________
const Properties& System::get_system_properties() const {
  return system_properties;
}
// _____________________________________________________________________________
const std::vector<Properties>& System::get_species_properties() const {
  return species_properties;
}
// _____________________________________________________________________________
DataField<double>* System::get_density_profile_pointer() {
  return density_profile;
}
// _____________________________________________________________________________
void System::bulk() {
  double bulk_density{0.};
  size_t grid_count{0};
  system_properties.get_property("grid count", &grid_count);
  for (auto it = species_properties.begin(); it != species_properties.end();
      ++it) {
    it->get_property("bulk density", &bulk_density);
    for (size_t i = 0; i != grid_count; ++i) {
      density_profile->at(it - species_properties.begin(), i) = bulk_density;
    }
  }
}
// _____________________________________________________________________________
