// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "system.hpp"  // NOLINT
#include <stdexcept>
// Class type declarations
template class System<DataFrame<1, double>>;
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
System<T>::System() {
  //
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
System<T>::System(
    Properties system_properties,
    std::vector<Properties> species_properties)
  : system_properties(system_properties),
    species_properties(species_properties) {
  // Get the number of grid points that were specified in system_properties
  size_t grid_count{0};
  system_properties.get_property("grid count", &grid_count);
  // From the system properties create the density Profile
  for (size_t i = 0; i < species_properties.size(); i++) {
    density_profiles.push_back(T(grid_count));
  }
  set_bulk_densities();
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
System<T>::~System() {
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
const Properties& System<T>::get_system_properties() const {
  return system_properties;
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
const std::vector<Properties>& System<T>::get_species_properties() const {
  return species_properties;
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
const std::vector<T>* System<T>::get_density_profiles_pointer() const {
  return &density_profiles;
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
std::vector<T>& System<T>::get_density_profiles() {
  return density_profiles;
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
void System<T>::update_density_profiles(
    const std::vector<T>& other_density_profiles) {
  // Check for correct size of fugacities in comparison to the number of
  // species.
  if (this->density_profiles.size() != other_density_profiles.size())
    throw std::length_error("Size of vector density_profiles does not agree with number of species.");  // NOLINT
  // Iterate through all species and update the density profiles
  auto it_other = other_density_profiles.begin();
  for (auto it_this = this->density_profiles.begin();
      it_this != this->density_profiles.end(); ++it_this) {
    *it_this = *it_other;
    ++it_other;
  }
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
void System<T>::set_bulk_densities() {
  double bulk_density{0.};
  size_t grid_count{0};
  system_properties.get_property("grid count", &grid_count);
  for (auto it = species_properties.begin(); it != species_properties.end();
      ++it) {
    it->get_property("bulk density", &bulk_density);
    for (size_t i = 0; i != grid_count; ++i) {
      density_profiles.at(it - species_properties.begin()).at(i) = bulk_density;
    }
  }
}
// _____________________________________________________________________________
template <typename T>  // Template for different data frames DF*
void System<T>::set_fugacities(const std::vector<double>& fugacities) {
  // TODO(Andreas): input-variablen -> Referenzen
  // TODO(Andreas): output-variablen -> Pointer
  // Check for correct size of fugacities in comparison to the number of
  // species.
  if (species_properties.size() != fugacities.size()) {
    throw std::length_error("Size of vector chempot does not agree with number of species.");  // NOLINT
  }
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // TODO(Andreas): It would be more save to give an index to each species and
  //                to connect the fugacities to these indices.
  // AH, 24.02.2022
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Iterate through all species and update the fugacities
  // std::vector<double>::iterator it2 = fugacities.begin();
  // for (std::vector<Properties>::iterator it = species_properties.begin();
  //      it != species_properties.end(); ++it) {
  //  it->update_property<double>("fugacity", *it2);
  //  ++it2;
  //}
}
// _____________________________________________________________________________
