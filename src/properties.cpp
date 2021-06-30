// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "properties.hpp"
// _____________________________________________________________________________
Properties::Properties() {
  clear();
}
// _____________________________________________________________________________
Properties::Properties(const Properties& other) {
  clear();
  properties = other.properties;
}
// _____________________________________________________________________________
Properties::~Properties() {
  //
}
// _____________________________________________________________________________
void Properties::clear() {
  properties.clear();
}
// _____________________________________________________________________________
bool Properties::contains_property(const std::string& property_name) const {
  auto search = properties.find(property_name);
  if (search != properties.end()) {
    return true;
  } else {
    return false;
  }
}
// _____________________________________________________________________________
