// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "properties.hpp"  // NOLINT
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
// Some template functions are implemented in the header.
// _____________________________________________________________________________
