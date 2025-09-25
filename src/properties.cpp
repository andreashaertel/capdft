// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <algorithm>
#include "properties.hpp"  // NOLINT
#include "../../parameter_handler/src/parameter_handler.hpp"
// _____________________________________________________________________________
Properties::Properties() {
  clear();
}
// _____________________________________________________________________________
Properties::Properties(ParameterHandler* parameters) {
  clear();
  parameter_handler = parameters;
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
//// _____________________________________________________________________________
//void Properties::extract_system_properties(std::string& filename) {
//  // extract system lengths, PBC, bjerrum or similar, grid counts
//  // calculate remaining bjerrum or similar
//}
//// _____________________________________________________________________________
//void Properties::extract_species_properties(std::string& filename,
//		size_t index) {
//  // extract diameter, valency, bulk density of species <index>
//}
//// _____________________________________________________________________________
//void Properties::extract_properties(std::string& filename) {
//  // decide on data type: check whether all are numbers or there's a comma or
//  // there are letters...
//}
//template <typename T>
//bool Properties::extract_property<T>(const std::string& property_name) {
//  std::cerr << "Properties::extract_property: invalid data type\n";
//  exit(1);
//  return false;
//}
template <>
bool Properties::extract_property<double>(const std::string& property_name, 
		double* property_value) {
  try {
    *property_value = parameter_handler->get_double(property_name);
    return true;
  } catch (const ParameterHandler::BadParamException*) {
    return false;
  }
}
template <>
bool Properties::extract_property<int>(const std::string& property_name,
		int* property_value) {
  try {
    *property_value = parameter_handler->get_int(property_name);
    return true;
  } catch (const ParameterHandler::BadParamException*) {
    return false;
  }
}
template <>
bool Properties::extract_property<bool>(const std::string& property_name,
		bool* property_value) {
  try {
    *property_value = parameter_handler->get_bool(property_name);
    return true;
  } catch (const ParameterHandler::BadParamException*) {
    return false;
  }
}
template <>
bool Properties::extract_property<std::string>(const std::string& property_name,
		std::string* property_value) {
  try {
    *property_value = parameter_handler->get_string(property_name);
    return true;
  } catch (const ParameterHandler::BadParamException*) {
    return false;
  }
}
template <>
bool Properties::extract_property<size_t>(const std::string& property_name,
		size_t* property_value) {
  try {
    int value = parameter_handler->get_int(property_name);
    if (value < 0) {
      std::cerr << "Properties::extract_property<size_t>: requested parameter "
	      "has negative value, cannot be converted to size_t\n";
      exit(1);
      return false;
    } else {
      *property_value = static_cast<size_t>(value);
      return true;
    }
  } catch (const ParameterHandler::BadParamException*) {
    return false;
  }
}
// _____________________________________________________________________________
std::string Properties::indexed(const std::string& property_name, size_t index)
	const {
  return std::to_string(index) + "_" + property_name;
}
// _____________________________________________________________________________
std::string Properties::parameter_name(const std::string& property_name) const {
  std::string alternative = property_name;
  // replace all spaces by underscores
  std::ranges::replace(alternative.begin(), alternative.end(), ' ', '_');
  return alternative;
}
// _____________________________________________________________________________
// Some template functions are implemented in the header.
// _____________________________________________________________________________
