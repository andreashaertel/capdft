// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_PROPERTIES_HPP_
#define SRC_PROPERTIES_HPP_
/** \file properties.hpp
 *  \brief Header file for the Properties class.
 *
 *  The file contains the class declarations of the Properties class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <iostream>
#include <string>
#include <unordered_map>
#include "data.hpp"
#include "template_data.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief Properties class is a container for species or system properties
 *
 */
class Properties {
 public:
  /** \brief Constructors
   *
   */
  Properties();
  Properties(const Properties& other);
  /** \brief Destructor
   *
   */
  ~Properties();
  /** \brief Removes all properties
   *
   */
  void clear();
  /** \brief Returns true if property is contained
   *
   */
  bool contains_property(const std::string& property_name) const;
  /** \brief Add a property with an arbitrary data type
   *
   */
  template<typename T>
  void add_property(
      const std::string& property_name, T property_value) {
    TemplateData<T>* new_property = new TemplateData<T>(property_value);
    properties[property_name] = new_property;
  }
  /** \brief Returns a property with an arbitrary data type
   *
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) const {
    if (contains_property(property_name)) {
      if (typeid(T) != *(properties.at(property_name)->type)) {
        std::cerr << "Properties::get_property(): \"";
        std::cerr << "Error: Requested type and property type not the same.\"";
        std::cerr << std::endl;
        exit(1);
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties.at(property_name)))->value;
      return true;
    } else {
      return false;
    }
  }

 private:
  /** \brief Contains all properties.
   *
   */
  std::unordered_map<std::string, Data*> properties;

 protected:
};
#endif  // SRC_PROPERTIES_HPP_
