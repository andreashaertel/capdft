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
#include <stdexcept>
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
  bool contains_property(const std::string& property_name);
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
  bool get_property(const std::string& property_name, T* property_value) {
    if (contains_property(property_name)) {
      if (typeid(T) != *(properties[property_name]->type)) {
        throw std::invalid_argument(
            "Requested type and property type not the same.");
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          properties[property_name]))->value;
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
