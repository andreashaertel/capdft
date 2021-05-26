/*******************************************************************************
* Copyright 2021 Moritz Bueltmann
* Authors: Moritz Bueltmann <moritz.bueltmann@gmx.de>
* Physics Department Albert-Ludwigs-Universitaet
*******************************************************************************/
#ifndef SRC_SPECIES_PROPERTIES_HPP_
#define SRC_SPECIES_PROPERTIES_HPP_
/** \file species_properties.hpp
 *  \brief Header file for the SpeciesProperties class.
 *
 *  The file contains the class declarations of the SpeciesProperties class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <unordered_map>
#include "data.hpp"
#include "template_data.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief SpeciesProperties class contains all species-dependent properties
 *
 */
class SpeciesProperties {
 public:
  /** \brief Constructors
   *
   */
  SpeciesProperties();
  /** \brief Destructor
   *
   */
  ~SpeciesProperties();
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
    species_properties[property_name] = new_property;
  }
  /** \brief Returns a property with an arbitrary data type
   *
   */
  template<typename T>
  bool get_property(const std::string& property_name, T* property_value) {
    if (contains_property(property_name)) {
      if (typeid(T) != *(species_properties[property_name]->type)) {
        exit(1);  // TODO: throw error
      }
      *property_value = (dynamic_cast<TemplateData<T>*>(
          species_properties[property_name]))->value;
      return true;
    } else {
      return false;
    }
  }

 private:
  /** \brief Contains all species properties.
   *
   */
  std::unordered_map<std::string, Data*> species_properties;

 protected:
};
#include "species_properties.hpp"
#endif  // SRC_SPECIES_PROPERTIES_HPP_
