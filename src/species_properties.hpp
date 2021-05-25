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
#include <cstddef>
#include <iostream>
#include <string>
#include <variant>
#include <vector>
#include <unordered_map>
// Class forward declarations
// _____________________________________________________________________________
/** \brief SpeciesProperties class contains all species-dependent properties
 *
 */
class SpeciesProperties {
  /** \brief This data type contains all standard types that can be used
   *
   */
  using Variant = std::variant<  // replaces typedef
      int,
      double,
      bool,
      std::string,
      std::vector<int>,
      std::vector<double>,
      std::vector<bool>,
      std::vector<std::string>>;

 public:
  /** \brief Constructor
   *
   */
  SpeciesProperties();
  explicit SpeciesProperties(size_t species_count);
  /** \brief Destructor
   *
   */
  ~SpeciesProperties();
  /** \brief Add property
   *
   */
  void add_property(
      std::string property_name,
      Variant property_values);
  /** \brief Get specific porperty
   *
   */
  bool get_property(
      const std::string& property_name, int* property_value);
  bool get_property(
      const std::string& property_name, double* property_value);
  bool get_property(
      const std::string& property_name, bool* property_value);
  bool get_property(
      const std::string& property_name, std::string* property_value);
  bool get_property(
      const std::string& property_name, std::vector<int>* property_value);
  bool get_property(
      const std::string& property_name, std::vector<double>* property_value);
  bool get_property(
      const std::string& property_name, std::vector<bool>* property_value);
  bool get_property(
      const std::string& property_name,
      std::vector<std::string>* property_value);
  /** \brief Clear all properties
   *
   */
  void clear();

 private:
  /** \brief The species properties stored in a dictionary
   *
   */
  std::unordered_map<std::string, Variant> species_properties;
  /** \brief Returns whether the speciefied property exists
   *
   */
  bool contains_property(const std::string& property_name);

 protected:
};
#endif  // SRC_SPECIES_PROPERTIES_HPP_
