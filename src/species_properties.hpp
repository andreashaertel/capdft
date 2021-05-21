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
  using VectorVariant = std::variant<  // replaces typedef
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
      VectorVariant property_values);
  /** \brief Get specific porperty
   *
   */
  void get_property(
      const std::string& property_name, std::vector<int>* poperty_values);
  void get_property(
      const std::string& property_name, std::vector<double>* poperty_values);
  void get_property(
      const std::string& property_name, std::vector<bool>* poperty_values);
  void get_property(
      const std::string& property_name,
      std::vector<std::string>* poperty_values);

 private:
  /** \brief Number of species
   *
   *  This variable keeps track of the species count, such that all properties
   *  have the same number of entries.
   */
  size_t species_count;
  /** \brief Stores the species properties in a dictionary
   *
   */
  std::unordered_map<std::string, VectorVariant> species_properties;
  /** \brief Calculates the species count by getting the vector size
   *
   */
  size_t get_vector_size(const VectorVariant& v);

 protected:
};
#endif  // SRC_SPECIES_PROPERTIES_HPP_
