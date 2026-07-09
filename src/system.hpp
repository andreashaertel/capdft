// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_SYSTEM_HPP_
#define SRC_SYSTEM_HPP_
/** \file system.hpp
 * \brief Header file for the System class.
 *
 *  The file contains the declarations of the System class.
 */
// Includes
#include "../../parameter_handler/src/parameter_handler.hpp"
#include "properties.hpp"
#include "data_frame.hpp"
#include "discretization_handler.hpp"
#include <vector>
#include <iostream>
#include <fstream>

/** \brief This class handles the system properties.
 * 
 * This class explicitly defines the system dimensions, electrical properties
 * and species properties. Derived from the DiscretizationHandler, it also has
 * methods to handle discretization issues.
 * Finally, it can implicitly contain any number of other parameters which can
 * be read out using the get_property function.
 */
template <Geometry geometry>
class System : public DiscretizationHandler<geometry> {
public:
  /** \brief Empty Constructor
   */
  System() {};
  /** \brief Automatic constructor with Properties object
   * 
   * If no species indices are specified, all existing species are read out
   * (as long as they are labeled by consecutive indices starting at 0).
   */
  System(Properties& properties);
  /** \brief Manual constructor with Properties object
   *
   * If the species indices are specified in the constructor, be aware that the
   * System class will
   * internally relabel them in consecutive order (i.e. the affected_species 
   * and the species_indices are in general based on a different labeling system than the
   * species_indices handed to the constructor). Therefore, don't access the
   * System::species_properties with the same indices you used to use, but just
   * with indices 0,1,2,...
   */
  System(Properties& properties, std::vector<size_t>& species_indices);
  /** \brief Automatic constructor with ParameterHandler
   */
  System(ParameterHandler* parameter_handler);
  /** \brief Manual constructor with ParameterHandler
   */
  System(ParameterHandler* parameter_handler, std::vector<size_t>& species_indices);
  /** \brief Electrical properties
   */
  double bjerrum;
  double dielectric;
  double temperature;
  /** \brief Species properties as vector
   */
  std::vector<Properties> species_properties;
  /** \brief Hard sphere species
   */
  std::vector<size_t> affected_species_fmt;
  /** \brief Electrically charged species
   */
  std::vector<size_t> affected_species_es;
  // System dimensions and discretization tools are inherited from
  // DiscretizationHandler. 
private:
  /** \brief Explicitly initialize the system properties as member variables
   */
  void extract_electrical_properties();
  /** \brief Construct a vector of Properties objects for the different species
   * of particles in the system
   *
   * \param species_indices: indices of the required species
   */
  void extract_species_properties();
  void extract_species_properties(std::vector<size_t>& species_indices);
};
#endif  // SRC_SYSTEM_HPP_
