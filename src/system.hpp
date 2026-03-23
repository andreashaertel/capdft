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
#include "properties.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"
#include <vector>
#include "data_frame.hpp"
#include <iostream>
#include <fstream>

/** \brief This class handles the system properties.
 * 
 * This class explicitly defines the system dimensions, electrical properties
 * and species properties. It also has methods to handle discretization
 * issues (e.g. conversion between DataFrame-index and corresponding position).
 * Finally, it can implicitly contain any number of other parameters which can
 * be read out using the get_property function.
 */
template <size_t dim>
class System : public Properties {
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
  /** \brief Conversion between different coordinate formats
   *
   * \param position (std::vector<double>): physical position, same unit as system_lengths
   * \param coordinates (std::vector<size_t>): discretized coordinates (DataFrame indices)
   * \param index (size_t): one-dimensional super-index (vector / flattened out DataFrame)
   *
   * (TODO?: These functions could take PBC into account if the given coordinate or
   * position is out of system bounds (and throw an std::out_of_range-exception
   * if it is out of bounds without PBC).)
   */
  std::vector<size_t> index_to_coordinates(size_t index) const;
  std::vector<double> coordinates_to_position(std::vector<size_t>& coordinates) const;
  std::vector<size_t> position_to_coordinates(std::vector<double>& position) const;
  size_t coordinates_to_index(std::vector<size_t>& coordinates) const;
  // ...
  /** \brief In-/decrease coordinates, taking PBC into account
   *
   * \param index: DataFrame index in the given direction
   * \param direction: along which axis (x=0, y=1, ...) to take the step
   * \param step: stepsize (Negative step in increase means decrease by the absolute value.)
   */
  size_t increase(size_t index, size_t direction, size_t step) const;
  size_t decrease(size_t index, size_t direction, size_t step) const;
  size_t increase(size_t index, size_t direction, int step) const;
  /** \brief Interpolate a field (DataFrame) of values at a given position
   */
  double interpolate(DataFrame<dim, double>& values,
		  std::vector<double>& position);
  /** \brief Print DataFrame to specified output
   */
  void print_data(DataFrame<dim, double>& values, std::ostream& out);
  void print_data(std::vector<DataFrame<dim, double>>& values, std::ostream& out);
  void print_data(DataFrame<dim, bool>& values, std::ostream& out);
  void print_data(std::vector<DataFrame<dim, bool>>& values, std::ostream& out);
  /** \brief Get DataFrame from data file (multiple columns)
   *
   * Return value says whether the process was successful.
   *
   * \param col_nums: indices of the relevant data columns 
   * \param data: pointer to DataFrame where the extracted values shall be stored
   */
  bool load_data(std::string filename, std::vector<size_t>& col_nums,
	std::vector<DataFrame<dim,double>>* data);
  /** \brief Get DataFrame from data file (single column)
   */
  bool load_data(std::string filename,
	size_t col_num, DataFrame<dim,double>* data);
  /** \brief System dimensions
   *
   * Note that the grid could be defined in slightly different ways. Here we
   * set, in each direction, the first grid point (index 0) at position 0 and
   * the last one (index grid_counts - 1) at position system_lengths - bin_sizes.
   */
  std::vector<double> bin_sizes;
  std::vector<size_t> grid_counts;
  std::vector<double> system_lengths;
  std::vector<bool> periodic_boundaries;
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
private:
  /** \brief Explicitly initialize the system properties as member variables
   */
  void extract_system_dimensions();
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
