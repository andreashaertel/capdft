// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DISCRETIZATION_HANDLER_HPP_
#define SRC_DISCRETIZATION_HANDLER_HPP_
/** \file discretization_handler.hpp
 * \brief Header file for the DiscretizationHandler class.
 *
 *  The file contains the declarations of the DiscretizationHandler class.
 */
// Includes
#include "../../parameter_handler/src/parameter_handler.hpp"
#include "properties.hpp"
#include "data_frame.hpp"
#include "integration.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

/** \brief Different geometries
 *
 * Depending on the system's symmetries, different geometries can be chosen:
 * - Cartesian: description in terms of 3d-cartesian coordinates x,y,z
 * - Planar: description in terms of a single coordinate z
 * - Spherical: description in terms of a single coordinate r
 */
enum Geometry {Cartesian, Planar, Spherical};
/** \brief Define effective dimension for each geometry
 */
constexpr size_t dimension(Geometry geometry) {
  switch (geometry) {
    case Cartesian: return 3;
    case Planar: return 1;
    case Spherical: return 1;
  }
}

/** \brief This class provides tools for handling the spatial discretization in
 * 	different geometries.
 * 
 * This includes for instance
 * - the conversion between different coordinate formats (DataFrame index,
 * 	physical position, ...),
 * - the identification of neighboring grid points
 * - in-/output of data from/into data files
 * - numerical integration, interpolation or differentiation of data
 *
 * The DiscretizationHandler is initialized based on a Properties or
 * ParameterHandler object, which must specify the properties length, bin_size, 
 * grid_count for each of the coordinate axes, which are labeled by the indices
 * i=0,...,dimension-1.
 *
 * Technical remark: We use here a multi-level inheritance chain 
 * Properties <- DiscretizationHandler <- System,
 * so that each of these classes provides access to all defined properties via
 * (Properties::)get_property(...).
 */
template <Geometry geometry>
class DiscretizationHandler : public Properties {
public:
  /** \brief Empty Constructor
   */
  DiscretizationHandler() {};
  /** \brief Automatic constructor with Properties object
   */
  DiscretizationHandler(Properties& properties);
  /** \brief Automatic constructor with ParameterHandler
   */
  DiscretizationHandler(ParameterHandler* parameter_handler);
  /** \brief Conversion between different coordinate formats
   *
   * \param position (std::vector<double>): physical position, same unit as discretization_handler_lengths
   * \param coordinates (std::vector<size_t>): discretized coordinates (DataFrame indices)
   * \param index (size_t): one-dimensional super-index (vector / flattened out DataFrame)
   *
   * (TODO?: These functions could take PBC into account if the given coordinate or
   * position is out of discretization_handler bounds (and throw an std::out_of_range-exception
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
  double interpolate(DataFrame<dimension(geometry), double>& values,
		  std::vector<double>& position);
  /** \brief Integrate a field (DataFrame) over the system volume
   */
  double integrate(DataFrame<dimension(geometry), double>& values); 
  /** \brief Print DataFrame to specified output
   */
  void print_data(DataFrame<dimension(geometry), double>& values, std::ostream& out);
  void print_data(std::vector<DataFrame<dimension(geometry), double>>& values, std::ostream& out);
  void print_data(DataFrame<dimension(geometry), bool>& values, std::ostream& out);
  void print_data(std::vector<DataFrame<dimension(geometry), bool>>& values, std::ostream& out);
  /** \brief Get DataFrame from data file (multiple columns)
   *
   * Return value says whether the process was successful.
   *
   * \param col_nums: indices of the relevant data columns 
   * \param data: pointer to DataFrame where the extracted values shall be stored
   */
  bool load_data(std::string filename, std::vector<size_t>& col_nums,
	std::vector<DataFrame<dimension(geometry),double>>* data);
  /** \brief Get DataFrame from data file (single column)
   */
  bool load_data(std::string filename,
	size_t col_num, DataFrame<dimension(geometry),double>* data);
  /** \brief System dimensions
   *
   * Note that the grid could be defined in slightly different ways. Here we
   * set, in each direction, the first grid point (index 0) at position 0 and
   * the last one (index grid_counts - 1) at position system_lengths - bin_sizes.
   */
  size_t dim;
  std::vector<double> bin_sizes;
  std::vector<size_t> grid_counts;
  std::vector<double> system_lengths;
  std::vector<bool> periodic_boundaries;
private:
  /** \brief Explicitly initialize the required properties as member variables
   */
  void extract_system_dimensions();
};
#endif  // SRC_DISCRETIZATION_HANDLER_HPP_
