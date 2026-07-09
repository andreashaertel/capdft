// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file discretization_handler.cpp
 *  \brief Source file for the DiscretizationHandler class.
 *
 *  The file contains the definitions of the DiscretizationHandler class.
 */
#include "discretization_handler.hpp"
#include "constants.hpp"
#include "properties.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"
#include <vector>
#include "data_frame.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

// Explicit instantiation for different Geometries
template class DiscretizationHandler<Cartesian>;
template class DiscretizationHandler<Planar>;
template class DiscretizationHandler<Spherical>;

// _____________________________________________________________________________
template <Geometry geometry>
DiscretizationHandler<geometry>::DiscretizationHandler(Properties& properties) : 
	Properties(properties) {
  extract_system_dimensions();
}
// _____________________________________________________________________________
template <Geometry geometry>
DiscretizationHandler<geometry>::DiscretizationHandler(ParameterHandler* parameters) : 
	Properties(parameters) {
  extract_system_dimensions();
}
// _____________________________________________________________________________
template <Geometry geometry>
void DiscretizationHandler<geometry>::extract_system_dimensions() {
  dim = dimension(geometry); 
  system_lengths.clear();
  grid_counts.clear();
  periodic_boundaries.clear();
  bin_sizes.clear();
  double length;
  size_t grid_count;
  bool PBC;
  if (dim == 1) { // parameters without indexing
    get_property("length", &length);
    get_property("grid_count", &grid_count);
    get_property("PBC", &PBC);
    system_lengths.push_back(length);
    grid_counts.push_back(grid_count);
    periodic_boundaries.push_back(PBC);
    bin_sizes.push_back(length / static_cast<double>(grid_count));
  } else { // indices indicate direction
    for (size_t dir = 0; dir < dim; dir++) {
      get_property("length", &length, dir);
      get_property("grid_count", &grid_count, dir);
      get_property("PBC", &PBC, dir);
      system_lengths.push_back(length);
      grid_counts.push_back(grid_count);
      periodic_boundaries.push_back(PBC);
      bin_sizes.push_back(length / static_cast<double>(grid_count));
    }
  }
}
// _____________________________________________________________________________
template <Geometry geometry> // TODO: generalize
std::vector<size_t> DiscretizationHandler<geometry>::index_to_coordinates(size_t index) const {
  std::cerr << "DiscretizationHandler::index_to_coordinates: This function is only implemented in 3d.\n";
  exit(1);
}
template <>
std::vector<size_t> DiscretizationHandler<Cartesian>::index_to_coordinates(size_t index) const {
  std::vector<size_t> coordinates(3);
  coordinates.at(2) = index % grid_counts.at(2);
  coordinates.at(1) = static_cast<size_t>(index / grid_counts.at(2)) % grid_counts.at(1);
  coordinates.at(0) = static_cast<size_t>(index / (grid_counts.at(2) * grid_counts.at(1)));
  return coordinates;
}
template <Geometry geometry>
std::vector<double> DiscretizationHandler<geometry>::coordinates_to_position(
		std::vector<size_t>& coordinates) const {
  std::vector<double> pos(dim);
  for (size_t i = 0; i < dim; i++) {
    pos.at(i) = static_cast<double>(coordinates.at(i)) * bin_sizes.at(i);
  }
  return pos;
}
template <Geometry geometry>
std::vector<size_t> DiscretizationHandler<geometry>::position_to_coordinates(
		std::vector<double>& position) const {
  std::vector<size_t> coord(dim);
  for (size_t i = 0; i < dim; i++) {
    // check whether given position is within discretization_handler bounds
    if (position.at(i) < 0. || position.at(i) > system_lengths.at(i)) {
      if (periodic_boundaries.at(i)) {
	position.at(i) = std::fmod(position.at(i), system_lengths.at(i));
      } else {
        throw std::out_of_range("DiscretizationHandler<geometry>::position_to_coordinates: position "
			"out of range\n");
      }
    }
    // calculate corresponding coordinates
    coord.at(i) = static_cast<size_t>(position.at(i) / bin_sizes.at(i) + 0.5);
    coord.at(i) %= grid_counts.at(i);
  }
  return coord;
}
template <Geometry geometry>
size_t DiscretizationHandler<geometry>::coordinates_to_index(std::vector<size_t>& coordinates) const {
  size_t index = 0;
  size_t factor = 1;
  for (size_t i = dim; i > 0; i--) {
    index += coordinates.at(i - 1) * factor;
    factor *= grid_counts.at(i - 1);
  }
  return index;
}
// _____________________________________________________________________________
template <Geometry geometry>
size_t DiscretizationHandler<geometry>::increase(size_t index, size_t direction,
		int step) const {
  if (step >= 0) {
    return increase(index, direction, static_cast<size_t>(step));
  } else {
    return decrease(index, direction, static_cast<size_t>(-step));
  }
}
template <Geometry geometry>
size_t DiscretizationHandler<geometry>::increase(size_t index, size_t direction,
		size_t step) const {
  size_t grid_count = grid_counts.at(direction);
  if (index < grid_count - step) { return index + step; }
  else {
    if (periodic_boundaries.at(direction)) { return (index+step)%grid_count; }
    else {
      throw std::out_of_range("DiscretizationHandler<geometry>::increase: index out "
		      	      "of range\n");
    }
  }
}
template <Geometry geometry>
size_t DiscretizationHandler<geometry>::decrease(size_t index, size_t direction,
		size_t step) const {
  size_t grid_count = grid_counts.at(direction);
  if (index >= step) { return index - step; }
  else {
    if (periodic_boundaries.at(direction)) {
      return grid_count - 1 - (step-index-1)%grid_count;
    } else {
      throw std::out_of_range("DiscretizationHandler<geometry>::decrease: index out "
		      	      "of range\n");
    }
  }
}
// _____________________________________________________________________________
template <Geometry geometry>
double DiscretizationHandler<geometry>::interpolate(
		DataFrame<dimension(geometry), double>& values,
	        std::vector<double>& position) {
  std::cerr << "DiscretizationHandler::interpolate: This function is only implemented in 3d.\n";
  exit(1);
}
template <>
double DiscretizationHandler<Cartesian>::interpolate(DataFrame<3, double>& values,
		  std::vector<double>& position) {
  // Initialize: dimension, result, mid- & interpolation-point coordinates
  size_t dim = 3;
  double result = 0.;
  std::vector<size_t> coordinates(dim);
  std::vector<double> neighbor_position(dim);
  // Find all neighboring grid points (usually 2^dim points, maybe less)
  std::vector<std::vector<double>> distances(dim);
  std::vector<double> widths(dim);
  double distance;
  for (size_t dir = 0; dir < dim; dir++) {
    /* The distances to the neighboring grid planes can be determined by taking
     * the modulus w.r.t. the bin size:
     */
    distance = std::fmod(position.at(dir), bin_sizes.at(dir));
    if (distance < bin_sizes.at(dir) * 1.e-4 || // workaround (TODO)
	bin_sizes.at(dir) - distance < bin_sizes.at(dir) * 1.e-4) {
      /* When the given position lies on a grid plane, only interpolate between
       * points within this plane. I.e. use zero offset in this direction.
       */
      distances.at(dir) = {0.};
    } else {
      /* When the given position lies between two grid planes, interpolate
       * between points on those two planes.
       */
      distances.at(dir) = {-distance, bin_sizes.at(dir) - distance};
    }
    widths.at(dir) = bin_sizes.at(dir);
  }
  // Sum over all these neighboring points with weights depending on their
  // distance to the given position
  double neighbor_value;
  double weight;
  double weight_sum = 0.;
  for (double dx : distances.at(0)) {
    neighbor_position.at(0) = position.at(0) + dx;
    for (double dy : distances.at(1)) {
      neighbor_position.at(1) = position.at(1) + dy;
      for (double dz: distances.at(2)) {
        neighbor_position.at(2) = position.at(2) + dz;
	try {
    	  // Calculate neighboring grid point position
          coordinates = position_to_coordinates(neighbor_position);
          neighbor_value = values.at(coordinates.at(0), coordinates.at(1), coordinates.at(2));
	  weight = (widths.at(0) - std::fabs(dx)) *
		   (widths.at(1) - std::fabs(dy)) *
		   (widths.at(2) - std::fabs(dz));
          result += neighbor_value * weight;
	  weight_sum += weight;
	} catch (...) { //const std::out_of_range*
	  std::cout << "DiscretizationHandler::interpolate: cannot interpolate to position ";
	  for (double pos : neighbor_position) {
	    std::cout << pos << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
  if (weight_sum == 0.) {
    std::cerr << "DiscretizationHandler::interpolate: no interpolation points available at position ";
    for (double pos : position) {
      std::cerr << pos << " ";
    }
    std::cerr << std::endl;
    exit(1);
  }
  return result / weight_sum;
}
// _____________________________________________________________________________
template <Geometry geometry>
double DiscretizationHandler<geometry>::integrate(
		DataFrame<dimension(geometry), double>& values) {
  std::cerr << "DiscretizationHandler::integrate: No implementation available ";
  std::cerr << "in geometry " << geometry << ".\n";
  exit(1);
}
// Specialization for different geometries: Refer to the implementation in
// integration.hpp.
template <>
double DiscretizationHandler<Cartesian>::integrate(
		DataFrame<3, double>& values) {
  return integration_3d_closed(values, bin_sizes);
}
template <>
double DiscretizationHandler<Planar>::integrate(
		DataFrame<1, double>& values) {
  return integration_1d_closed(values, bin_sizes.at(0));
}
template <>
double DiscretizationHandler<Spherical>::integrate(
		DataFrame<1, double>& values) {
  return integration_1d_radial_open_closed(values, bin_sizes.at(0));
}
// _____________________________________________________________________________
template <Geometry geometry>
void DiscretizationHandler<geometry>::print_data(DataFrame<dimension(geometry), double>& values, 
		std::ostream& out_stream) {
  std::vector<DataFrame<dimension(geometry), double>> vector = {values};
  print_data(vector, out_stream);
}
template <Geometry geometry>
void DiscretizationHandler<geometry>::print_data(std::vector<DataFrame<dimension(geometry), double>>& values,
		std::ostream& out_stream) {
  std::cerr << "DiscretizationHandler::print_data: This function is only implemented in 3d.\n";
  exit(1);
}
template <>
void DiscretizationHandler<Cartesian>::print_data(std::vector<DataFrame<3, double>>& values,
		std::ostream& out_stream) {
  // print grid dimensions as header
  out_stream << "# grid_counts=";
  for (size_t gc : grid_counts) {
    out_stream << gc << ",";
  }
  out_stream << std::endl;
  // print data
  std::vector<double> position(3);
  std::vector<size_t> coordinates(3);
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
	coordinates = {i, j, k};
	position = coordinates_to_position(coordinates);
        out_stream << position.at(0) << " " << position.at(1) << " " << 
		position.at(2) << " ";
	for (size_t n = 0; n < values.size(); n++) {
          out_stream << values.at(n).at(i, j, k) << " ";
	}
        out_stream << std::endl;
      }
    }
  }
}
template <Geometry geometry>
void DiscretizationHandler<geometry>::print_data(DataFrame<dimension(geometry), bool>& values, 
		std::ostream& out_stream) {
  std::vector<DataFrame<dimension(geometry), bool>> vector = {values};
  print_data(vector, out_stream);
}
template <Geometry geometry>
void DiscretizationHandler<geometry>::print_data(std::vector<DataFrame<dimension(geometry), bool>>& values,
		std::ostream& out_stream) {
  std::cerr << "DiscretizationHandler::print_data: This function is only implemented in 3d.\n";
  exit(1);
}
template <>
void DiscretizationHandler<Cartesian>::print_data(std::vector<DataFrame<3, bool>>& values,
		std::ostream& out_stream) {
  // print grid dimensions as header
  out_stream << "# grid_counts=";
  for (size_t gc : grid_counts) {
    out_stream << gc << ",";
  }
  out_stream << std::endl;
  // print data
  std::vector<double> position(3);
  std::vector<size_t> coordinates(3);
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
	coordinates = {i, j, k};
	position = coordinates_to_position(coordinates);
        out_stream << position.at(0) << " " << position.at(1) << " " << 
		position.at(2) << " ";
	for (size_t n = 0; n < values.size(); n++) {
	  if (values.at(n).at(i, j, k)) {
            out_stream << 1 << " ";
	  } else {
            out_stream << 0 << " ";
	  }
	}
        out_stream << std::endl;
      }
    }
  }
}
//_____________________________________________________________________________
template <Geometry geometry>
bool DiscretizationHandler<geometry>::load_data(std::string filename, std::vector<size_t>& col_nums,
	std::vector<DataFrame<dimension(geometry),double>>* data) {
  std::cerr << "DiscretizationHandler::load_data: This function is only implemented in 3d.\n";
  exit(1);
}
template <>
bool DiscretizationHandler<Cartesian>::load_data(std::string filename, std::vector<size_t>& col_nums,
	std::vector<DataFrame<3,double>>* data) {
  const size_t dim = 3;
  std::cout << "load data from " << filename << std::endl;
  for (size_t gc : grid_counts) { std::cout << gc << ","; }
  std::cout << std::endl;
  std::string line;
  size_t row = 0; // line number (ignoring empty lines and comments)
  std::vector<size_t> coord(dim); // corresponding coordinate indices i,j,k
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "DiscretizationHandler::load_data: file not found: " << filename << std::endl;
    return false;
  }
  while (std::getline(file, line)) {
    if (line.size() > 1 && line.at(0) != '#') { // filter out empty lines / comments
      coord = index_to_coordinates(row);
      size_t pos = 0; // position in line
      size_t col = 0; // corresponding column number
      size_t col_index = 0; // corresponding index in col_nums
			    // - assuming col_nums is sorted
      while (col <= col_nums.back()) {
	if (pos > line.size()) {
	  std::cout << "DiscretizationHandler::load_data: found less values than expected in data line " << row;
	  std::cout << " (" << col << " values in " << pos << " characters)\n";
	  exit(1);
	}
        if (col_nums.at(col_index) == col) { // add to *data
	  data->at(col_index).at(coord.at(0),coord.at(1),coord.at(2)) =
	  	  std::atof(line.data() + pos);
	  col_index++;
	}
	// find next value in line
	pos = line.find(' ', pos) + 1;
	col++;
      }
      row++;
    } else { // do nothing with empty / comment lines
      std::cout << line << std::endl;
    }
  }
  std::cout << "DiscretizationHandler::load_data: found " << row << " data lines.\n";
  return true;
}

//_____________________________________________________________________________
template <Geometry geometry>
bool DiscretizationHandler<geometry>::load_data(std::string filename,
	size_t col_num, DataFrame<dimension(geometry),double>* data) {
  std::vector<size_t> col_nums = {col_num};
  std::vector<DataFrame<dimension(geometry),double>> data_vector = {*data};
  bool success = load_data(filename, col_nums, &data_vector);
  *data = data_vector.at(0);
  return success;
}

