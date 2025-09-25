// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file system.cpp
 *  \brief Source file for the System class.
 *
 *  The file contains the definitions of the System class.
 */
#include "system.hpp"
#include "constants.hpp"
#include "properties.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"
#include <vector>
#include "data_frame.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

// Explicit instantiation for 3d
template class System<3>;

// _____________________________________________________________________________
template <size_t dim>
System<dim>::System(Properties& properties) : Properties(properties) {
  extract_system_dimensions();
  extract_electrical_properties();
  extract_species_properties();
}
// _____________________________________________________________________________
template <size_t dim>
System<dim>::System(ParameterHandler* parameters) : Properties(parameters) {
  extract_system_dimensions();
  extract_electrical_properties();
  extract_species_properties();
}
// _____________________________________________________________________________
template <size_t dim>
System<dim>::System(Properties& properties, 
		std::vector<size_t>& species_indices) : Properties(properties) {
  extract_system_dimensions();
  extract_electrical_properties();
  extract_species_properties(species_indices);
}
// _____________________________________________________________________________
template <size_t dim>
System<dim>::System(ParameterHandler* parameters, 
		std::vector<size_t>& species_indices) : Properties(parameters) {
  extract_system_dimensions();
  extract_electrical_properties();
  extract_species_properties(species_indices);
}
// _____________________________________________________________________________
template <size_t dim>
void System<dim>::extract_system_dimensions() {
  system_lengths.clear();
  grid_counts.clear();
  periodic_boundaries.clear();
  bin_sizes.clear();
  double length;
  size_t grid_count;
  bool PBC;
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
// _____________________________________________________________________________
// Copied over from functional_es_mf::extract_electrical_properties() because
// this is also needed outside of the functional's method.
template <size_t dim>
void System<dim>::extract_electrical_properties() {
  // This variable saves a number representing, which of the el. variables are
  // available.
  int cases{0};
  // Get the variables
  try {
    if (get_property("temperature", &temperature)) {
      cases += 1;
    }
  } catch (...) {}
  try {
    if (get_property("bjerrum length", &bjerrum)) {
      cases += 2;
    }
  } catch (...) {}
  try {
    if (get_property("dielectric constant", &dielectric)) {
      cases += 4;
    }
  } catch (...) {}
  // Auxiliary variable e^2 / (4 * \epsilon_0 * \pi * k_B)
  double e2_by_4pi_boltz = 1e-3 * ELECTRON_CHARGE * ELECTRON_CHARGE /
      (4. * M_PI * VACUUM_PERMITIVITY * BOLTZMANN);
  // Calculate missing variables
  switch (cases) {
    case 3:  // dielectric constant missing
      dielectric = e2_by_4pi_boltz / (bjerrum * 1e-9 * temperature);
      break;
    case 5:  // bjerrum length missing
      bjerrum = 1e9 * e2_by_4pi_boltz / (dielectric * temperature);
      break;
    case 6:  // temperature missing
      temperature = e2_by_4pi_boltz / (dielectric * bjerrum * 1e-9);;
      break;
    case 7:  // nothing missing
      std::cerr << "extract_electrical_properties(): ";
      std::cerr << "\"Warning: All three electrical system properties are ";
      std::cerr << "specified. Consistency will not be checked.";
      std::cerr << std::endl;
      break;
    default:
      std::cerr << "extract_electrical_properties(): ";
      std::cerr << "\"Error: Specify two properties of the following: ";
      std::cerr << "temperature, bjerrum length, dielectric constant!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
template <size_t dim>
void System<dim>::extract_species_properties(
		std::vector<size_t>& species_indices) {
  // initialize
  size_t number = species_indices.size();
  species_properties.resize(number);
  affected_species_fmt.clear();
  affected_species_es.clear();
  double bulk_density;
  double diameter;
  double valency;
  // extract properties
  for (size_t i = 0; i < number; i++) {
    species_properties.at(i).clear();
    get_property("bulk density", &bulk_density, species_indices.at(i));
    species_properties.at(i).add_property("bulk density", bulk_density);
    try {
      get_property("diameter", &diameter, species_indices.at(i));
      species_properties.at(i).add_property("diameter", diameter);
      affected_species_fmt.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    try {
      get_property("valency", &valency, species_indices.at(i));
      species_properties.at(i).add_property("valency", valency);
      affected_species_es.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
  }
}
// _____________________________________________________________________________
template <size_t dim>
void System<dim>::extract_species_properties() {
  // initialize
  species_properties.clear();
  affected_species_fmt.clear();
  affected_species_es.clear();
  double bulk_density;
  double diameter;
  double valency;
  // extract properties (indices i = 0,1,... as long as bulk density is defined)
  size_t i = 0;
  Properties properties;
  while (true) {
    properties.clear();
    try {
      get_property("bulk density", &bulk_density, i);
      properties.add_property("bulk density", bulk_density);
    } catch (const Properties::MissingPropertyException*) {
      break;
    }
    try {
      get_property("diameter", &diameter, i);
      properties.add_property("diameter", diameter);
      affected_species_fmt.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    try {
      get_property("valency", &valency, i);
      properties.add_property("valency", valency);
      affected_species_es.push_back(i);
    } catch (const Properties::MissingPropertyException*) { }
    species_properties.push_back(properties);
    i++;
  }
}
// _____________________________________________________________________________
template <> // TODO: generalize
std::vector<size_t> System<3>::index_to_coordinates(size_t index) const {
  std::vector<size_t> coordinates(3);
  coordinates.at(2) = index % grid_counts.at(2);
  coordinates.at(1) = static_cast<size_t>(index / grid_counts.at(2)) % grid_counts.at(1);
  coordinates.at(0) = static_cast<size_t>(index / (grid_counts.at(2) * grid_counts.at(1)));
  return coordinates;
}
template <size_t dim>
std::vector<double> System<dim>::coordinates_to_position(
		std::vector<size_t>& coordinates) const {
  std::vector<double> pos(dim);
  for (size_t i = 0; i < dim; i++) {
    pos.at(i) = static_cast<double>(coordinates.at(i)) * bin_sizes.at(i);
  }
  return pos;
}
template <size_t dim>
std::vector<size_t> System<dim>::position_to_coordinates(
		std::vector<double>& position) const {
  std::vector<size_t> coord(dim);
  for (size_t i = 0; i < dim; i++) {
    // check whether given position is within system bounds
    if (position.at(i) < 0. || position.at(i) > system_lengths.at(i)) {
      if (periodic_boundaries.at(i)) {
	position.at(i) = std::fmod(position.at(i), system_lengths.at(i));
      } else {
        throw std::out_of_range("System<dim>::position_to_coordinates: position "
			"out of range\n");
      }
    }
    // calculate corresponding coordinates
    coord.at(i) = static_cast<size_t>(position.at(i) / bin_sizes.at(i) + 0.5);
    coord.at(i) %= grid_counts.at(i);
  }
  return coord;
}
template <size_t dim>
size_t System<dim>::coordinates_to_index(std::vector<size_t>& coordinates) const {
  size_t index = 0;
  size_t factor = 1;
  for (size_t i = dim; i > 0; i--) {
    index += coordinates.at(i - 1) * factor;
    factor *= grid_counts.at(i - 1);
  }
  return index;
}
// _____________________________________________________________________________
template <size_t dim>
size_t System<dim>::increase(size_t index, size_t direction,
		int step) const {
  if (step >= 0) {
    return increase(index, direction, static_cast<size_t>(step));
  } else {
    return decrease(index, direction, static_cast<size_t>(-step));
  }
}
template <size_t dim>
size_t System<dim>::increase(size_t index, size_t direction,
		size_t step) const {
  size_t grid_count = grid_counts.at(direction);
  if (index < grid_count - step) { return index + step; }
  else {
    if (periodic_boundaries.at(direction)) { return (index+step)%grid_count; }
    else {
      throw std::out_of_range("System<dim>::increase: index out "
		      	      "of range\n");
    }
  }
}
template <size_t dim>
size_t System<dim>::decrease(size_t index, size_t direction,
		size_t step) const {
  size_t grid_count = grid_counts.at(direction);
  if (index >= step) { return index - step; }
  else {
    if (periodic_boundaries.at(direction)) {
      return grid_count - 1 - (step-index-1)%grid_count;
    } else {
      throw std::out_of_range("System<dim>::decrease: index out "
		      	      "of range\n");
    }
  }
}
// _____________________________________________________________________________
template <>
double System<3>::interpolate(DataFrame<3, double>& values,
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
	  std::cout << "System::interpolate: cannot interpolate to position ";
	  for (double pos : neighbor_position) {
	    std::cout << pos << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }
  }
  if (weight_sum == 0.) {
    std::cerr << "System::interpolate: no interpolation points available at position ";
    for (double pos : position) {
      std::cerr << pos << " ";
    }
    std::cerr << std::endl;
    exit(1);
  }
//  // global rescaling for the interpolation weights:
//  double weight_sum = 1.;
//  for (size_t i = 0; i < dim; i++) {
//    weight_sum *= widths.at(i);
//  }
  return result / weight_sum;
}
//// _____________________________________________________________________________
//template <>
//double System<3>::laplace(DataFrame<3, double>& values,
//		  std::vector<double>& position, double mid_value) const {
//  size_t dim = 3;
//  std::vector<size_t> midpoint = position_to_coordinates(position);
//  std::vector<size_t> coordinates(dim);
//  std::vector<double> neighbor_position(dim);
//  /** Calculate usual discretized Laplacian, but with modified distances:
//   */
//  // Initialize result
//  // Calculate distances to the eight nearest neighbors
//  std::vector<std::vector<double>> distances(3);
//  std::vector<double> widths(3);
//  double interpolation_rescaling = 1.;
//  double laplacian_weight = 0.;
//  double distance;
//  for (size_t dir = 0; dir < 3; dir++) {
//    distance = std::fmod(position.at(dir), bin_sizes.at(dir));
//  //  std::cout << distance << std::endl;
//    if (distance < bin_sizes.at(dir) * 1.e-4 || // workaround (TODO)
//	bin_sizes.at(dir) - distance < bin_sizes.at(dir) * 1.e-4) {
//      distances.at(dir) = {-bin_sizes.at(dir), bin_sizes.at(dir)};
//    } else {
//      distances.at(dir) = {-distance, bin_sizes.at(dir) - distance};
//    }
//    widths.at(dir) = distances.at(dir).at(1) - distances.at(dir).at(0);
//    interpolation_rescaling /= widths.at(dir);
//    laplacian_weight += -1. / (distances.at(dir).at(0) * distances.at(dir).at(1));
//  }
//  // Sum over all neighbors with appropriate weights
//  double result = 0.;
//  double neighbor_value;
//  for (double dx : distances.at(0)) {
//    neighbor_position.at(0) = position.at(0) + dx;
//    for (double dy : distances.at(1)) {
//      neighbor_position.at(1) = position.at(1) + dy;
//      for (double dz: distances.at(2)) {
//        neighbor_position.at(2) = position.at(2) + dz;
//    	// Calculate neighboring grid point position
//        coordinates = position_to_coordinates(neighbor_position);
//	// Add neighbor contribution to result
//        neighbor_value =
//    		values.at(coordinates.at(0), coordinates.at(1), coordinates.at(2));
//        result += (widths.at(0) - std::fabs(dx)) * (widths.at(1) - std::fabs(dy)) * (widths.at(2) - std::fabs(dz))
//		* neighbor_value;
//      }
//    }
//  }
//  // Rescale interpolation weights
//  result *= interpolation_rescaling;
//  // Add midpoint contribution to result
//  result -= 2. * mid_value;
//  // Multiply all contributions with Laplacian weight
//  result *= laplacian_weight;
//  return result;
//}
// _____________________________________________________________________________
template <size_t dim>
void System<dim>::print_data(DataFrame<dim, double>& values, 
		std::ostream& out_stream) {
  std::vector<DataFrame<dim, double>> vector = {values};
  print_data(vector, out_stream);
}
template <>
void System<3>::print_data(std::vector<DataFrame<3, double>>& values,
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
template <size_t dim>
void System<dim>::print_data(DataFrame<dim, bool>& values, 
		std::ostream& out_stream) {
  std::vector<DataFrame<dim, bool>> vector = {values};
  print_data(vector, out_stream);
}
template <>
void System<3>::print_data(std::vector<DataFrame<3, bool>>& values,
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
template <>
bool System<3>::load_data(std::string filename, std::vector<size_t>& col_nums,
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
    std::cout << "load_data: file not found: " << filename << std::endl;
    return false;
  }
  while (std::getline(file, line)) {
    if (line.size() > 1 && line.at(0) != '#') { // filter
      coord = index_to_coordinates(row);
      size_t pos = 0; // position in line
      size_t col = 0; // corresponding column number
      size_t col_index = 0; // corresponding index in col_nums
			    // - assuming col_nums is sorted
      while (col <= col_nums.back()) {
	if (pos > line.size()) {
	  std::cout << "load_data: found less values than expected in data line " << row;
	  std::cout << " (" << col << " values in " << pos << " characters)\n";
	  exit(1);
	}
        if (col_nums.at(col_index) == col) { // add to *data
	  data->at(col_index).at(coord.at(0),coord.at(1),coord.at(2)) =
	  	  std::atof(line.data() + pos);
	  col_index++;
	}
	// find next value in line
	pos = line.find(' ', pos) + 1; // - syntax?
	col++;
      }
      row++;
    } else { // do nothing with empty / comment lines
      std::cout << line << std::endl;
    }
  }
  std::cout << "load_data: found " << row << " data lines.\n";
  return true;
}

//_____________________________________________________________________________
template <size_t dim>
bool System<dim>::load_data(std::string filename,
	size_t col_num, DataFrame<dim,double>* data) {
  std::vector<size_t> col_nums = {col_num};
  std::vector<DataFrame<dim,double>> data_vector = {*data};
  bool success = load_data(filename, col_nums, &data_vector);
  *data = data_vector.at(0);
  return success;
}

