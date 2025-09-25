// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface_planar.cpp
 *  \brief Source file for the BoundarySurfacePlanar class.
 *
 *  The file contains the definitions of the BoundarySurfacePlanar class.
 */
#include "boundary_surface_planar.hpp"
#include "system.hpp"
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>

// Template class instantiation of 3d version
template class BoundarySurfacePlanar<3>;

// _____________________________________________________________________________
template <size_t dim>
BoundarySurfacePlanar<dim>::BoundarySurfacePlanar(System<dim>& system, size_t index) :
	BoundarySurface<dim>(system, index) {
  extract_special_properties(system, index);
}

// _____________________________________________________________________________
template <> // could be template
void BoundarySurfacePlanar<3>::extract_special_properties(System<3>& system, 
		size_t index) {
  size_t dim = 3;
  sides.clear();
  system_heights.clear();
  size_t side;
  bool left, right;
  for (size_t dir = 0; dir < dim; dir++) {
    // determine which system walls are part of the surface: there is one
    // parameter for each axis which specifies the relevant sides (encoded in
    // binary because we cannot use vector-type properties)
    system.get_property(std::to_string(dir)+"_sides", &side, index);
    if (side > 3) {
      std::cerr << "BoundarySurfacePlanar::extract_special_properties: ";
      std::cerr << "unexpected parameter value sides = " << side << std::endl;
      exit(1);
    }
    left = (side % 2 == 1) ? true : false;
    right = (side >= 2) ? true : false;
    sides.push_back({left, right});
    // define position of the right side surface (left side is at zero)
    system_heights.push_back(
	system.system_lengths.at(dir) - system.bin_sizes.at(dir));
  }
  // periodic boundary conditions are ignored
}

// _____________________________________________________________________________
template <> // could be template
bool BoundarySurfacePlanar<3>::is_within_boundary(std::vector<double>& position)
	const {
  size_t dim = 3;
  bool result = true;
  for (size_t dir = 0; dir < dim; dir++) {
    // set result to false if the position lies on any of the surfaces
    if (sides.at(dir).at(0) && position.at(dir) <= 0.) {
      result = false;
    }
    if (sides.at(dir).at(1) && position.at(dir) >= system_heights.at(dir)) {
      result = false;
    }
  }
  return result;
}

// _____________________________________________________________________________
template <> // could be template
double BoundarySurfacePlanar<3>::distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const {
  if (direction >= 3) {
      throw(std::invalid_argument("BoundarySurfacePlanar::distance_directed: "
	"invalid direction argument " + std::to_string(direction)));
  }
  if (!is_within_boundary(position)) { return 0.; }
  if (forward) {
    if (sides.at(direction).at(1)) { // distance to right surface
      return system_heights.at(direction) - position.at(direction);
    } else { // no surface on the right
      return -1.;
    }
  } else {
    if (sides.at(direction).at(0)) { // distance to the left surface
      return position.at(direction);
    } else { // no surface on the left
      return -1.;
    }
  }
}

// _____________________________________________________________________________
template <>
void BoundarySurfacePlanar<3>::discretize_surface(std::vector<std::vector<double>>* 
		positions, double resolution) const {
  positions->clear();
  // discrete points distributed equidistantly in x and y direction
  std::vector<double> bin_sizes = {resolution, resolution, resolution};
  std::vector<size_t> grid_counts(3);
  for (size_t dir = 0; dir < 3; dir++) {
    grid_counts.at(dir) = static_cast<size_t>(system_lengths.at(dir) / resolution);
  }
  std::vector<double> position(3);
  // iterate over all three directions
  // be aware that this includes certain edge positions twice
  for (size_t dir = 0; dir < 3; dir++) {
    if (sides.at(dir).at(0) || sides.at(dir).at(1)) {
      // label the two axes orthogonal to dir by x and y (not corresponding to
      // the actual x and y axes, just aliases for a more concise code!)
      size_t dir_x = (dir + 1) % 3;
      size_t dir_y = (dir + 2) % 3;
      for (size_t i = 0; i < grid_counts.at(dir_x); i++) {
        position.at(dir_x) = bin_sizes.at(dir_x) * static_cast<double>(i);
        for (size_t j = 0; j < grid_counts.at(dir_y); j++) {
          position.at(dir_y) = bin_sizes.at(dir_y) * static_cast<double>(j);
          if (sides.at(dir).at(0)) {
	    position.at(dir) = 0.;
            positions->push_back(position);
	  }
          if (sides.at(dir).at(1)) {
	    position.at(dir) = system_heights.at(dir);
            positions->push_back(position);
	  }
	}
      }
    }
  }
}

// _____________________________________________________________________________
template <>
std::vector<double> BoundarySurfacePlanar<3>::surface_normal(
	std::vector<double>& position) const {
  // Take care with edge positions! Does recognize them as part of only one
  // surface, doesn't realize it lies in two surfaces at once.
  std::vector<double> normal = {0., 0., 0.};
  for (size_t dir = 0; dir < 3; dir++) {
    if (position.at(dir) == 0) {
      normal.at(dir) = 1.;
      return normal;
    } else if (position.at(dir) == system_heights.at(dir)) {
      normal.at(dir) = -1.;
      return normal;
    }
  }
  // If the positions lies on none of the surfaces, raise an error.
  std::cerr << "BoundarySurfacePlanar::surface_normal: given position could ";
  std::cerr << "not be associated with any surface point.\n";
  exit(1);
}

// _____________________________________________________________________________
template <>
void BoundarySurfacePlanar<3>::discretize_surface(
    std::vector<std::vector<std::vector<double>>>* distribution,
    double resolution) const {
  distribution->clear();
  std::vector<std::vector<double>> positions;
  discretize_surface(&positions, resolution);
  double area = resolution * resolution;
  std::vector<std::vector<double>> point(2);
  std::vector<double> normal(3);
  for (std::vector<double>& position : positions) {
    point.at(0) = position;
    normal = surface_normal(position);
    for (double& n : normal) {
      n *= area;
    }
    point.at(1) = normal;
    distribution->push_back(point);
  }
}

// _____________________________________________________________________________
template <>
double BoundarySurfacePlanar<3>::distance_minimal(std::vector<double>& position, 
		double upper_limit) const {
  if (!is_within_boundary(position)) {
    return 0.;
  }
  double minimum = upper_limit;
  double distance;
  for (size_t dir = 0; dir < 3; dir++) {
    if (sides.at(dir).at(0)) {
      distance = position.at(dir);
      if (distance < minimum) {
	minimum = distance;
      }
    }
    if (sides.at(dir).at(1)) {
      distance = system_heights.at(dir) - position.at(dir);
      if (distance < minimum) {
	minimum = distance;
      }
    }
  }
  return minimum;
}
// _____________________________________________________________________________
template <>
void BoundarySurfacePlanar<3>::minimal_distances(DataFrame<3, double>* result,
		double resolution, double upper_limit) {
  // Initialize
  std::vector<size_t> grid_counts = result->size_dim();
  std::vector<double> bin_sizes = bins(grid_counts);
  if (upper_limit < 0.) {
    upper_limit = 0.;
    for (double length : system_lengths) {
      upper_limit += length * length;
    }
    upper_limit = sqrt(upper_limit);
  }
  // Iterate over all grid points
  std::vector<double> position;
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      for (size_t k = 0; k != grid_counts.at(2); ++k) {
	position = {i * bin_sizes.at(0), j * bin_sizes.at(1),
		  k * bin_sizes.at(2)};
	result->at(i, j, k) = distance_minimal(position, upper_limit);
      }
    }
  }
}
