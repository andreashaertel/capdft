// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface_sphere.cpp
 *  \brief Source file for the BoundarySurfaceSphere class.
 *
 *  The file contains the definitions of the BoundarySurfaceSphere class.
 */
#include "boundary_surface_sphere.hpp"
#include "system.hpp"
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>

// _____________________________________________________________________________
BoundarySurfaceSphere::BoundarySurfaceSphere(System<3>& system, size_t index) :
	BoundarySurface<3>(system, index) {
  extract_special_properties(system, index);
}

// _____________________________________________________________________________
void BoundarySurfaceSphere::extract_special_properties(System<3>& system, 
		size_t index) {
  system.get_property("radius", &radius, index);
  try {
    system.get_property("x_mid", &x_mid, index);
  } catch (const Properties::MissingPropertyException*) {
    x_mid = system_lengths.at(0) / 2.;
  }
  try {
    system.get_property("y_mid", &y_mid, index);
  } catch (const Properties::MissingPropertyException*) {
    y_mid = system_lengths.at(1) / 2.;
  }
  try {
    system.get_property("z_mid", &z_mid, index);
  } catch (const Properties::MissingPropertyException*) {
    z_mid = system_lengths.at(2) / 2.;
  }
  // Periodic boundary conditions are ignored.
}

// _____________________________________________________________________________
void BoundarySurfaceSphere::shifted_coordinates(
		std::vector<double>& cartesian_coordinates,
		std::vector<double>* shifted) const {
  shifted->at(0) = cartesian_coordinates.at(0) - x_mid;
  shifted->at(1) = cartesian_coordinates.at(1) - y_mid;
  shifted->at(2) = cartesian_coordinates.at(2) - z_mid;
}

// _____________________________________________________________________________
bool BoundarySurfaceSphere::is_within_boundary(std::vector<double>& position)
	const {
  // shift coordinates to sphere centre
  std::vector<double> pos(3);
  shifted_coordinates(position, &pos);
  // calculate distance to sphere centre
  double distance_centre_squared = 0.;
  for (double dist : pos) {
    distance_centre_squared += dist * dist;
  }
  // decide whether point is outside sphere
  if (distance_centre_squared < radius * radius) {
    return false;
  } else {
    return true;
  }
}

// _____________________________________________________________________________
double BoundarySurfaceSphere::distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const {
  if (direction >= 3) {
    throw(std::invalid_argument("BoundarySurfaceSphere::distance_directed: "
	"invalid direction argument " + std::to_string(direction)));
  }
  if (!is_within_boundary(position)) { return 0.; }
  // shift coordinates to sphere centre
  std::vector<double> pos(3);
  shifted_coordinates(position, &pos);
  // define index_1 and index_2 as the two axes perpendicular to the given
  // direction
  size_t index_1 = (direction + 1) % 3;
  size_t index_2 = (direction + 2) % 3;
  // calculate distance to given direction's axis
  double distance_axis_squared = pos.at(index_1) * pos.at(index_1) +
  			 	 pos.at(index_2) * pos.at(index_2);
  // return -1 if there's no boundary in the given direction
  if (distance_axis_squared > radius * radius) {
    return -1.;
  }
  if ((pos.at(direction) > 0. && forward) ||
      (pos.at(direction) < 0. && !forward)) {
    return -1.;
  }
  // calculate distance to boundary
  double boundary_position = std::sqrt(radius * radius - distance_axis_squared);
  return std::fabs(pos.at(direction)) - boundary_position;
}

// _____________________________________________________________________________
void BoundarySurfaceSphere::discretize_surface(std::vector<std::vector<double>>* 
		positions, double resolution) const {
  positions->clear();
  // TODO
  std::cerr << "BoundarySurfaceSphere::discretize_surface: ";
  std::cerr << "Incomplete implementation. Check out source code.\n";
}

// _____________________________________________________________________________
void BoundarySurfaceSphere::discretize_surface(
    std::vector<std::vector<std::vector<double>>>* distribution,
    double resolution) const {
  distribution->clear();
  std::vector<std::vector<double>> positions;
  discretize_surface(&positions, resolution);
  std::vector<std::vector<double>> point(2);
  for (std::vector<double>& position : positions) {
    point.at(0) = position;
    // shift coordinates to sphere centre to obtain surface-orthogonal vector
    shifted_coordinates(position, &point.at(1));
    // Rescale surface normal such that its norm corresponds to the associated
    // surface area.
    // TODO
    distribution->push_back(point);
  }
}

// _____________________________________________________________________________
double BoundarySurfaceSphere::distance_minimal(std::vector<double>& position)
	const {
  if (!is_within_boundary(position)) {
    return 0.;
  }
  // shift coordinates to sphere centre
  std::vector<double> pos(3);
  shifted_coordinates(position, &pos);
  // calculate minimal boundary distance
  double distance_centre_squared = 0.;
  for (double dist : pos) {
    distance_centre_squared += dist * dist;
  }
  return std::sqrt(distance_centre_squared) - radius;
}
// _____________________________________________________________________________
void BoundarySurfaceSphere::minimal_distances(DataFrame<3, double>* result,
		double resolution, double upper_limit) {
  // Initialize
  std::vector<size_t> grid_counts = result->size_dim();
  std::vector<double> bin_sizes = bins(grid_counts);
  // Iterate over all grid points
  std::vector<double> position;
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      for (size_t k = 0; k != grid_counts.at(2); ++k) {
	position = {i * bin_sizes.at(0), j * bin_sizes.at(1),
		  k * bin_sizes.at(2)};
	result->at(i, j, k) = distance_minimal(position);
      }
    }
  }
}
