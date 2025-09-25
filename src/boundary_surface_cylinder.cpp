// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface_cylinder.cpp
 *  \brief Source file for the BoundarySurfaceCylinder class.
 *
 *  The file contains the definitions of the BoundarySurfaceCylinder class.
 */
#include "boundary_surface_cylinder.hpp"
#include "system.hpp"
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>

// _____________________________________________________________________________
BoundarySurfaceCylinder::BoundarySurfaceCylinder(System<3>& system, size_t index) :
	BoundarySurface<3>(system, index) {
  extract_special_properties(system, index);
}

// _____________________________________________________________________________
void BoundarySurfaceCylinder::extract_special_properties(System<3>& system, 
		size_t index) {
  system.get_property("height", &half_height, index);
  half_height *= 0.5;
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
void BoundarySurfaceCylinder::cylinder_coordinates(
		std::vector<double>& cartesian_coordinates,
		std::vector<double>* cylinder_coordinates) const {
  // shift coordinate system's origin to system center
  double x = cartesian_coordinates.at(0) - x_mid;
  double y = cartesian_coordinates.at(1) - y_mid;
  double z = cartesian_coordinates.at(2) - z_mid;
  // calculate cylinder coordinates
  // radius
  cylinder_coordinates->at(0) = sqrt(x * x + y * y);
  // angle (note that std::atan returns values from [-pi/2, pi/2])
  if (x >= 0. && y >= 0.) {
    cylinder_coordinates->at(1) = atan(y / x);
  } else if (x >= 0. && y < 0.) {
    cylinder_coordinates->at(1) = 2. * M_PI + atan(y / x);
  } else {
    cylinder_coordinates->at(1) = M_PI + atan(y / x);
  }
  // height
  cylinder_coordinates->at(2) = z;
}

// _____________________________________________________________________________
bool BoundarySurfaceCylinder::is_within_boundary(std::vector<double>& position)
	const {
  // convert to cylinder coordinates
  std::vector<double> pos(3);
  cylinder_coordinates(position, &pos);
  // check whether position is in- or outside the cylinder
  if (pos.at(0) < radius && std::fabs(pos.at(2)) < half_height) {
    return false;
  } else {
    return true;
  }
}

// _____________________________________________________________________________
double BoundarySurfaceCylinder::distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const {
  if (!is_within_boundary(position)) { return 0.; }
  // convert to cylinder coordinates
  std::vector<double> pos(3);
  cylinder_coordinates(position, &pos);
  // determine nearest boundary distance
  switch (direction) {
    case 2: { // distance in z-direction
      if (pos.at(0) < radius) { // position vertically above or below cylinder
	if (pos.at(2) >= half_height && !forward) { // above cylinder
	  return pos.at(2) - half_height;
        } else if (pos.at(2) <= - half_height && forward) { // below cylinder
	  return -pos.at(2) - half_height;
	} else { // specified direction points away from cylinder
	  return -1.;
	}
      } else { // not aligned with cylinder in z-direction
	return -1.;
      }
    }
    case 1: { // distance in y-direction
      double x = position.at(0) - x_mid;
      double y = position.at(1) - y_mid;
      if (std::fabs(pos.at(2)) < half_height && std::fabs(x) < radius) {
	if (y < 0. && forward) {
	  return radius * sin(pos.at(1)) - y;
        } else if (y > 0. && !forward) {
	  return y - radius * sin(pos.at(1));
	} else { // specified direction points away from cylinder
	  return -1.;
	}
      } else { // not aligned with cylinder in y-direction
        return -1.;
      }
    }
    case 0 : { // distance in x-direction
      double x = position.at(0) - x_mid;
      double y = position.at(1) - y_mid;
      if (std::fabs(pos.at(2)) < half_height && std::fabs(y) < radius) {
	if (x < 0. && forward) {
	  return radius * cos(pos.at(1)) - x;
        } else if (x > 0. && !forward) {
	  return x - radius * cos(pos.at(1));
	} else { // specified direction points away from cylinder
	  return -1.;
	}
      } else { // not aligned with cylinder in x-direction
        return -1.;
      }
    }
    default:
      throw(std::invalid_argument("BoundarySurfaceCylinder::distance_directed: "
	"invalid direction argument " + std::to_string(direction)));
  }
}

// _____________________________________________________________________________
void BoundarySurfaceCylinder::discretize_surface(std::vector<std::vector<double>>* 
		positions, double resolution) const {
  positions->clear();
  // discrete points distributed equidistantly in x and y direction
//  double dx = system_lengths.at(0) / static_cast<double>(resolution);
//  double dy = system_lengths.at(1) / static_cast<double>(resolution);
  double dx = resolution;
  double dy = resolution;
  double dz = resolution;
  double dphi = resolution / radius;
  // grid dimensions
  size_t grid_counts_x = static_cast<size_t>(system_lengths.at(0) / dx);
  size_t grid_counts_y = static_cast<size_t>(system_lengths.at(1) / dy);
  // surface positions
  std::vector<double> position(3);
  double x_squared, y_squared;
  double z_top = z_mid + half_height;
  double z_bottom = z_mid - half_height;
  // cylinder caps (could be implemented faster by iterating only over relevant
  // region, if necessary)
  for (size_t i = 0; i < grid_counts_x; i++) {
    position.at(0) = dx * static_cast<double>(i);
    x_squared = pow(x_mid - position.at(0), 2);
    for (size_t j = 0; j < grid_counts_y; j++) {
      position.at(1) = dy * static_cast<double>(j);
      y_squared = pow(y_mid - position.at(1), 2);
      if (x_squared + y_squared < radius * radius) {
	position.at(2) = z_top;
	positions->push_back(position);
	position.at(2) = z_bottom;
	positions->push_back(position);
      }
    }
  }
  // cylinder mantles
  double phi = 0.;
  while (phi < 2. * M_PI) {
    position.at(0) = x_mid + radius * cos(phi);
    position.at(1) = y_mid + radius * sin(phi);
    double z = z_bottom;
    while (z < z_top) {
      position.at(2) = z;
      positions->push_back(position);
      z += dz;
    }
    phi += dphi;
  }
}

// _____________________________________________________________________________
void BoundarySurfaceCylinder::discretize_surface(
    std::vector<std::vector<std::vector<double>>>* distribution,
    double resolution) const {
  distribution->clear();
  std::vector<std::vector<double>> positions;
  discretize_surface(&positions, resolution);
  double area = resolution * resolution;
  std::vector<std::vector<double>> point(2);
  for (std::vector<double>& position : positions) {
    point.at(0) = position;
    point.at(1) = surface_normal(position);
    // Rescale surface normal such that its norm corresponds to the associated
    // surface area. For that, just multiply with the area dx*dy resp.
    // radius * dphi * dz
    for (size_t i = 0; i < 3; i++) {
      point.at(1).at(i) *= area;
    }
    distribution->push_back(point);
  }
}

// _____________________________________________________________________________
std::vector<double> BoundarySurfaceCylinder::surface_normal(
	std::vector<double>& position) const {
  // calculate surface normal vector (normalised to 1.)
  std::vector<double> normal(3);
  if (position.at(2) == z_mid + half_height) { // position on upper cylinder cap
    normal = {0., 0., 1.};
  } else if (position.at(2) == z_mid - half_height) { // lower cylinder cap
    normal = {0., 0., -1.};
  } else { // position on cylinder mantle
    double x = position.at(0) - x_mid;
    double y = position.at(1) - y_mid;
    normal = {x / radius, y / radius, 0.};
    // check that position is actually on surface
    double dev = radius * radius - x * x - y * y;
    if (dev > radius * 1.e-6) { // workaround (TODO)
      std::cerr << "BoundarySurfaceCylinder::surface_normal: position ";
      for (double pos : position) {
	std::cerr << pos << " ";
      }
      std::cerr << "not on surface (deviation " << dev << ")\n";
      exit(1);
    }
  }
  return normal;
}

// _____________________________________________________________________________
double BoundarySurfaceCylinder::distance_minimal(std::vector<double>& position)
	const {
  if (!is_within_boundary(position)) {
    return 0.;
  }
  // convert to cylinder coordinates
  std::vector<double> pos(3);
  cylinder_coordinates(position, &pos);
  // calculate minimal boundary distance
  if (pos.at(0) < radius) { // position above or below cylinder
    return std::fabs(pos.at(2)) - half_height;
  } else if (std::fabs(pos.at(2)) < half_height) { // on same height as cylinder
    return pos.at(0) - radius;
  } else { // somewhere in the edges
    double distance_squared = pow(std::fabs(pos.at(2)) - half_height, 2) +
	       		      pow(pos.at(0) - radius, 2);
    return sqrt(distance_squared);
  }
}
// _____________________________________________________________________________
void BoundarySurfaceCylinder::minimal_distances(DataFrame<3, double>* result,
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
