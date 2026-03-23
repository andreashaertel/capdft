// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface_sine.cpp
 *  \brief Source file for the BoundarySurfaceSine class.
 *
 *  The file contains the definitions of the BoundarySurfaceSine class.
 */
#include "boundary_surface_sine.hpp"
#include "system.hpp"
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>

// _____________________________________________________________________________
BoundarySurfaceSine::BoundarySurfaceSine(System<3>& system, size_t index) :
	BoundarySurface<3>(system, index) {
  extract_special_properties(system, index);
}

// _____________________________________________________________________________
void BoundarySurfaceSine::extract_special_properties(System<3>& system, 
		size_t index) {
  // waveform parameters
  system.get_property("side_inversed", &side_inversed, index);
  system.get_property("amplitude", &amplitude, index);
  size_t maxima_count;
  system.get_property("maxima_count", &maxima_count, index);
  wave_vector = 2. * M_PI * static_cast<double>(maxima_count) / 
	  system_lengths.at(0);
  // general system parameters
  system_height = system.system_lengths.at(2) - system.bin_sizes.at(2);
//  periodic_boundaries = {true, true, false};
  for (size_t dir : {0, 1, 2}) {
    if (system.contains_property(system.indexed("PBC", dir))) {
      std::cout << "BoundarySurfaceSine: warning: periodic boundaries specified "
      	"explicitly in properties even though they are fixed to be {true, true, "
  	"false}.\n";
      break;
    }
  }
}

// _____________________________________________________________________________
double BoundarySurfaceSine::surface_height(std::vector<double>& position)
	const {
  double height = 0.5 * amplitude * (1. - cos(position.at(0) * wave_vector));
  if (side_inversed) {
    height = system_height - height;
  }
  return height;
}

// _____________________________________________________________________________
bool BoundarySurfaceSine::is_within_boundary(std::vector<double>& position)
	const {
  if (!side_inversed) {
    if (position.at(2) > surface_height(position)) {
      return true;
    } else {
      return false;
    }
  } else {
    if (position.at(2) < surface_height(position)) {
      return true;
    } else {
      return false;
    }
  }
}

// _____________________________________________________________________________
double BoundarySurfaceSine::distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const {
  if (!is_within_boundary(position)) { return 0.; }
  switch (direction) {
    case 2:
      if (forward) {
        if (side_inversed) { return surface_height(position) - position.at(2); }
        else { return -1.; }
      } else {
        if (side_inversed) { return -1.; }
        else { return position.at(2) - surface_height(position); }
      }
    case 1:
      return -1.;
    case 0 : {
      // calculate distance of given point from baseline z=0 or z=system_height
      double height;
      if (side_inversed) { height = system_height - position.at(2); }
      else { height = position.at(2); }
      if (height > amplitude) { return -1.; }
      // calculate cosine phase corresponding to the given point
      double phase_point = // in the interval [-M_PI, M_PI]
		std::fmod(position.at(0) * wave_vector + M_PI, 2. * M_PI) - M_PI;
      // calculate position of nearest surface in forward / backward direction
      // and the corresponding phase difference
      if (forward) {
        double phase_surface = acos(1. - 2. * height / amplitude); // [0,M_PI]
	return (phase_surface - phase_point) / wave_vector;
      } else {
        double phase_surface = - acos(1. - 2. * height / amplitude); // [-M_PI,0]
	return (phase_point - phase_surface) / wave_vector;
      }
    }
    default:
      throw(std::invalid_argument("BoundarySurfaceSine::distance_directed: "
	"invalid direction argument " + std::to_string(direction)));
  }
}

// _____________________________________________________________________________
void BoundarySurfaceSine::discretize_surface(std::vector<std::vector<double>>* 
		positions, double resolution) const {
  positions->clear();
  /* Discrete points distributed equidistantly in x and y direction:
   * We want a maximum distance between any two neighboring points given by
   * 'resolution'. Therefore, we choose dx such that in the steepest surface
   * region, the distance ds = sqrt(dx² + dz²) between two surface points is
   * equal to the given resolution.
   */
  double dx = resolution / std::sqrt(1 + pow(amplitude * wave_vector / 2., 2));
  double dy = resolution;
  size_t grid_count_x = static_cast<size_t>(system_lengths.at(0) / dx);
  size_t grid_count_y = static_cast<size_t>(system_lengths.at(1) / dy);
  std::vector<double> position(3);
  for (size_t i = 0; i < grid_count_x; i++) {
    position.at(0) = dx * static_cast<double>(i);
    for (size_t j = 0; j < grid_count_y; j++) {
      position.at(1) = dy * static_cast<double>(j);
      position.at(2) = surface_height(position);
      positions->push_back(position);
    }
  }
}

// _____________________________________________________________________________
void BoundarySurfaceSine::discretize_surface(
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
    // surface area. For that, just multiply with the projected area dx*dy; the
    // factor for conversion to the actual (curved) area is already included in
    // the return value of surface_normal.
    for (size_t i = 0; i < 3; i++) {
      point.at(1).at(i) *=
	      area / std::sqrt(1 + pow(amplitude * wave_vector / 2., 2));
    }
    distribution->push_back(point);
  }
}

// _____________________________________________________________________________
std::vector<double> BoundarySurfaceSine::surface_normal(
	std::vector<double>& position) const {
  double x = position.at(0);
  // calculate surface normal vector (not yet normalised)
  double derivative_x = 0.5 * amplitude * wave_vector * sin(x * wave_vector);
  std::vector<double> normal(3);
  if (!side_inversed) {
    normal = {- derivative_x, 0., 1.};
  } else {
    normal = {- derivative_x, 0., -1.};
  }
  return normal;
}

