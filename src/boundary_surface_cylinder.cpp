// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface_cylinder.cpp
 *  \brief Source file for the BoundarySurfaceCylinder class.
 *
 *  The file contains the definitions of the BoundarySurfaceCylinder class.
 */
#include "boundary_surface_cylinder.hpp"
#include "system.hpp"
#include "discretization_handler.hpp"
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>

// _____________________________________________________________________________
BoundarySurfaceCylinder::BoundarySurfaceCylinder(System<Cartesian>& system, size_t index) :
	BoundarySurface<Cartesian>(system, index) {
  extract_special_properties(system, index);
}

// _____________________________________________________________________________
void BoundarySurfaceCylinder::extract_special_properties(System<Cartesian>& system, 
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
  try {
    system.get_property("caps", &caps, index);
  } catch (const Properties::MissingPropertyException*) {
    caps = "flat";
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
double BoundarySurfaceCylinder::distance_to_spherical_cap(
		std::vector<double>& cylinder_coordinates) const {
  double rho = cylinder_coordinates.at(0);
  double z = cylinder_coordinates.at(2);
  // shift the coordinate system's origin to the cap center, choosing the cap 
  // which is closest to the given position
  if (z > 0.) {
    z = cylinder_coordinates.at(2) - half_height;
  } else {
    z = cylinder_coordinates.at(2) + half_height;
 }
 // calculate distance to cap surface
 return sqrt(rho * rho + z * z) - radius;
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
    if (caps == "spherical") {
      if (distance_to_spherical_cap(pos) < 0.) {
        return false;
      }
    }
  }
  return true;
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
	double true_height = half_height;
	if (caps == "spherical") {
	  // add cap height to obtain true cylinder height
	  true_height += sqrt(radius * radius - pos.at(0) * pos.at(0));
	}
	if (pos.at(2) >= true_height && !forward) { // above cylinder
	  return pos.at(2) - true_height;
        } else if (pos.at(2) <= - true_height && forward) { // below cylinder
	  return -pos.at(2) - true_height;
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
      } else {
        if (caps == "spherical") {
	    double z_to_cap = std::fabs(pos.at(2)) - half_height;
	    double x_to_cap = position.at(0) - x_mid;
	    double r_xz_sq = x_to_cap * x_to_cap + z_to_cap * z_to_cap;
	    if (r_xz_sq < radius * radius) { // aligned with spherical cap
	      double y = position.at(1) - y_mid;
	      double y_cap = sqrt(radius * radius - r_xz_sq);
	      if (y < 0 && forward) {
                return - y - y_cap;
	      } else if (y > 0 && !forward) {
                return y - y_cap;
	      }
	    }
	} // else: not aligned with cylinder in y-direction
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
      } else {
        if (caps == "spherical") {
	    double z_to_cap = std::fabs(pos.at(2)) - half_height;
	    double y_to_cap = position.at(1) - y_mid;
	    double r_yz_sq = y_to_cap * y_to_cap + z_to_cap * z_to_cap;
	    if (r_yz_sq < radius * radius) { // aligned with spherical cap
	      double x = position.at(0) - x_mid;
	      double x_cap = sqrt(radius * radius - r_yz_sq);
	      if (x < 0 && forward) {
                return - x - x_cap;
	      } else if (x > 0 && !forward) {
                return x - x_cap;
	      }
	    }
	} // else: not aligned with cylinder in x-direction
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
  // discrete points distributed equidistantly with distance = resolution 
  double dz = resolution;
  double dphi = resolution / radius;
  // position vector
  std::vector<double> position(3);
  // Discretize cylinder mantle:
  double z_top = z_mid + half_height;
  double z_bottom = z_mid - half_height;
  double phi = 0.;
  while (phi < 2. * M_PI) {
    position.at(0) = x_mid + radius * cos(phi);
    position.at(1) = y_mid + radius * sin(phi);
    double z = z_bottom + 0.5 * dz;
    while (z < z_top) {
      position.at(2) = z;
      positions->push_back(position);
      z += dz;
    }
    phi += dphi;
  }
  // Discretize cylinder caps:
  if (caps == "flat") {
    // Use discrete points distributed equidistantly in x and y direction.
    double dx = resolution;
    double dy = resolution;
    // grid dimensions
    size_t grid_counts_x = static_cast<size_t>(system_lengths.at(0) / dx);
    size_t grid_counts_y = static_cast<size_t>(system_lengths.at(1) / dy);
    // iterate over x-y-grid (could be implemented faster by iterating only
    // over relevant region, if necessary)
    double x_squared, y_squared;
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
  }
  if (caps == "spherical") {
    // Use discrete points distributed equidistantly along each line of latitude
    // (parametrized by phi). The angular spacing dphi depends on the longitude
    // theta such that the distances R*sin(theta)*dphi=resolution=R*dtheta both
    // equal (approximately) the desired resolution. We slightly adapt the spacing
    // to obtain a symmetric distribution over the full range [0,PI/2] / [0,2PI).
    double z_top, z_bottom;
    // For longitudinal distribution: Add 0.5*resolution to have the correct
    // distance to the points on the cylinder mantle.
    double length_theta = 0.5 * M_PI * radius - 0.5 * resolution;
    size_t grid_count_theta = static_cast<size_t>(length_theta / resolution) + 1;
    double dtheta = length_theta / radius / static_cast<double>(grid_count_theta - 1);
    double theta = 0.;
    // Latitudinal distribution:
    size_t grid_count_phi;
    double rho, dphi, phi;
    for (size_t i = 0; i < grid_count_theta; i++) {
      z_top = z_mid + half_height + radius * cos(theta);
      z_bottom = z_mid - half_height - radius * cos(theta);
      // For latitudinal distribution:
      rho = radius * sin(theta);
      phi = 2. * M_PI * rho; // maximum
      grid_count_phi = static_cast<size_t>(phi / resolution);
      if (grid_count_phi == 0) grid_count_phi += 1;
      dphi = 2. * M_PI / static_cast<double>(grid_count_phi);
      phi = 0.; // minimum, starting point for iteration
      for (size_t j = 0; j < grid_count_phi; j++) {
        position.at(0) = x_mid + rho * cos(phi);
        position.at(1) = y_mid + rho * sin(phi);
        position.at(2) = z_top;
        positions->push_back(position); // upper cap
        position.at(2) = z_bottom;
        positions->push_back(position); // lower cap
        phi += dphi;
      }
      theta += dtheta;
    }
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
    // surface area. For that, just multiply with the area dx*dy or
    // radius² * sin(theta) * dphi * dtheta or radius * dphi * dz, respectively.
    // The position discretization has been chosen such that these areas are all
    // equal to area = resolution².
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
  if (position.at(2) >= z_mid + half_height) { // position on upper cylinder cap
    if (caps == "spherical") {
      normal = {(position.at(0) - x_mid) / radius,
	        (position.at(1) - y_mid) / radius, 
	        (position.at(2) - z_mid - half_height) / radius};
    } else { 
      normal = {0., 0., 1.};
    }
  } else if (position.at(2) <= z_mid - half_height) { // lower cylinder cap
    if (caps == "spherical") {
      normal = {(position.at(0) - x_mid) / radius,
	        (position.at(1) - y_mid) / radius, 
	        (position.at(2) - z_mid + half_height) / radius};
    } else { 
      normal = {0., 0., -1.};
    }
  } else { // position on cylinder mantle
    double x = position.at(0) - x_mid;
    double y = position.at(1) - y_mid;
    normal = {x / radius, y / radius, 0.};
//    // Check that position is actually on surface?
//    // But the allowed deviations depend on the used bin size...
//    double dev = radius * radius - x * x - y * y;
//    if (dev > radius * 1.e-6) {
//      std::cerr << "BoundarySurfaceCylinder::surface_normal: position ";
//      for (double pos : position) {
//	std::cerr << pos << " ";
//      }
//      std::cerr << "not on surface (deviation " << dev << ")\n";
//      exit(1);
//    }
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
  if (std::fabs(pos.at(2)) < half_height) { // on same height as cylinder
    return pos.at(0) - radius;
  } else if (caps == "spherical") {
    return distance_to_spherical_cap(pos);
  } else { // flat caps
    if (pos.at(0) < radius) { // position vertically above or below cylinder
      return std::fabs(pos.at(2)) - half_height;
    } else { // somewhere in the edges
      double distance_squared = pow(std::fabs(pos.at(2)) - half_height, 2) +
    		      pow(pos.at(0) - radius, 2);
      return sqrt(distance_squared);
    }
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
