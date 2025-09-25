// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file boundary_surface.cpp
 *  \brief Source file for the BoundarySurface class.
 *
 *  The file contains the definitions of the BoundarySurface class.
 */
#include <vector>
#include <cmath>
#include <algorithm>
#include "boundary_surface.hpp"
#include "system.hpp"
#include "properties.hpp"
#include "data_frame.hpp"
#include <string>

// Template class instantiation
template class BoundarySurface<3>;

// _____________________________________________________________________________
template <size_t dim>
BoundarySurface<dim>::BoundarySurface(System<dim>& system, size_t index) {
  extract_standard_properties(system, index);
}

// _____________________________________________________________________________
template <size_t dim>
void BoundarySurface<dim>::set_boundary_value(double value, Type boundary_type) {
  boundary_value = value;
  type = boundary_type;
}

// _____________________________________________________________________________
template <size_t dim>
double BoundarySurface<dim>::get_boundary_value() const {
  if (type == None) {
    std::cerr << "BoundarySurface::get_boundary_value: no boundary values\n";
    exit(1);
  }
  return boundary_value;
}

// _____________________________________________________________________________
template <size_t dim>
double BoundarySurface<dim>::distance_minimal(std::vector<double>& position, 
		double resolution, double upper_limit) const {
  if (!is_within_boundary(position)) {
    return 0.;
  }
  std::vector<std::vector<double>> surface_points;
  discretize_surface(&surface_points, resolution);
  double minimum_squared = upper_limit * upper_limit;
  for (std::vector<double>& point : surface_points) {
    double distance_squared = 0.;
    for (size_t dir = 0; dir < dim; dir++) {
      distance_squared += pow(point.at(dir) - position.at(dir), 2);
    }
    if (distance_squared < minimum_squared) {
      minimum_squared = distance_squared;
    }
  }
  return sqrt(minimum_squared);
}
// _____________________________________________________________________________
template <size_t dim>
double BoundarySurface<dim>::distance_minimal(std::vector<double>& position, 
	std::vector<std::vector<double>>& surface_points, double upper_limit) const {
  if (!is_within_boundary(position)) {
    return 0.;
  }
  double minimum_squared = upper_limit * upper_limit;
  for (std::vector<double>& point : surface_points) {
    double distance_squared = 0.;
    for (size_t dir = 0; dir < dim; dir++) {
      distance_squared += pow(point.at(dir) - position.at(dir), 2);
    }
    if (distance_squared < minimum_squared) {
      minimum_squared = distance_squared;
    }
  }
  return sqrt(minimum_squared);
}
// _____________________________________________________________________________
template <size_t dim>
void BoundarySurface<dim>::extract_standard_properties(System<dim>& system,
		size_t index) {
  system_lengths = system.system_lengths;
  try { // optional parameter: boundary value
    system.get_property("boundary_value", &boundary_value, index);
    type = PotentialES;
  } catch (const Properties::MissingPropertyException*) {
    type = None;
  }
  try { // optional parameter: hard walls?
    system.get_property("hard_walls", &hard_walls, index);
  } catch (const Properties::MissingPropertyException*) {
    hard_walls = true;
  }
} 
 
// _____________________________________________________________________________
template <>
void BoundarySurface<3>::discretize_volume(DataFrame<3, bool>* grid) {
  std::vector<size_t> grid_counts = grid->size_dim();
  std::vector<double> bin_sizes = bins(grid_counts);
  std::vector<double> position(3);
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    position.at(0) = bin_sizes.at(0) * static_cast<double>(i);
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      position.at(1) = bin_sizes.at(1) * static_cast<double>(j);
      for (size_t k = 0; k != grid_counts.at(2); ++k) {
        position.at(2) = bin_sizes.at(2) * static_cast<double>(k);
	grid->at(i, j, k) = is_within_boundary(position);
      }
    }
  }
}
// _____________________________________________________________________________
template <>
void BoundarySurface<3>::minimal_distances(DataFrame<3, double>* result,
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
  // Construct a dataframe containing which points in the grid are within the 
  // boundary, such that this has not to be calculated repeatedly in the 
  // following.
  DataFrame<3, bool> grid(grid_counts);
  discretize_volume(&grid);
  std::vector<std::vector<double>> surface_points;
  discretize_surface(&surface_points, resolution);
  // Iterate over all grid points
  std::vector<double> position;
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      for (size_t k = 0; k != grid_counts.at(2); ++k) {
        if (!grid.at(i, j, k)) {
	  result->at(i, j, k) = 0.;
        } else {
	  position = {i * bin_sizes.at(0), j * bin_sizes.at(1),
		  k * bin_sizes.at(2)};
	  result->at(i, j, k) = distance_minimal(position, surface_points, 
			  upper_limit);
	}
      }
    }
  }
}
// _____________________________________________________________________________
template <>
void BoundarySurface<3>::exp_external_potential_hs(System<3>& system,
		std::vector<DataFrame<3, double>>* potential) {
  // Initialize
  std::vector<size_t> grid_counts = potential->at(0).size_dim();
  potential->resize(system.species_properties.size());
  // set all to 1 for unaffected species or if hard_walls is false:
  for (size_t s = 0; s < potential->size(); s++) {
    potential->at(s).set_all_elements_to(1.);
  }
  if (!hard_walls) {
    return;
  }
  // Get species properties
  std::vector<size_t> affected_species = system.affected_species_fmt;
  std::vector<double> diameters(system.species_properties.size());
  double diameter;
  for (size_t s : affected_species) {
    system.species_properties.at(s).get_property("diameter", &diameter);
    diameters.at(s) = diameter;
  }
  // Calculate distance to boundary from each grid point in the system
  DataFrame<3, double> distances(grid_counts);
  // TODO: surface resolution arbitrarily hardcoded?
  double resolution = 0.25 * 
	  *std::min_element(diameters.begin(), diameters.end());
  // upper limit: all distances above the largest ion's radius are irrelevant
  double upper_limit = 0.5 * 
	  *std::max_element(diameters.begin(), diameters.end());
  minimal_distances(&distances, resolution, upper_limit);
  // Iterate over all points in the discretized grid, and compare the boundary
  // distance to the particle radius.
  for (size_t s : affected_species) {
    diameter = diameters.at(s);
    for (size_t i = 0; i != grid_counts.at(0); ++i) {
      for (size_t j = 0; j != grid_counts.at(1); ++j) {
        for (size_t k = 0; k != grid_counts.at(2); ++k) {
          if (distances.at(i, j, k) < diameter * 0.5) {
            potential->at(s).at(i,j,k) = 0.;
          } else {
            potential->at(s).at(i,j,k) = 1.;
          }	
	}
      } 
    } // end of (i, j, k)-loops
  } // end of species loop
}

// _____________________________________________________________________________
template <size_t dim>
bool BoundarySurface<dim>::is_electrostatic() const {
  switch (type) {
  case PotentialES: return true;
  case None: return false;
  }
  std::cerr << "BoundarySurface::is_electrostatic: unknown boundary value type ";
  std::cerr << type << std::endl;
  exit(1);
}

// _____________________________________________________________________________
template <size_t dim>
std::vector<double> BoundarySurface<dim>::bins(std::vector<size_t>&
		grid_counts) const {
  std::vector<double> result(dim);
  for (size_t i = 0; i < dim; i++) { 
    result.at(i) = system_lengths.at(i) / static_cast<double>(grid_counts.at(i));
  }
  return result;
}

