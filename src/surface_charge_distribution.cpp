// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file surface_charge_distribution.cpp
 * \brief Source file for the surface charge calculation
 */
#include "system.hpp"
#include <vector>
#include "surface_charge_distribution.hpp"
#include <utility>
#include "boundary_surface.hpp"
#include <cmath>
#include "system.hpp"
#include "stl_algorithms.hpp"

// _____________________________________________________________________________
template <>
double charge_density(BoundarySurface<3>& boundary, System<3>& system,
	DataFrame<3, double>& potential,
	std::vector<double>& position_boundary, std::vector<double>& normal) {
  /** The local charge density at any surface point can be determined from the derivative of the
   * electrostatic potential in direction orthogonal to the surface. To
   * calculate this normal derivative, we consider a point which is at a small (~ bin size)
   * distance from the surface point (in direction of the bulk, therefore called 'position_bulk').
   * The potential value at this 'bulk' point is determined via interpolation. From the two
   * potential values at the surface and the 'bulk' point, we calculate the derivative.
   */
  double prefactor = - 1. / (4. * M_PI * system.bjerrum);
  double boundary_distance = system.bin_sizes.at(0); // interpolation distances ~ bin size
  double boundary_potential = boundary.get_boundary_value();
  // define 'bulk' point
  std::vector<double> position_bulk(3);
  for (size_t i = 0; i < 3; i++) {
    position_bulk.at(i) = position_boundary.at(i) +
		normal.at(i) * boundary_distance;
  }
  // calculate normal derivative and from that the charge density
//  double normal_derivative = (system.interpolate(potential, position_bulk) 
//		    	 - boundary_potential) / boundary_distance;
  double bulk_potential = system.interpolate(potential, position_bulk);
  double normal_derivative = std::log(bulk_potential / boundary_potential) * 
			     boundary_potential / boundary_distance;
  return prefactor * normal_derivative;
}
// _____________________________________________________________________________
template <size_t dim>
void charge_distribution(BoundarySurface<dim>& boundary, System<dim>& system,
	DataFrame<dim, double>& potential,
	std::vector<std::pair<std::vector<double>, double>>* distribution,
   	double* total_charge, double resolution) {
  std::cerr << "charge_distribution: this is only implemented in 3d\n";
  exit(1);
}
// _____________________________________________________________________________
template <>
void charge_distribution(BoundarySurface<3>& boundary, System<3>& system,
	DataFrame<3, double>& potential,
	std::vector<std::pair<std::vector<double>, double>>* distribution,
   	double* total_charge, double resolution) {
  // initialize
  std::cout << "init\n";
  distribution->clear();
  *total_charge = 0.;
  std::vector<double> bins = system.bin_sizes;
  std::vector<std::vector<std::vector<double>>> surface_points;
  boundary.discretize_surface(&surface_points, resolution);
  std::vector<double> position_boundary(3);
  std::vector<double> normal(3);
  std::pair<std::vector<double>, double> result;
  std::cout << "calc\n";
  double norm;
  double charge;
  double area = 0.;
  for (std::vector<std::vector<double>>& point : surface_points) {
    // define surface point
    position_boundary = point.at(0);
    normal = point.at(1);
    norm = stl_algorithm::vector_norm(normal);
    for (double& value : normal) {
      value /= norm;
    }
    // calculate charge density
    charge = charge_density(boundary, system, potential, position_boundary, normal);
    // save results
    result.first = position_boundary;
    result.second = charge;
    distribution->push_back(result);
    /** The total charge is calculated by numerical integration over all surface points.
     * The integration measure (area attributed to the surface point) is given by the
     * norm of the normal vectors returned by BoundarySurface::discretize_surface. 
     */
    *total_charge += charge * norm;
    area += norm;
  }
  std::cout << "total area: " << area << std::endl;
}
// alternative method:
//// _____________________________________________________________________________
//template <>
//void charge_distribution_poisson(BoundarySurface<3>& boundary, System<3>& system,
//	DataFrame<3, double>& potential,
//	std::vector<std::pair<std::vector<double>, double>>* distribution,
//   	double* total_charge, double resolution) {
//  // initialize
//  std::cout << "init\n";
//  distribution->clear();
//  *total_charge = 0.;
//  double boundary_value = boundary.get_boundary_value();
//  std::vector<double> bins = system.bin_sizes;
//  double voxel_volume = bins.at(0) * bins.at(1) * bins.at(2);
//  double prefactor = - voxel_volume / (4. * M_PI * system.bjerrum);
//  std::vector<std::vector<std::vector<double>>> surface_points;
//  boundary.discretize_surface(&surface_points, resolution);
//  std::vector<double> position_boundary(3);
//  std::vector<double> normal(3);
//  std::pair<std::vector<double>, double> result;
//  // calculate Laplacian of potential at each surface point
//  std::cout << "calc\n";
//  double norm;
//  double charge;
//  for (std::vector<std::vector<double>>& point : surface_points) {
//    position_boundary = point.at(0);
//    normal = point.at(1);
//    norm = stl_algorithm::vector_norm(normal);
//    charge = prefactor * system.laplace(
//	      potential, position_boundary, boundary_value);
//    result.first = position_boundary;
//    result.second = charge / norm;
//    distribution->push_back(result);
//    *total_charge += charge;
//  }
//}
