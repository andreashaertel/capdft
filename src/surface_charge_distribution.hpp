// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_SURFACE_CHARGE_DISTRIBUTION_HPP_
#define SRC_SURFACE_CHARGE_DISTRIBUTION_HPP_
/** \file surface_charge_distribution.hpp
 * \brief Header file for surface charge calculation methods
 */
// _____________________________________________________________________________
// Includes
#include "boundary_surface.hpp"
#include "data_frame.hpp"
#include "system.hpp"
#include <vector>
// _____________________________________________________________________________
/** \brief Calculate the surface charge density at a given position from a
 * potential profile
 *
 * This function uses the potential's normal gradient to determine the surface
 * charge density. Warning: This function doesn't check whether the given
 * position actually lies on the boundary surface; this lies in the
 * responsibility of the user.
 *
 * \param boundary: object describing the boundary surface
 * \param system: object containing the system properties
 * \param potential: potential profile, on an equidistant grid that spans the
 * 	system lengths
 * \param position_boundary: position where the charge density is to be calculated
 * \param normal: surface normal vector at the given position (norm 1)
 */
template <size_t dim>
double charge_density(BoundarySurface<dim>& boundary, System<dim>& system,
	DataFrame<dim, double>& potential,
	std::vector<double>& position_boundary, std::vector<double>& normal);
/** \brief Calculate the surface charge distribution from a potential profile
 *
 * This function uses the potential's normal gradient to determine the surface
 * charge density.
 *
 * \param boundary: object describing the boundary surface
 * \param system: object containing the system properties
 * \param potential: potential profile, on an equidistant grid that spans the
 * 	system lengths
 * \param distribution: pointer to return value for the charge distribution
 * \param total_charge: pointer to return value for the total surface charge
 * \param resolution: spatial resolution of the surface discretization (as in
 * 		BoundarySurface::discretize_surface)
 * \param total_area: pointer to return value for the total surface area
 */
template <size_t dim>
void charge_distribution(BoundarySurface<dim>& boundary, System<dim>& system,
   DataFrame<dim, double>& potential,
   std::vector<std::pair<std::vector<double>, double>>* distribution,
   double* total_charge, double resolution, double* total_area = nullptr);

#endif // SRC_SURFACE_CHARGE_DISTRIBUTION_HPP_
