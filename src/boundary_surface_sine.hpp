// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARY_SURFACE_SINE_HPP_
#define SRC_BOUNDARY_SURFACE_SINE_HPP_
/** \file boundary_surface_sine.hpp
 * \brief Header file for the BoundarySurfaceSine class.
 *
 *  The file contains the declarations of the BoundarySurfaceSine class.
 */
// Includes
#include "boundary_surface.hpp"
#include "system.hpp"
#include <vector>
#include "data_frame.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"

/** \brief Class for a sinusoidal boundary structure in a 3d system.
 *
 * The surface lies either at z=0 or at the opposite side of the system and is
 * invariant in y-direction, i.e. the sinusoidal structure is in the x-z-plane.
 * At both ends of the x-axis, the surface height is zero. The number of maxima
 * inbetween those points as well as their amplitude are chosen by the user.
 * Periodic boundary conditions in x and y direction are assumed.
 *
 * Required parameters (given to constructor via System object):
 * \param amplitude
 * \param maxima_count
 * \param side_inversed
 * ... plus the standard parameters listed in the documentation of the base
 * class BoundarySurface.
 */
class BoundarySurfaceSine : public BoundarySurface<3> {
  public:
    /** \brief Constructor
     */
    BoundarySurfaceSine(System<3>& system, size_t index);
    /** \brief Is this position within the boundaries?
     */
    bool is_within_boundary(std::vector<double>& position) const override;
    /** \brief Return distance to the boundary in the given direction.
     */
    double distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const override;
  private:
    /** \brief amplitude of the surface structure (maximal minus minimal height)
     */
    double amplitude;
    /** \brief wave number of the surface structure
     * 	(2 Pi / distance between nearest minima)
     */
    double wave_vector;
    /** \brief effective height of the system in z direction
     *
     * This differs from the system_lengths by one bin size.
     */
    double system_height;
    /** \brief false if surface is at z = 0 (default); true if it is opposite to
     * that (i.e. at z = system_height)
     */
    bool side_inversed;
   /** \brief Return a set of points (and normal vectors) on the surface
    *
    * The points are distributed equidistantly in x and y direction.
    *
    * \param resolution: bin size in x and y direction
    */
   void discretize_surface(
	    std::vector<std::vector<std::vector<double>>>* distribution,
	    double resolution) const override;
   void discretize_surface(
		std::vector<std::vector<double>>* points, double resolution)
	   	const override;
   /** \brief Return boundaries' z position at given x and y
    *
    * \param position: specifies x and y position (may also contain a z
    * 		coordinate, but then that is ignored)
    */
   double surface_height(std::vector<double>& position) const;
   /** \brief Return surface normal vector at given x and y
    *
    * The returned vector is orthogonal to the surface, points away from the
    * wall into the bulk. Its norm is the ratio between actual surface area
    * and the projected area bin_size.at(0)*bin_size.at(1) associated with
    * the grid point (x,y).
    *
    * \param position: specifies x and y position (may also contain a z
    * 		coordinate, but then that is ignored)
    */
   std::vector<double> surface_normal(std::vector<double>& position) const;
   /** \brief Initialize the properties specific to this subclass
    */
   void extract_special_properties(System<3>& system, size_t index);
};
#endif // SRC_BOUNDARY_SURFACE_SINE_HPP_
