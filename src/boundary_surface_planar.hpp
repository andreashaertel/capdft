// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARY_SURFACE_PLANAR_HPP_
#define SRC_BOUNDARY_SURFACE_PLANAR_HPP_
/** \file boundary_surface_planar.hpp
 * \brief Header file for the BoundarySurfacePlanar class.
 *
 *  The file contains the declarations of the BoundarySurfacePlanar class.
 */
// Includes
#include "boundary_surface.hpp"
#include "system.hpp"
#include <vector>
#include "data_frame.hpp"

/** \brief Class template for a planar boundary surface.
 *
 * Can be defined for any dimension; only the 3d version is implemented so far.
 * The surface can consist of any combination of any of the (in 3d case) six system walls.
 *
 * Required parameters (given to constructor via System object):
 * \param 0_sides
 * \param 1_sides
 * \param 2_sides
 * ... plus the standard parameters listed in the documentation of the base
 * class BoundarySurface.
 */
template <size_t dim>
class BoundarySurfacePlanar : public BoundarySurface<dim> {
  public:
    /** \brief Construct structure from given System object.
     *
     * Special parameters required for this surface object: For each coordinate
     * axis, one parameter 'sides' encodes which of the two walls along that axis
     * are to be included in the surface object.
     */
    BoundarySurfacePlanar(System<dim>& system, size_t index);
    /** \brief Is this position within the boundaries?
     */
    bool is_within_boundary(std::vector<double>& position) const override;
    /** \brief Return distance to the boundary in the given direction.
     */
    double distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const override;
    /** \brief Return distance to the nearest boundary point.
     *
     * Returns zero if given point is outside boundaries. Returns upper_limit if
     * the distance to the boundary is larger than that value.
     *
     * \param upper_limit: use this if only distances below a certain limit are of interest
     */
    double distance_minimal(std::vector<double>& position, double upper_limit) const;
    /** \brief Return the minimal boundary distances of each point in a
     * discretized grid.
     *
     * \param result: pointer to return value; needs to be initialized with
     * 		desired grid dimensions
     * \param resolution: spatial resolution of surface (see distance_minimal)
     * \param upper_limit: upper limit on distance (see distance_minimal) - By
     * 		default, there is no limit except for the finite system size.
     */
    void minimal_distances(DataFrame<dim, double>* result,
		    double resolution, double upper_limit = -1.) override;
  private:
   /** \brief Specifies which of the six system walls are part of the surface
    */
   std::vector<std::vector<size_t>> sides;
   /** \brief Effective height of the system between opposite surfaces
    *
    * This differs from the system_lengths by one bin size.
    */
   std::vector<double> system_heights;
   /** \brief Calculate surface normal vector at any point on the surface
    */
   std::vector<double> surface_normal(std::vector<double>& position) const;
   /** \brief Return a set of points (and normal vectors) on the surface
    *
    * The points are distributed equidistantly on a square grid.
    * \param resolution: bin size in all directions
    */
   void discretize_surface(
	    std::vector<std::vector<std::vector<double>>>* distribution,
	    double resolution) const override;
   void discretize_surface(
		std::vector<std::vector<double>>* points, double resolution)
	   	const override;
   /** \brief Initialize the properties specific to this subclass
    */
   void extract_special_properties(System<dim>& system, size_t index);
};
#endif // SRC_BOUNDARY_SURFACE_PLANAR_HPP_
