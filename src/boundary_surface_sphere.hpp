// SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARY_SURFACE_SPHERE_HPP_
#define SRC_BOUNDARY_SURFACE_SPHERE_HPP_
/** \file boundary_surface_sphere.hpp
 * \brief Header file for the BoundarySurfaceSphere class.
 *
 *  The file contains the declarations of the BoundarySurfaceSphere class.
 */
// Includes
#include "boundary_surface.hpp"
#include "system.hpp"
#include <vector>
#include "data_frame.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"

/** \brief Class for a spherical boundary surface.
 *
 * Per default, the sphere is located in the middle of the system,
 * but optional parameters allow to set it at any other
 * position. The radius is specified by the user.
 *
 * Required parameters (given to constructor via System object):
 * \param radius
 * Optional parameters (given to constructor via System object):
 * \param x_mid, y_mid, z_mid
 * ... plus the standard parameters listed in the documentation of the base
 * class BoundarySurface.
 */
class BoundarySurfaceSphere : public BoundarySurface<3> {
  public:
    /** \brief Constructor
     */
    BoundarySurfaceSphere(System<3>& system, size_t index);
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
     */
    double distance_minimal(std::vector<double>& position) const;
    /** \brief Return the minimal boundary distances of each point in a
     * discretized grid.
     *
     * \param result: pointer to return value; needs to be initialized with
     * 		desired grid dimensions
     * \param resolution: spatial resolution of surface (see distance_minimal)
     * \param upper_limit: upper limit on distance (see distance_minimal) - By
     * 		default, there is no limit except for the finite system size.
     */
    void minimal_distances(DataFrame<3, double>* result,
		    double resolution, double upper_limit = -1.) override;
    /** \brief Shift coordinate system such that sphere centre is at zero 
     */
    void shifted_coordinates(std::vector<double>& cartesian,
		    std::vector<double>* shifted) const;
  private:
    /** \brief Sphere radius
     */
    double radius;
    /** \brief Position of sphere centre
     */
    double x_mid, y_mid, z_mid;
    /** \brief Return a set of points (and normal vectors) on the surface
    * 
    * \param resolution: bin size of the discretization
    */
    void discretize_surface(
	    std::vector<std::vector<std::vector<double>>>* distribution,
	    double resolution) const override;
    void discretize_surface(
		std::vector<std::vector<double>>* points, double resolution)
	   	const override;
    /** \brief Return surface normal vector at given x and y
    *
    * The returned vector is orthogonal to the surface, points away from the
    * wall into the bulk.
    */
   std::vector<double> surface_normal(std::vector<double>& position) const;
   /** \brief Initialize the properties specific to this subclass
    */
   void extract_special_properties(System<3>& system, size_t index);
};
#endif // SRC_BOUNDARY_SURFACE_SPHERE_HPP_
