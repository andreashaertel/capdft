// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARY_SURFACE_CYLINDER_HPP_
#define SRC_BOUNDARY_SURFACE_CYLINDER_HPP_
/** \file boundary_surface_cylinder.hpp
 * \brief Header file for the BoundarySurfaceCylinder class.
 *
 *  The file contains the declarations of the BoundarySurfaceCylinder class.
 */
// Includes
#include "boundary_surface.hpp"
#include "system.hpp"
#include <vector>
#include "data_frame.hpp"
#include "../../parameter_handler/src/parameter_handler.hpp"

/** \brief Class for a cylinder-shaped boundary surface.
 *
 * The cylinder is oriented in z-direction. Per default, it is located in the
 * middle of the system, but optional parameters allow to set it at any other
 * position. The radius and height of the cylinder are specified by the user.
 *
 * Required parameters (given to constructor via System object):
 * \param radius
 * \param height
 * Optional parameters (given to constructor via System object):
 * \param x_mid, y_mid, z_mid
 * ... plus the standard parameters listed in the documentation of the base
 * class BoundarySurface.
 */
class BoundarySurfaceCylinder : public BoundarySurface<3> {
  public:
    /** \brief Constructor
     */
    BoundarySurfaceCylinder(System<3>& system, size_t index);
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
    /** \brief Convert cartesian coordinates to cylinder coordinates
     */
    void cylinder_coordinates(std::vector<double>& cartesian,
		    std::vector<double>* cylinder) const;
  private:
    /** \brief Half of the cylinder height
     */
    double half_height;
    /** \brief Cylinder radius
     */
    double radius;
    /** \brief Position of cylinder centre
     */
    double x_mid, y_mid, z_mid;
   /** \brief Return a set of points (and normal vectors) on the surface
    *
    * The points that lie on the cylinder caps are distributed on a square grid
    * in the x-y-plane. Those that lie on the cylinder mantle are distributed
    * on a square grid wrapped around the mantle.
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
    * wall into the bulk. Its norm is the ratio between actual surface area
    * and the projected area bin_size.at(0)*bin_size.at(1) associated with
    * the grid point (x,y).
    */
   std::vector<double> surface_normal(std::vector<double>& position) const;
   /** \brief Initialize the properties specific to this subclass
    */
   void extract_special_properties(System<3>& system, size_t index);
};

/** \example ../examples/boundary_surface_cylinder/src/main.cpp
 * This is an example for the BSCylinder.
 * A cylinder with fixed potential value...
 */
#endif // SRC_BOUNDARY_SURFACE_CYLINDER_HPP_
