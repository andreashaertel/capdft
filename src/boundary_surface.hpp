// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARY_SURFACE_HPP_
#define SRC_BOUNDARY_SURFACE_HPP_
/** \file boundary_surface.hpp
 *  \brief Header file for the BoundarySurface class.
 *
 *  The file contains the declarations of the BoundarySurface class.
 */
// Includes
#include <vector>
#include "data_frame.hpp"
#include "system.hpp"
#include <string>

/** \brief Abstract class for the boundary surface of a system.
 *
 * Each BoundarySurface object represents the wall or part of the walls of a
 * system. This can be a hard wall, and it can also have an associated
 * (uniform) boundary value, e.g. a constant potential on the surface.
 * Subclasses represent different surface structures. Multiple surface objects
 * can be combined via the Boundaries class to define the entire boundary of a
 * system.
 *
 * Required parameters (given to constructor via System object):
 * \param system_lengths (member variable of System, no index required)
 * Optional parameters (given to constructor via System object):
 * \param boundary_value
 * \param hard_walls
 * ... plus the special parameters listed in the documentation of the
 * subclass of interest.
 * Unless specified otherwise, all these parameters need a common index which
 * distinguishes different boundary objects from each other (cf. constructor).
 * This indexing is only relevant for the BoundarySurface- (and -subclass-)
 * constructors, it does not have to follow any specific order and is in
 * particular independent from the ordering of these objects within a Boundaries
 * container.
 */
template <size_t dim>
class BoundarySurface {
  public:
    /** \brief General constructor
     *
     * Since the System object can contain parameters for multiple
     * BoundarySurface objects, these properties need to be indexed (indices are
     * usually put in front of parameter names, cf. Properties::indexed). The 
     * parameter 'index' specifies which of these parameters the Constructor
     * should use.
     *
     * \param system: specifies system dimensions (system lengths, maybe PBC),
     * 		boundary parameters (e.g. sine amplitude / sphere radius...),
     * 		and optional boundary value
     * \param index: specifies which set of boundary parameters in the system
     * 		object should be used
     */
    BoundarySurface(System<dim>& system, size_t index);
    /** \brief Different types of boundary values
     */
    enum Type {PotentialES, None};
    /** \brief Set boundary value (uniform over whole surface)
     */
    void set_boundary_value(double value, Type boundary_type = PotentialES);
    /** \brief Get boundary value
     */
    double get_boundary_value() const;
    /** \brief Get boundary type
     */
    Type get_type() const;
    bool is_electrostatic() const;
    /** \brief Is this position within the boundaries?
     */
    virtual bool is_within_boundary(std::vector<double>& position) const = 0;
    /** \brief Calculate for each position in a discretized grid whether it is
     * in- or outside the boundaries
     *
     * \param grid: pointer to return value; needs to be initialized with the
     * 		desired grid dimensions
     */
    void discretize_volume(DataFrame<dim, bool>* grid);
    /** \brief Discretize the surface into a list of points
     *
     * This function returns an (unordered) list of points on the surface.
     *
     * \param positions: pointer to return value
     * \param resolution: spatial resolution; exact meaning might depend on
     * 		specific subclass implementation
     */
    virtual void discretize_surface(
	    std::vector<std::vector<double>>* positions,
	    double resolution) const = 0;
    /** \brief Discretize the surface into a list of points with normal vectors
     *
     * This function returns an (unordered) list of points, including their
     * exact positions and a local surface-orthogonal vector: for instance
     * {{{x1,y1,z1}, {o1,p1,q1}}, {{x2,y2,z2}, {o2,p2,q2}}, ...},
     * where {x,y,z} are position vectors and {o,p,q} describes the direction
     * orthogonal to the surface at this position (pointing away from the wall).
     * The norm of {o,p,q} is the associated area / weight that should be taken
     * into account in an integral over this surface.
     *
     * \param distribution: pointer to return value
     * \param resolution: spatial resolution; exact meaning might depend on
     * 		specific subclass implementation
     */
    virtual void discretize_surface(
	    std::vector<std::vector<std::vector<double>>>* distribution,
	    double resolution) const = 0;
    /** \brief Calculate distance to the boundary in the given direction.
     *
     * \return -1 if there is no boundary in the given direction, 0 if the given
     * position itself is outside the boundaries.
     *
     * \param direction: index of axis (0,1,...,dim-1) along which the distance
     * 		should be calculated
     * \param forward: if true, the search goes in positive direction along the
     * 		specified axis, else in opposite direction
     */
    virtual double distance_directed(std::vector<double>& position,
		    size_t direction, bool forward) const = 0;
    /** \brief Return distance to the nearest boundary point.
     *
     * \return zero if given point is outside boundaries. Returns upper_limit if
     * the distance to the boundary is larger than that value.
     *
     * \param resolution (version 1): In case that the distance cannot be
     * 		calculated exactly, the surface is approximated by a discrete set of points. This
     * 		parameter defines the spatial resolution.
     * \param surface_points (version 2): user-defined set of surface points to
     * 		be used for the calculation
     * \param upper_limit: might be used to reduce the amount of computation
     * 		needed if only distances below a certain limit are of interest
     */
    double distance_minimal(std::vector<double>& position, double resolution,
		    double upper_limit) const;
    double distance_minimal(std::vector<double>& position,
		    std::vector<std::vector<double>>& surface_points,
	            double upper_limit) const;
    /** \brief Return the minimal boundary distances of each point in a
     * discretized grid.
     *
     * \param result: pointer to return value; needs to be initialized with
     * 		desired grid dimensions
     * \param resolution: spatial resolution of surface (see distance_minimal)
     * \param upper_limit: upper limit on distance (see distance_minimal) - By
     * 		default, there is no limit except for the finite system size or
     * 		periodic boundary conditions.
     */
    virtual void minimal_distances(DataFrame<dim, double>* result,
		    double resolution, double upper_limit = -1.);
    /** \brief Calculate the (exponentiated) external potential for hard spheres
     *
     * \param system: specifies all particle species' properties
     * \param potential: pointer to return value; needs to be initialized with
     * 		desired grid dimensions
     * \param resolution: spatial resolution of surface discretization, in units
     * 		of the smallest ion diameter
     */
    void exp_external_potential_hs(System<dim>& system,
    	std::vector<DataFrame<dim, double>>* potential, double resolution = 0.25);
  protected:
    /** \brief System dimensions
     */
    std::vector<double> system_lengths;
    // periodic boundary conditions not always relevant?
//    /** \brief Periodic boundary conditions
//     */
//    std::vector<bool> periodic_boundaries;
    /** \brief Type of the boundary value
     * 
     * Per default, there is no boundary value. This is changed only when
     * a parameter boundary_value is defined in the Constructor's Properties or when
     * set_boundary_values is called. It is relevant for the Poisson solver.
     */
    Type type = None;
    /** \brief Type of the boundary surface: Hard wall or permeable for particles?
     *
     * Per default, boundary surfaces are defined as hard walls. This can be
     * deactivated via an optional parameter in the constructor. It is relevant
     * for the exp_external_potential_hs function.
     */
    bool hard_walls = true;
    /** \brief Boundary value on the surface (electrostatic potential)
     */
    double boundary_value;
    /** \brief Initialize properties from System object
     *
     * This only extracts the standard properties (system dimensions, boundary
     * value...), not the subclass-specific properties.
     *
     * \param index: specifies which boundary value should be chosen, if several
     * 		different surfaces exist in the system
     */
    void extract_standard_properties(System<dim>& system, size_t index);
    /** \brief Calculate bin sizes for a certain discretization
     */
    std::vector<double> bins(std::vector<size_t>& grid_counts) const;
};
#endif // SRC_BOUNDARY_SURFACE_HPP_
