// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_BOUNDARIES_HPP_
#define SRC_BOUNDARIES_HPP_
/** \file boundaries.hpp
 * \brief Header file for the Boundaries class.
 *
 *  The file contains the declarations of the CartesianPoissonSolverAny class.
 */
// Includes
#include <vector>
#include <list>
#include "properties.hpp"
#include "boundary_surface.hpp"
#include "system.hpp"

/** \brief This class is a container for multiple surfaces in one system.
 *
 * This class is derived from an std::list, such that it can be handled like a
 * list for instance when accessing the elements. Note that the elements are
 * stored by pointers, not by value.
 *
 * Making sure that the surface objects put into this container are consistent
 * (e.g. having the same system lengths and periodic boundary conditions) lies
 * in the responsibility of the user.
 */
template <size_t dim>
class Boundaries : public std::list<BoundarySurface<dim>*> {
  public:
    /** \brief Construct empty container
     */
    Boundaries();
//    /** \brief Construct and initialize from System object
//     *
//     * \param system: must contain indexed properties for the different
//     * 		BoundarySurface objects to be initialized here, as well as an
//     * 		indexed std::string-property "type" specifying the corresponding
//     * 		BoundarySurface-subclass
//     * \param indices: which sets of boundary properties to choose
//     */
//    Boundaries(System<dim>& system,
//		    std::vector<size_t>& indices); // implementation missing
//    /** \brief Destructor
//     */
//    ~Boundaries();
    /** \brief Add a surface object to the container
     *
     * Note that, since the object is given as a pointer, it has to be defined
     * globally - i.e. one can't construct it locally, hand it to the Boundaries
     * object via add_surface and then forget about it. The surface object must 
     * remain intact in the main code as long as the Boundaries object is used.
     */
    void add_surface(BoundarySurface<dim>* surface);
//    /** \brief Add a surface object to the container
//     *
//     * Parameters analogous to the constructor
//     */
//    void add_surface(System<dim>& system, size_t index);
    /** \brief Set uniform boundary value over all surfaces
     */
    void set_all_boundary_values(double value);
    /** \brief Calculate total external potential
     *
     * \param resolution: If value is larger than zero, this specifies the
     * 		surface resolution (in units of the smallest ion diameter) as
     * 		defined in BoundarySurface::exp_external_potential. Else, the
     * 		default resolution settings are used.
     */
    void exp_external_potential_hs(System<dim>& system,
		    std::vector<DataFrame<dim, double>>* result,
		    double resolution = 0.);
};
#endif // SRC_BOUNDARIES_HPP_
