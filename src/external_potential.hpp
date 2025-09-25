//
//
//#ifndef SRC_EXTERNAL_POTENTIAL_HPP_
//#define SRC_EXTERNAL_POTENTIAL_HPP_
///** \file external_potential.hpp
// * \brief External potential functions
// *
// * This  file contains functions to calculate external potentials
// * (hard sphere repulsion / electrostatic) for any BoundaryStructure.
// */
//#include "boundary_structure.hpp"
//#include "poisson_solver.hpp"
//#include <vector>
//#include "properties.hpp"
//#include "data_frame.hpp"
//
///** \brief Calculate the (exponentiated) external potential for hard spheres
// *
// * This function uses the minimal_distances method of the BoundaryStructure
// * class to find the excluded volume.
// */
//template <size_t dim>
//void exp_external_potential_hs(BoundaryStructure<dim>& boundary_structure,
//		std::vector<Properties>& species_properties,
//		std::vector<size_t>& affected_species,
//		std::vector<DataFrame<dim, double>>* potential);
///** \brief Exponentiated external HS potential from multiple surface objects
// */
//template <size_t dim>
//void exp_external_potential_hs(
//		std::vector<BoundaryStructure<dim>>& boundary_structures,
//		std::vector<Properties>& species_properties,
//		std::vector<size_t>& affected_species,
//		std::vector<DataFrame<dim, double>>* potential);
//
///** \brief Calculate the electrostatic potential
// *
// * Here the boundary structure and surface potential are defined within the
// * PoissonSolver which is given as function argument.
// */
//template <size_t dim>
//void electrostatic_potential(
//		PoissonSolver& poisson_solver,
//		DataFrame<dim, double>* potential);
//#endif // SRC_EXTERNAL_POTENTIAL_HPP_
