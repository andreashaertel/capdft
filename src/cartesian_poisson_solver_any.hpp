// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>,
// 			   2026 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CARTESIAN_POISSON_SOLVER_ANY_HPP_
#define SRC_CARTESIAN_POISSON_SOLVER_ANY_HPP_
/** \file cartesian_poisson_solver_any.hpp
 *  \brief Header file for the CartesianPoissonSolverAny class.
 *
 *  The file contains the declarations of the CartesianPoissonSolverAny class.
 */
// Includes
#include <cstddef>
#include <vector>
#include "sparse_matrix.hpp"  // NOLINT
#include "boundary_surface.hpp"
#include "poisson_solver_cartesian.hpp"
#include "boundaries.hpp"
#include "system.hpp"
/** \brief This class contains tools to solve the cartesian (3D)
 *         poisson equation with structured boundaries
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a sparse matrix, that mostly contains zeros.
 *  These kind of matrices can be solved via the GMRES algorithm rather
 *  efficiently.
 */
class CartesianPoissonSolverAny : public PoissonSolverCartesian,
	public SparseMatrix {
 public:
  /** \brief Empty Constructor
   */
  CartesianPoissonSolverAny();
  /** \brief Constructor with one BoundarySurface object
   *
   *  \param system defines the spatial dimensions and periodic boundary cond.
   *  \param boundary_surface specifies the boundary values and positions
   */
  CartesianPoissonSolverAny(System<3>& system,
      BoundarySurface<3>* boundary_surface);
  /** \brief Constructor with multiple BoundarySurface objects
   *
   *  \param system defines the spatial dimensions and periodic boundary cond.
   *  \param boundaries specifies the boundary values and positions
   *
   *  Warning: Errors might occur if the BoundarySurface objects in boundaries
   *  have overlaps. For correct usage, make sure they are spatially separated.
   */
  CartesianPoissonSolverAny(System<3>& system, Boundaries<3>* boundaries);
  /** \brief Destructor
   */
  ~CartesianPoissonSolverAny();
  /** \brief Solve the linear equation system
   *  
   *  The GMRES solver of the SparseMatrix class is used to find an iterative
   *  solution.
   *
   *  \param rhs: right hand side of the Poisson equation: charge density
   *  		profile, rescaled by -4*PI*bjerrum length
   *  \param solution: pointer to return value. Initial state is used as initial
   *  		guess.
   */
  void solve(std::vector<double>& rhs, std::vector<double>& solution) override;
  /** \brief Solve the linear equation system
   *  
   *  Analogous to above version, but with DataFrames instead of vectors.
   */
  void solve(DataFrame<3, double>& rhs, DataFrame<3, double>& solution);

 private:
  /** \brief System object, specifying e.g. the dimensions and PBCs.
   */
  System<3> system;
  /** \brief Number of bins in every dimension
   */
  std::vector<size_t> bin_count;
  /** \brief Bin size
   */
  std::vector<double> bin_size;
  /** \brief Bin size squared
   */
  std::vector<double> bin_size_squared;
  /** \brief Holds the specified boundary conditions
   */
  std::vector<bool> periodic_boundaries;
  /** \brief Holds the indices of rows & cols to remove due to the boundaries
   *
   *  If boundaries inside the system are used, certain rows in the matrix
   *  become redundant and are hence removed; the corresponding columns also
   *  need to be removed because their potential value is already known.
   *  These positions coincide with the boundary positions.
   *  This vector remembers their positions.
   */
  std::vector<size_t> boundary_positions;
  /** \brief Holds the indices to remove and their associated boundary values
   *
   *  This vector remembers not only the removed row and column indices like
   *  boundary_position, but also the corresponding boundary value.
   */
  std::vector<std::pair<size_t, double>> boundary_points;
  /** \brief Holds the values that need to be added to the rhs due to the
   *         boundary conditions
   */
  std::vector<double> rhs_addition;
  /** \brief Set the Laplace matrix with boundary conditions
   */
  void set_laplacian();
  /** \brief Set the Laplace matrix without boundary conditions
   */
  void set_bare_laplacian();
  /** \brief Modify the Laplace matrix according to the periodic_boundaries
   */
  void set_boundary_conditions();
  /** \brief From coordinate indices calculate the equation index
   *
   *  \return false if index out of bounds (even if there are PBC)
   */
  bool coordinates_to_index(size_t i, size_t j, size_t k, size_t* index);
  bool coordinates_to_index(std::vector<size_t> pos, size_t* index);
  /** \brief From coordinate indices calculate the position
   */
  void coordinates_to_position(std::vector<size_t> pos, std::vector<double>* position);
  /** \brief Set the boundary position and values
   */
  void set_boundary_values(BoundarySurface<3>* boundary_surface);
  void set_boundary_values(Boundaries<3>* boundaries);
  /** \brief Calculate the modifications of right hand side and matrix elements
   * needed to set the boundary positions and values of one boundary object
   */
  void calc_boundary_values(BoundarySurface<3>* boundary_surface);
  /** \brief Calculate the modification induced by a given boundary point
   *
   * This calculates all modifications of the right hand side and matrix
   * elements induced by the given boundary point.
   */
  void modify_boundary_point(BoundarySurface<3>* boundary_surface,
		std::vector<size_t>& index_position);	
  /** \brief Calculate the modification affecting a given grid point
   *
   * This calculates the modifications of all right hand side and matrix
   * elements associated with the given grid point.
   */
  void modify_grid_point(BoundarySurface<3>* boundary_surface,
		std::vector<size_t>& index_position);	
  /** \brief Add the boundary values to a given right-hand side vector
   */
  void add_boundary_values(std::vector<double>& rhs);
  /** \brief Subtract the boundary values from a given right-hand side vector
   */
  void remove_boundary_values(std::vector<double>& rhs);
};
#endif  // SRC_CARTESIAN_POISSON_SOLVER_ANY_HPP_
