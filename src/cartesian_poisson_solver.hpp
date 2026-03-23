// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CARTESIAN_POISSON_SOLVER_HPP_
#define SRC_CARTESIAN_POISSON_SOLVER_HPP_
/** \file cartesian_poisson_solver.hpp
 *  \brief Header file for the CartesianPoissonSolver class.
 *
 *  The file contains the declarations of the CartesianPoissonSolver class.
 */
// Includes
#include <cstddef>
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "sparse_matrix.hpp"  // NOLINT
#include "poisson_solver_cartesian.hpp"  // NOLINT
/** \brief This class contains tools to solve the cartesian (3D)
 *         Poisson equation with planar boundaries
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a sparse matrix, that mostly contains zeros.
 *  These kind of matrices can be solved via the GMRES algorithm rather
 *  efficiently.
 *
 *  Note: Since the introduction of the CartesianPoissonSolverAny, the usage of
 *  this class is not recommended except for demonstration and testing purposes.
 *  For solving problems with planar boundaries, the faster PlanarPoissonSolver
 *  is recommended. For other boundary shapes, the CartesianPoissonSolverAny
 *  or RadialPoissonSolver must be used.
 */
class CartesianPoissonSolver : public PoissonSolverCartesian,
	public SparseMatrix {
 public:
  /** \brief Empty Constructor
   */
  CartesianPoissonSolver();
  /** \brief Proper Constructor
   *
   *  \param bin_count grid points in all three dimensions
   *  \param bin_size is the size of the bins in the three dimensions
   *  \param periodic_boundaries contains wether a dimension has PBC
   */
  CartesianPoissonSolver(
      std::vector<size_t> bin_count, std::vector<double> bin_size,
      std::vector<bool> periodic_boundaries);
  /** \brief Destructor
   */
  ~CartesianPoissonSolver();
  /** \brief Solve the linear equation system
   *  
   *  The GMRES solver of the SparseMatrix class is used to find an iterative
   *  solution. The boundary values are by default zero.
   */
  void solve(
      std::vector<double>& rhs,
      std::vector<double>& solution) override;
  /** \brief Solve the linear equation system with given boundary values
   *
   * There are six boundary values on the six faces of the system. They are
   * stored in vector format as following:
   * boundary_values = {{x_left, x_right}, {y_left, y_right}, {z_left, z_right}}
   */
  void solve(
      std::vector<double>& rhs,
      std::vector<std::vector<double>> boundary_values,
      std::vector<double>& solution);
  void solve(
      DataFrame<3, double>& rhs,
      std::vector<std::vector<double>> boundary_values,
      DataFrame<3, double>& solution);

 private:
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
  /** \brief Set the Laplace matrix with boundary conditions
   */
  void set_laplacian();
  /** \brief Set the Laplace matrix without boundary conditions
   */
  void set_bare_laplacian();
  /** \brief Modify the Laplace matrix according to the periodic_boundaries
   */
  void set_boundary_conditions();
  /** \brief From a coordinate index calculate the equation index
   *
   *  \return false if index out of bounds
   */
  bool coordinates_to_index(size_t i, size_t j, size_t k, size_t* index);
  /** \brief Add the boundary values to a given right-hand side vector
   */
  void add_boundary_values(
      std::vector<std::vector<double>>& boundary_values,
      std::vector<double>& rhs);
  /** \brief Subtract the boundary values from a given right-hand side vector
   */
  void remove_boundary_values(
      std::vector<std::vector<double>>& boundary_values,
      std::vector<double>& rhs);
};
#endif  // SRC_CARTESIAN_POISSON_SOLVER_HPP_
