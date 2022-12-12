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
#include "sparse_matrix.hpp"  // NOLINT
/** \brief This class contians tools to solve the cartesian (3D)
 *         poisson equation
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a sparse matrix, that mostly contains zeros.
 *  These kind of matrices can be solved via the GMRES algorithm rather
 *  efficiently.
 */
class CartesianPoissonSolver : public SparseMatrix {
 public:
  /** \brief Empty Constructor
   */
  CartesianPoissonSolver();
  /** \brief Proper Constructor
   * 
   * Passes the "size" of the array to be solved, the grid spacing "dr" of
   * the array.
   */
  CartesianPoissonSolver(
      std::vector<size_t> size, double bin_size,
      std::vector<bool> periodic_boundaries);
  /** \brief Destructor
   */
  ~CartesianPoissonSolver();
  /** \brief Solve the linear equation system
   *  
   *  The GMRES solver of the SparseMatrix class is used to find an iterative
   *  solution.
   */
  void solve(
      std::vector<double>& rhs,
      std::vector<std::vector<double>> boundary_values,
      std::vector<double>& solution);

 private:
  /** \brief Number of bins in every dimension
   */
  std::vector<size_t> bin_count;
  /** \brief Bin size
   */
  double bin_size;
  /** \brief Bin size squared
   */
  double bin_size_squared;
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
