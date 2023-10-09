// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
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
/** \brief This class contians tools to solve the cartesian (3D)
 *         poisson equation
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a sparse matrix, that mostly contains zeros.
 *  These kind of matrices can be solved via the GMRES algorithm rather
 *  efficiently.
 */
class CartesianPoissonSolverAny : public SparseMatrix {
 public:
  /** \brief Empty Constructor
   */
  CartesianPoissonSolverAny();
  /** \brief Proper Constructor
   *
   *  \param bin_count grid points in all three dimensions
   *  \param bin_size is the size of the bins in the three dimensions
   *  \param pariodic_boundaries contains wether a dimension has PBC
   *  \param boundary_positions see set_boundary_values()
   *  \param boundary_values see set_boundary_values()
   */
  CartesianPoissonSolverAny(
      std::vector<size_t> bin_count,
      std::vector<double> bin_size,
      std::vector<bool> periodic_boundaries,
      std::vector<std::vector<std::vector<std::vector<size_t>>>>& boundary_positions,  // NOLINT
      std::vector<std::vector<std::vector<double>>>& boundary_values);
  /** \brief Destructor
   */
  ~CartesianPoissonSolverAny();
  /** \brief Solve the linear equation system
   *  
   *  The GMRES solver of the SparseMatrix class is used to find an iterative
   *  solution.
   */
  void solve(std::vector<double>& rhs, std::vector<double>& solution);

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
  /** \brief Holds the rows to remove due to the boundary conditions
   *
   *  If boundaries inside the system are used, certain rows in the matrix
   *  become redundant and are hnce removed. This vector remembers them.
   */
  std::vector<size_t> rows_to_remove;
  /** \brief Holds the columns to remove due to the boundary conditions
   *
   *  If boundaries inside the system are used, certain columns in the matrix
   *  become need to be removed, because their potential value is already known.
   *  Column position and boundary position coincide.
   *  This vector remembers the removed column positions and the corresponding
   *  boundary value.
   */
  std::vector<std::pair<size_t, double>> cols_to_remove;
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
  /** \brief From a coordinate index calculate the equation index
   *
   *  \return false if index out of bounds
   */
  bool coordinates_to_index(size_t i, size_t j, size_t k, size_t* index);
  bool coordinates_to_index(std::vector<size_t> pos, size_t* index);
  /** \brief Set the boundary position and values
   *
   *  Since non-flat (any) boundaries are quite complicated to set, it is done
   *  only if this function is called externally.
   *
   *  The 3D box contains 6 faces. For each face there is a 2D grid. It is
   *  possible to set one boundary value per face per 2D grid point.
   *  This is the same as setting two boundary values per direction per 2D grid
   *  point.
   *
   *  \param boundary_positions is a vector of a vector of pairs of position
   *         vectors e.g.
   *         (
   *         ((())),
   *         (((l1,m1,n1), (l1,m2,n1)), ..., ((l10,m19,n10), (l10,m10,n20))),
   *         (((i1,j1,k1), (i1,j1,k2)), ..., ((i10,j10,k19), (i10,j10,k20)))
   *         )
   *         would set 20 boundary values in the y and z-direction.
   *         The outermost vector stands for the three directions, the next
   *         vector contains all the pairs of boundary positions in that
   *         direction.
   */
  void set_boundary_values(
      std::vector<std::vector<std::vector<std::vector<size_t>>>>& boundary_positions,  // NOLINT
      std::vector<std::vector<std::vector<double>>>& boundary_values);
  /** \brief Add the boundary values to a given right-hand side vector
   */
  void add_boundary_values(std::vector<double>& rhs);
  /** \brief Subtract the boundary values from a given right-hand side vector
   */
  void remove_boundary_values(std::vector<double>& rhs);
};
#endif  // SRC_CARTESIAN_POISSON_SOLVER_ANY_HPP_
