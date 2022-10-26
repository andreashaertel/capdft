// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_PLANAR_POISSON_SOLVER_HPP_
#define SRC_PLANAR_POISSON_SOLVER_HPP_
/** \file planar_poisson_solver.hpp
 *  \brief Header file for the PlanarPoissonSolver class.
 *
 *  The file contains the declarations of the PlanarPoissonSolver class.
 */
// Includes
#include <cstddef>
#include <vector>
/** \brief Flag definition for boundary conditions
 */
enum BoundaryFlag {
  NEUMANN_NEUMANN = 0,
  DIRICHLET_DIRICHLET = 1,
  NEUMANN_DIRICHLET = 2,
  DIRICHLET_NEUMANN = 3
};
/** \brief This class contians tools to solve the planar (1D) poisson equation
 * 
 *  The numerical Poisson equation can be rewritten into a matrix equation
 *  containing a tridiagonal matrix. These matrices are called "sparse",
 *  because they mostly contain zeros.
 *  The only non-zero elements are: diagonal elements A(i,i), upper diagonal
 *  elements A(i,i+1), lower diagonal elements A(i,i-1).
 *  Hence this matrix can be solved easily and stored in three arrays.
 */
class PlanarPoissonSolver {
 public:
  /** \brief Empty Constructor
   */
  PlanarPoissonSolver();
  /** \brief Proper Constructor
   * 
   * Passes the size "dim" of the array to be solved, the grid spacing "dr" of
   * the array.
   */
  PlanarPoissonSolver(size_t dim, double dz, BoundaryFlag flag);
  /** \brief Destructor
   */
  ~PlanarPoissonSolver();
  /** \brief Solve the linear equation system: the boundary must be specified
   *  
   *  For example the DIRICHLET_NEUMANN case expects a left boundary value
   *  and a right first derivative value of the sought function.
   *  The DIRICHLET_DIRICHLET case expects a left and right boundary value of
   *  the sought function.
   *
   *  The boundary is assumed to be half a grid point away from the first/last
   *  grid point (z_{left boundary}=-0.5 * dz).
   *
   *  Solving the tridiagonal matrix equation is done in the following way.
   *  One starts by eliminating solution[0] by subtracting the first from the
   *  second line. One ends up with the same situation and does the whole thing
   *  again with solution[1] and the second and the third line. And so on ...
   *  This obviously fails for NEUMANN_NEUMANN, because the system is determined
   *  up to a constant and thus the system of equation is underdetermined.
   */
  void solve(double left_boundary_value, double right_boundary_value,
      double* rhs, double* solution);
  /** \brief Returns the size of the tridiagonal matrix
   */
  size_t size();

 private:
  /** \brief Number of rows/columns of the square tridiagonal matrix
   */
  size_t dim;
  /** \brief Bin size
   */
  double dz;
  /** \brief Bin size squared
   */
  double dzdz;
  /** \brief Holds the specified boundary conditions
   */
  BoundaryFlag flag;
  /** \brief Upper off-diagonal elements of the tridiagonal matrix
   */
  std::vector<double> upper;
  /** \brief Diagonal elements of the tridiagonal matrix
   */
  std::vector<double> diag;
  /** \brief Lower off-diagonal elements of the tridiagonal matrix
   */
  std::vector<double> lower;
  /** \brief Check if the boundary positions are valid
   */
  bool check_boundary_positions();
  /** \brief Set the Laplace matrix with certain boundary conditions
   */
  void set_laplacian_NN();
  void set_laplacian_DD();
  void set_laplacian_ND();
  void set_laplacian_DN();
  /** \brief Set up arrays for needed to describe the tridiagonal matrix
   */
  void initialize_matrix();
  /** \brief Set the Laplace matrix with boundary conditions
   */
  void set_laplacian();
  /** \brief Set the Laplace matrix without boundary conditions
   */
  void set_bare_laplacian();
  /** \brief Add the boundary values to the right-hand side
   */
  void set_boundary_values_laplace(double left_boundary_value,
      double right_boundary_value, double* rhs);
  /** \brief Solution algorithm for the tridiagonal matrix equation
   */
  void solve(double* rhs, double* solution);
};
#endif  // SRC_PLANAR_POISSON_SOLVER_HPP_
