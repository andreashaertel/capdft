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
// Flag definition for boundary conditions
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
   * the array and the left boundary condition position "inner_distance".
   */
  PlanarPoissonSolver(size_t dim, double dz);
  /** \brief Destructor
   */
  ~PlanarPoissonSolver();
  /** \brief Set the Laplace matrix with the chosen boundary conditions
   */
  void set_laplacian(BoundaryFlag flag);
  /** \brief Solve the linear equation system: the boundary must be specified
   *  
   *  For example the DIRICHLET_NEUMANN case expects a left boundary value
   *  and a right first derivative value of the sought function.
   *  The DIRICHLET_DIRICHLET case expects a left and right boundary value of
   *  the sought function.
   */
  void solve(double leftBoundaryValue, double rightBoundaryValue, double* rhs,
      double* solution);
  /** \brief Returns the size of the tridiagonal matrix
   */
  size_t size();

 private:
  /** \brief Number of rows/columns of the square tridiagonal matrix
   */
  size_t dim;
  /** Radial bin size (not needed in carthesian case)
   */
  double dz;
  /** \brief Holds the specified boundary conditions
   */
  BoundaryFlag flag;
  /** \brief Boundary condition dependent positions of the first/last bin
   * 
   *  It is either (dr or dr/2) and determined automatically
   */
  double shiftLeft, shiftRight;
  /** \brief Upper off-diagonal elements of the tridiagonal matrix
   */
  double* upper;
  /** \brief Diagonal elements of the tridiagonal matrix
   */
  double* diag;
  /** \brief Lower off-diagonal elements of the tridiagonal matrix
   */
  double* lower;
  /** \brief Set the Laplace matrix with certain boundary conditions
   *
   *  For NN the boundaries should be half a grid point away from the start/end
   *  point (x_{Wall}=x_{-0.5}|x_0|x_1|...|x_N|x_{N+0.5}=x_{Wall}).
   *  For DD the boundaries should be a full grid point away from the start/end
   *  point (x_{Wall}=x_{-1}|x_0|x_1|...|x_N|x_{N+1}=x_Wall).
   */
  void set_laplacian_NN();
  void set_laplacian_DD();
  void set_laplacian_ND();
  void set_laplacian_DN();
  /** \brief Set the Laplace matrix without boundary conditions
   */
  void set_laplacian();
  /** \brief Add the boundary values to the right-hand side
   */
  void set_boundary_values_laplace_radial(double leftBoundaryValue,
      double rightBoundaryValue, double* rhs);
};
#endif  // SRC_PLANAR_POISSON_SOLVER_HPP_
