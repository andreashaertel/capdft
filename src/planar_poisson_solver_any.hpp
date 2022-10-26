// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_PLANAR_POISSON_SOLVER_ANY_HPP_
#define SRC_PLANAR_POISSON_SOLVER_ANY_HPP_
/** \file planar_poisson_solver_any.hpp
 *  \brief Header file for the PlanarPoissonSolverAny class.
 *
 *  The file contains the class declarations of the PlanarPoissonSolverAny
 *  class.
 */
// Includes
#include "banded_matrix.hpp"
#include <vector>
/** \brief This class represents a special banded matrix of a Laplace operator
 *         with boundary conditions that can by ANYwhere.
 *
 *  It is derived from the BandedMatrix class and represents a pentadiagonal
 *  matrix.
 */
class PlanarPoissonSolverAny : public BandedMatrix {
 public:
  /** \brief Empty constructor
   */
  PlanarPoissonSolverAny();
  /** \brief Constructor
   *
   *  The constructor of BandedMatrix is called with the band widths set to 2
   *  (pentadiagonal matrix) and with 2 less bins. Since we use boundaries
   *  somewhere inside the system, we need to replace two bins with the
   *  boundaries, which leads to "size - 2" bins.
   *
   *  \param dz is the bin size.
   *  \param size is the size of the density array.
   */
  PlanarPoissonSolverAny(
      size_t size, double dz, std::vector<size_t> boundary_positions);
  /** \brief Destructor
   */
  ~PlanarPoissonSolverAny();
  /** \brief Solve the Poisson equation
   *
   *  Solve the Poisson equation using the underlying solving algorithm of
   *  banded matrices of the base class.
   *  \param rhs is a pointer to the data of the right hand side of the Poisson
   *         equation.
   *  \param solution is a pointer to the array to which the solution of the
   *         Poisson equation is written.
   */
   void solve(
      double* rhs, std::vector<double> boundary_values, double* solution);

 private: 
  /** \brief Boundary positions
   */
  std::vector<size_t> boundary_positions;
  /** \brief Bin size
   */
  double dz;
  /** \brief Bin size squared
   */
  double dzdz;
  /** \brief Set the matrix elements according to the laplacian with any
   *         boundary position
   */
  void set_laplacian();
  /** \brief The boundary conditions are added to the right hand side of the
   *         Poisson equation
   */
  void add_boundary_values(double* rhs, std::vector<double> boundary_values);
  /** \brief The boundary conditions are subtracted from the right hand side of
   *         the Poisson equation
   */
  void remove_boundary_values(double* rhs, std::vector<double> boundary_values);
  /** \brief Rearrange the elements of the solution array
   *
   *  The Poisson equation with n elements corresponds to a banded matrix
   *  equation system with 2 less rows. This function rearranges the solution
   *  array, such that it fits the Poisson equation solution.
   */
  void rearrange_solution(double* rhs, std::vector<double> boundary_values);
};
#endif  // SRC_PLANAR_POISSON_SOLVER_ANY_HPP_
