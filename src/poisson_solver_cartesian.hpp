// SPDX-FileCopyrightText: 2025 Fabienne Dressler <fab.dressler@web.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_POISSON_SOLVER_CARTESIAN_HPP_
#define SRC_POISSON_SOLVER_CARTESIAN_HPP_
/** \file poisson_solver_cartesian.hpp
 * \brief Header file for the PoissonSolverCartesian class.
 *
 *  The file contains the declarations of the PoissonSolverCartesian class.
 */
// Includes
#include <vector>
#include <cstddef>

/** \brief Abstract class for the numerical solution of the Poisson equation
 *
 * This class acts as a superclass for the CartesianPoissonSolver (planar
 * boundaries with homogeneous boundary potential) and the
 * CartesianPoissonSolverAny (more generally structured boundaries). This is 
 * useful e.g. for the generalized FunctionalESMFCartesian, which should work
 * with any of the two PoissonSolverCartesians.
 */
class PoissonSolverCartesian {
 public:
  // Constructor is specific to subclasses
  /** \brief Destructor
   */
  virtual ~PoissonSolverCartesian() {};
  /** \brief Solve the linear equation system
   *  
   *  The GMRES solver of the SparseMatrix class is used to find an iterative
   *  solution.
   */
  virtual void solve(std::vector<double>& rhs, std::vector<double>& solution) = 0;
  // Private member functions and variables are specific to subclasses
};
#endif // SRC_POISSON_SOLVER_CARTESIAN_HPP_
