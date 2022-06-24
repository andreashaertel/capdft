// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file radial_poisson_solver.cpp
 *  \brief Source file for the RadialPoissonSolver class.
 *
 *  The file contains the definitions of the RadialPoissonSolver class.
 */
#include "radial_poisson_solver.hpp"  // NOLINT
#include <iostream>
#include <vector>
// _____________________________________________________________________________
RadialPoissonSolver::RadialPoissonSolver() {
  //
}
// _____________________________________________________________________________
RadialPoissonSolver::RadialPoissonSolver(
    size_t dim, double dr, double inner_distance
): dim(dim), dr(dr), inner_distance(inner_distance) {
  this->drdr = dr * dr;
  // Allocate memory for the matrix elements
  upper = new double[dim];
  diag = new double[dim];
  lower = new double[dim];
}
// _____________________________________________________________________________
RadialPoissonSolver::~RadialPoissonSolver() {
  // Free memory
  delete [] upper;
  delete [] diag;
  delete [] lower;
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian(BoundaryFlag flag) {
  this->flag = flag;
  // Set boundary conditions
  switch (flag) {
    case NEUMANN_NEUMANN:
      this->shiftLeft = 0.5;
      this->shiftRight = 0.5;
      set_radial_laplacian_NN();
      break;
    case DIRICHLET_DIRICHLET:
      this->shiftLeft = 1.;
      this->shiftRight = 1.;
      set_radial_laplacian_DD();
      break;
    case NEUMANN_DIRICHLET:
      this->shiftLeft = 0.5;
      this->shiftRight = 1.;
      set_radial_laplacian_ND();
      break;
    case DIRICHLET_NEUMANN:
      this->shiftLeft = 1.;
      this->shiftRight = 0.5;
      set_radial_laplacian_DN();
      break;
    default:
      std::cerr << "set_radial_laplacian(): \"ERROR: Invalid flag!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian() {
  double r = 0.;
  // Just set the matrix elements without any boundary conditions
  // First the regular second derivative
  for (size_t i = 0; i < dim; ++i) {
    upper[i] = 1. / drdr;
    diag[i] = -2. / drdr;
    lower[i] = 1. / drdr;
  }
  // Then add the first derivative divided by r (spherical coordinates).
  for (size_t i = 0; i < dim; ++i) {
    r = inner_distance + (static_cast<double>(i) + shiftLeft) * dr;
    upper[i] += 1. / (r * dr);
    diag[i] += 0.;
    lower[i] += -1. / (r * dr);
  }
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian_NN() {
  double r = 0.;
  // Set Laplacian without boundary conditions
  set_radial_laplacian();
  // First set the regular second derivative at the boundaries.
  upper[0] = 1. / drdr;
  diag[0] = -1. / drdr;
  lower[dim-1] = 1. / drdr;
  diag[dim-1] = -1. / drdr;
  // Then add the first derivative divided by r at the boundaries.
  r = inner_distance + shiftLeft * dr;
  upper[0] += 1. / (r * dr);
  diag[0] += -1. / (r * dr);
  r = inner_distance + (static_cast<double>(dim-1) + shiftLeft) * dr;
  lower[dim-1] += -3. / (r * dr);
  diag[dim-1] += 3. / (r * dr);
  // These matrix element do not exist.
  lower[0] = 0.;
  upper[dim-1] = 0.;
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian_DD() {
  // Set Laplacian without boundary conditions
  set_radial_laplacian();
  // The standard Laplacian corresponds to the DD case.
  // These matrix element do not exist.
  lower[0] = 0.;
  upper[dim-1] = 0.;
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian_ND() {
  double r = 0.;
  // Set Laplacian without boundary conditions
  set_radial_laplacian();
  // First set the regular second derivative at the boundaries.
  upper[0] = 1. / drdr;
  diag[0] = -1. / drdr;
  // Then add the first derivative divided by r at the boundaries.
  r = inner_distance + shiftLeft * dr;
  upper[0] += 1. / (r * dr);
  diag[0] += -1. / (r * dr);
  // These matrix element do not exist.
  lower[0] = 0.;
  upper[dim-1] = 0.;
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_radial_laplacian_DN() {
  double r = 0.;
  // Set Laplacian without boundary conditions
  set_radial_laplacian();
  // First set the regular second derivative at the boundaries.
  lower[dim-1] = 1. / drdr;
  diag[dim-1] = -1. / drdr;
  // Then add the first derivative divided by r at the boundaries.
  r = inner_distance + (static_cast<double>(dim-1) + shiftLeft) * dr;
  lower[dim-1] += -3. / (r * dr);
  diag[dim-1] += 3. / (r * dr);
  // These matrix element do not exist.
  lower[0] = 0.;
  upper[dim-1] = 0.;
}
// _____________________________________________________________________________
void RadialPoissonSolver::solve(
    double leftBoundaryValue, double rightBoundaryValue, double* rhs,
    double* solution
) {
  // Work space vector
  std::vector<double> work(dim);
  double diagTemp;
  size_t index;
  // Set the boundary values
  set_boundary_values_laplace_radial(
      leftBoundaryValue, rightBoundaryValue, rhs);
  // Solving algorithm for tridiagonal matrices:
  // One starts by eliminating solution[0] by subtracting the first from the
  // second line. One ends up with the same situation and does the whole thing
  // again with solution[1] and the second and the third line. And so on ...
  // THIS OBVIOUSLY FAILS FOR NEUMANN-NEUMANN, BECAUSE THE SYSTEM IS DETERMINED
  // UP TO A CONSTANT. THUS THE LINEAR EQUATION SYSTEM IS UNDERDETERMINED.
  diagTemp = diag[0];
  if (diagTemp == 0.) {
    std::cerr << "solve(): \"ERROR: Reduce your matrix by one line!\"";
    std::cerr << std::endl;
    exit(1);
  }
  solution[0] = rhs[0] / diagTemp;
  for (size_t i = 1; i < dim; ++i) {
    work.at(i) = upper[i-1] / diagTemp;
    diagTemp = diag[i] - lower[i] * work.at(i);  // subtraction of the two lines
    if (diagTemp == 0.) {
      std::cerr << "solve(): \"ERROR: Algorithm failed in row ";
      std::cerr << i << "!\"" << std::endl;
      exit(1);
    }
    solution[i] = (rhs[i] - lower[i] * solution[i-1]) / diagTemp;
  }
  // The matrix solves in the last line, because there, there are only two
  // unknowns left instead of three. The solution is obtained by backward
  // substitution.
  for (size_t i = 0; i < dim-1 ; ++i) {
    index = dim - 2 - i;  // go through the arrays backwards
    solution[index] += -work.at(index+1) * solution[index+1];
  }
}
// _____________________________________________________________________________
void RadialPoissonSolver::set_boundary_values_laplace_radial(
    double leftBoundaryValue, double rightBoundaryValue, double* rhs
) {
  double rLeft, rRight;
  rLeft = inner_distance + shiftLeft * dr;
  rRight = inner_distance + (static_cast<double>(dim-1) + shiftLeft) * dr;
  // Add the boundary values to the first and the last element of the
  // right-hand side.
  switch (flag) {
    case NEUMANN_NEUMANN:
      rhs[0] += leftBoundaryValue / dr;
      rhs[0] += -leftBoundaryValue / rLeft;
      rhs[dim-1] += -rightBoundaryValue / dr;
      rhs[dim-1] += rightBoundaryValue / rRight;
      break;
    case DIRICHLET_DIRICHLET:
      rhs[0] += -leftBoundaryValue / drdr;
      rhs[0] += leftBoundaryValue / (rLeft * dr);
      rhs[dim-1] += -rightBoundaryValue / drdr;
      rhs[dim-1] += -rightBoundaryValue / (rRight * dr);
      break;
    case NEUMANN_DIRICHLET:
      rhs[0] += leftBoundaryValue / dr;
      rhs[0] += -leftBoundaryValue / rLeft;
      rhs[dim-1] += -rightBoundaryValue / drdr;
      rhs[dim-1] += -rightBoundaryValue / (rRight * dr);
      break;
    case DIRICHLET_NEUMANN:
      rhs[0] += -leftBoundaryValue / drdr;
      rhs[0] += leftBoundaryValue / (rLeft * dr);
      rhs[dim-1] += -rightBoundaryValue / dr;
      rhs[dim-1] += rightBoundaryValue / rRight;
      break;
    default:
      std::cerr << "set_boundary_values_laplace_radial(): ";
      std::cerr << "\"ERROR: Invalid flag!\"";
      std::cerr << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
size_t RadialPoissonSolver::size() {
  return dim;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
