// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file planar_poisson_solver.cpp
 *  \brief Source file for the PlanarPoissonSolver class.
 *
 *  The file contains the definitions of the PlanarPoissonSolver class.
 */
#include "planar_poisson_solver.hpp"  // NOLINT
#include <iostream>
#include <vector>
// _____________________________________________________________________________
PlanarPoissonSolver::PlanarPoissonSolver() {
  //
}
// _____________________________________________________________________________
PlanarPoissonSolver::PlanarPoissonSolver(
    size_t dim, double dz, BoundaryFlag flag)
  : dim(dim), dz(dz), dzdz(dz * dz), flag(flag) {
  initialize_matrix();
  set_laplacian();
}
// _____________________________________________________________________________
PlanarPoissonSolver::~PlanarPoissonSolver() {
  //
}
// _____________________________________________________________________________
void PlanarPoissonSolver::initialize_matrix() {
  upper = std::vector<double>(dim);
  diag = std::vector<double>(dim);
  lower = std::vector<double>(dim);
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian() {
  // Set boundary conditions
  switch (flag) {
    case NEUMANN_NEUMANN:
      set_laplacian_NN();
      break;
    case DIRICHLET_DIRICHLET:
      set_laplacian_DD();
      break;
    case NEUMANN_DIRICHLET:
      set_laplacian_ND();
      break;
    case DIRICHLET_NEUMANN:
      set_laplacian_DN();
      break;
    default:
      std::cerr << "PlanarPoissonSolver::set_laplacian(): ";
      std::cerr << "\"ERROR: Invalid flag!\"" << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_bare_laplacian() {
  // Just set the matrix elements without any boundary conditions
  for (size_t i = 0; i < dim; ++i) {
    upper[i] = 1. / dzdz;
    diag[i] = -2. / dzdz;
    lower[i] = 1. / dzdz;
  }
  // These matrix elements do not exist, since off-diagonals are of size dim-1.
  lower[0] = 0.;
  upper[dim-1] = 0.;
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_NN() {
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  diag[0] = -1. / dzdz;
  upper[0] = 1. / dzdz;
  diag[dim-1] = -1. / dzdz;
  lower[dim-1] = 1. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_DD() {
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  diag[0] = -4. / dzdz;
  upper[0] = 4. / (3. * dzdz);
  diag[dim-1] = -4. / dzdz;
  lower[dim-1] = 4. / (3. * dzdz);
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_ND() {
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  diag[0] = -1. / dzdz;
  upper[0] = 1. / dzdz;
  diag[dim-1] = -4. / dzdz;
  lower[dim-1] = 4. / (3. * dzdz);
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_DN() {
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  diag[0] = -4. / dzdz;
  upper[0] = 4. / (3. * dzdz);
  diag[dim-1] = -1. / dzdz;
  lower[dim-1] = 1. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::solve(
    double left_boundary_value, double right_boundary_value, double* rhs,
    double* solution) {
  // Set the boundary values
  set_boundary_values_laplace(left_boundary_value, right_boundary_value, rhs);
  // Solve matrix equation
  solve(rhs, solution);
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_boundary_values_laplace(
    double left_boundary_value, double right_boundary_value, double* rhs) {
  // Add the boundary values to the first and the last element of the
  // right-hand side.
  switch (flag) {
    case NEUMANN_NEUMANN:
      rhs[0] += left_boundary_value / dz;
      rhs[dim-1] += -right_boundary_value / dz;
      break;
    case DIRICHLET_DIRICHLET:
      rhs[0] += -8. * left_boundary_value / (3. * dzdz);
      rhs[dim-1] += -8. * right_boundary_value / (3. * dzdz);
      break;
    case NEUMANN_DIRICHLET:
      rhs[0] += left_boundary_value / dz;
      rhs[dim-1] += -8. * right_boundary_value / (3. * dzdz);
      break;
    case DIRICHLET_NEUMANN:
      rhs[0] += -8. * left_boundary_value / (3. * dzdz);
      rhs[dim-1] += -right_boundary_value / dz;
      break;
    default:
      std::cerr << "PlanarPoissonSolver::set_boundary_values_laplace(): ";
      std::cerr << "\"ERROR: Invalid flag!\"" << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void PlanarPoissonSolver::solve(double* rhs, double* solution) {
  std::vector<double> work(dim);  // work space vector
  double diag_temp;  // temporary variable of the diagonal value
  size_t index;  // index used for going through the lines in reverse
  diag_temp = diag[0];
  if (diag_temp == 0.) {
    std::cerr << "PlanarPoissonSolver::solve(): ";
    std::cerr << "\"ERROR: Reduce your matrix by one line!\"" << std::endl;
    exit(1);
  }
  solution[0] = rhs[0] / diag_temp;
  for (size_t i = 1; i < dim; ++i) {
    work.at(i) = upper[i-1] / diag_temp;
    diag_temp = diag[i] - lower[i] * work.at(i);  // subtraction of two lines
    if (diag_temp == 0.) {
      std::cerr << "PlanarPoissonSolver::solve(): \"ERROR: ";
      std::cerr << "Algorithm failed in row " << i << "!\"" << std::endl;
      exit(1);
    }
    solution[i] = (rhs[i] - lower[i] * solution[i-1]) / diag_temp;
  }
  // The matrix equation is solved in the last line, because there are only two
  // unknowns left instead of three. The solution is obtained by backward
  // substitution.
  for (size_t i = 0; i < dim-1 ; ++i) {
    index = dim - 2 - i;  // go through the arrays backwards
    solution[index] += -work.at(index+1) * solution[index+1];
  }
}
// _____________________________________________________________________________
size_t PlanarPoissonSolver::size() {
  return dim;
}
// _____________________________________________________________________________
