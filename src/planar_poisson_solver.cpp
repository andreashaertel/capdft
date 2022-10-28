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
    size_t size, double dz, BoundaryFlag flag)
  : BandedMatrix(size, 1, 1),
    dz(dz), dzdz(dz * dz), flag(flag) {
  set_laplacian();
  BandedMatrix::lu_decomposition();
}
// _____________________________________________________________________________
PlanarPoissonSolver::~PlanarPoissonSolver() {
  //
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
  size_t size{BandedMatrix::size()};
  for (size_t i = 0; i < size; ++i) {
    // Lower band has no first element
    if (i != 0) { BandedMatrix::at(0, i) = 1. / dzdz; }
    // Diagonal band
    BandedMatrix::at(1, i) = -2. / dzdz;
    // Upper band has no last element
    if (i != size - 1) { BandedMatrix::at(2, i) = 1. / dzdz; }
  }
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_NN() {
  size_t size{BandedMatrix::size()};
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  BandedMatrix::at(1, 0) = -1. / dzdz;
  BandedMatrix::at(2, 0) = 1. / dzdz;
  BandedMatrix::at(0, size - 1) = 1. / dzdz;
  BandedMatrix::at(1, size - 1) = -1. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_DD() {
  size_t size{BandedMatrix::size()};
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  BandedMatrix::at(1, 0) = -4. / dzdz;
  BandedMatrix::at(2, 0) = 4. / (3. * dzdz);
  BandedMatrix::at(0, size - 1) = 4. / (3. * dzdz);
  BandedMatrix::at(1, size - 1) = -4. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_ND() {
  size_t size{BandedMatrix::size()};
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  BandedMatrix::at(1, 0) = -1. / dzdz;
  BandedMatrix::at(2, 0) = 1. / dzdz;
  BandedMatrix::at(0, size - 1) = 4. / (3. * dzdz);
  BandedMatrix::at(1, size - 1) = -4. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::set_laplacian_DN() {
  size_t size{BandedMatrix::size()};
  // Set Laplacian without boundary conditions
  set_bare_laplacian();
  // Set the second derivative at the boundaries
  BandedMatrix::at(1, 0) = -4. / dzdz;
  BandedMatrix::at(2, 0) = 4. / (3. * dzdz);
  BandedMatrix::at(0, size - 1) = 1. / dzdz;
  BandedMatrix::at(1, size - 1) = -1. / dzdz;
  // See the "Finite Difference Coefficients Calculator" from
  // [https://web.media.mit.edu/~crtaylor/calculator.html].
}
// _____________________________________________________________________________
void PlanarPoissonSolver::solve(
    double left_boundary_value, double right_boundary_value, double* rhs,
    double* solution) {
  add_boundary_values(rhs, left_boundary_value, right_boundary_value);
  BandedMatrix::lu_solve(rhs, solution);
  remove_boundary_values(rhs, left_boundary_value, right_boundary_value);
}
// _____________________________________________________________________________
void PlanarPoissonSolver::add_boundary_values(
    double* rhs, double left_boundary_value, double right_boundary_value) {
  size_t size{BandedMatrix::size()};
  // Add the boundary values to the first and the last element of the
  // right-hand side.
  switch (flag) {
    case NEUMANN_NEUMANN:
      rhs[0] += left_boundary_value / dz;
      rhs[size-1] += -right_boundary_value / dz;
      break;
    case DIRICHLET_DIRICHLET:
      rhs[0] += -8. * left_boundary_value / (3. * dzdz);
      rhs[size-1] += -8. * right_boundary_value / (3. * dzdz);
      break;
    case NEUMANN_DIRICHLET:
      rhs[0] += left_boundary_value / dz;
      rhs[size-1] += -8. * right_boundary_value / (3. * dzdz);
      break;
    case DIRICHLET_NEUMANN:
      rhs[0] += -8. * left_boundary_value / (3. * dzdz);
      rhs[size-1] += -right_boundary_value / dz;
      break;
    default:
      std::cerr << "PlanarPoissonSolver::add_boundary_values(): ";
      std::cerr << "\"ERROR: Invalid flag!\"" << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
void PlanarPoissonSolver::remove_boundary_values(
    double* rhs, double left_boundary_value, double right_boundary_value) {
  size_t size{BandedMatrix::size()};
  // Add the boundary values to the first and the last element of the
  // right-hand side.
  switch (flag) {
    case NEUMANN_NEUMANN:
      rhs[0] -= left_boundary_value / dz;
      rhs[size-1] -= -right_boundary_value / dz;
      break;
    case DIRICHLET_DIRICHLET:
      rhs[0] -= -8. * left_boundary_value / (3. * dzdz);
      rhs[size-1] -= -8. * right_boundary_value / (3. * dzdz);
      break;
    case NEUMANN_DIRICHLET:
      rhs[0] -= left_boundary_value / dz;
      rhs[size-1] -= -8. * right_boundary_value / (3. * dzdz);
      break;
    case DIRICHLET_NEUMANN:
      rhs[0] -= -8. * left_boundary_value / (3. * dzdz);
      rhs[size-1] -= -right_boundary_value / dz;
      break;
    default:
      std::cerr << "PlanarPoissonSolver::remove_boundary_values(): ";
      std::cerr << "\"ERROR: Invalid flag!\"" << std::endl;
      exit(1);
  }
}
// _____________________________________________________________________________
