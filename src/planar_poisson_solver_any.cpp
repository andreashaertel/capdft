// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file planar_poisson_solver_any.cpp
 *  \brief Source file for the PlanarPoissonSolverAny class.
 *
 *  The file contains the class definitions of the PlanarPoissonSolverAny class.
 */
#include "planar_poisson_solver_any.hpp"
#include <algorithm>
#include <iostream>
// _____________________________________________________________________________
PlanarPoissonSolverAny::PlanarPoissonSolverAny() {
  //
}
// _____________________________________________________________________________
PlanarPoissonSolverAny::PlanarPoissonSolverAny(
    size_t size, double dz, std::vector<size_t> boundary_positions)
  : BandedMatrix::BandedMatrix(size - 2, 2, 2),
    boundary_positions(boundary_positions),
    dz(dz), dzdz(dz * dz) {
  std::sort(boundary_positions.begin(), boundary_positions.end());
  this->boundary_positions.at(1) += 1;
  set_laplacian();
  BandedMatrix::lu_decomposition();
}
// _____________________________________________________________________________
PlanarPoissonSolverAny::~PlanarPoissonSolverAny() {
}
// _____________________________________________________________________________
void PlanarPoissonSolverAny::set_laplacian() {
  size_t pos1{boundary_positions.at(0)};
  size_t pos2{boundary_positions.at(1)};
  for (size_t i = 0; i < pos1; ++i) {
    BandedMatrix::at(2, i) = 1. / dzdz;
    BandedMatrix::at(3, i) = -2. / dzdz;
    BandedMatrix::at(4, i) = 1. / dzdz;
  }
  for (size_t i = pos1; i < pos2; ++i) {
    BandedMatrix::at(1, i) = 1. / dzdz;
    BandedMatrix::at(2, i) = -2. / dzdz;
    BandedMatrix::at(3, i) = 1. / dzdz;
  }
  for (size_t i = pos2; i < BandedMatrix::size(); ++i) {
    BandedMatrix::at(0, i) = 1. / dzdz;
    BandedMatrix::at(1, i) = -2. / dzdz;
    BandedMatrix::at(2, i) = 1. / dzdz;
  }
  // First boundary position
  BandedMatrix::at(2, pos1 - 3) = 12. / (15. * dzdz);
  BandedMatrix::at(3, pos1 - 3) = -20. / (15. * dzdz);
  BandedMatrix::at(4, pos1 - 3) = 0.;
  BandedMatrix::at(2, pos1 - 2) = 2. / (3. * dzdz);
  BandedMatrix::at(3, pos1 - 2) = 6. / (3. * dzdz);
  BandedMatrix::at(4, pos1 - 2) = 0.;
  BandedMatrix::at(2, pos1 - 1) = -12. / (3. * dzdz);
  BandedMatrix::at(3, pos1 - 1) = 4. / (3. * dzdz);
  BandedMatrix::at(4, pos1 - 1) = 0.;
  // Second boundary position
  BandedMatrix::at(1, pos2 - 3) = 4. / (3. * dzdz);
  BandedMatrix::at(2, pos2 - 3) = -12. / (3. * dzdz);
  BandedMatrix::at(3, pos2 - 3) = 0.;
  BandedMatrix::at(1, pos2 - 2) = 6. / (3. * dzdz);
  BandedMatrix::at(2, pos2 - 2) = 2. / (3. * dzdz);
  BandedMatrix::at(3, pos2 - 2) = 0.;
  BandedMatrix::at(1, pos2 - 1) = -20. / (15. * dzdz);
  BandedMatrix::at(2, pos2 - 1) = 12. / (15. * dzdz);
  BandedMatrix::at(3, pos2 - 1) = 0.;
}
// _____________________________________________________________________________
void PlanarPoissonSolverAny::solve(
    double* rhs, std::vector<double> boundary_values, double* solution) {
  add_boundary_values(rhs, boundary_values);
  BandedMatrix::lu_solve(rhs + 1, solution + 2);
  remove_boundary_values(rhs, boundary_values);
  rearrange_solution(solution, boundary_values);
}
// _____________________________________________________________________________
void PlanarPoissonSolverAny::add_boundary_values(
    double* rhs, std::vector<double> boundary_values) {
  size_t pos1{boundary_positions.at(0)};
  size_t pos2{boundary_positions.at(1)};
  rhs[pos1 - 2] += -8. * boundary_values.at(0) / (15. * dzdz);
  rhs[pos1 - 1] += 8. * boundary_values.at(0) / (3. * dzdz);
  rhs[pos1 - 0] += -8. * boundary_values.at(0) / (3. * dzdz);
  rhs[pos2 - 2] += -8. * boundary_values.at(1) / (3. * dzdz);
  rhs[pos2 - 1] += 8. * boundary_values.at(1) / (3. * dzdz);
  rhs[pos2 - 0] += -8. * boundary_values.at(1) / (15. * dzdz);
}
// _____________________________________________________________________________
void PlanarPoissonSolverAny::remove_boundary_values(
    double* rhs, std::vector<double> boundary_values) {
  size_t pos1{boundary_positions.at(0)};
  size_t pos2{boundary_positions.at(1)};
  rhs[pos1 - 2] -= -8. * boundary_values.at(0) / (15. * dzdz);
  rhs[pos1 - 1] -= 8. * boundary_values.at(0) / (3. * dzdz);
  rhs[pos1 - 0] -= -8. * boundary_values.at(0) / (3. * dzdz);
  rhs[pos2 - 2] -= -8. * boundary_values.at(1) / (3. * dzdz);
  rhs[pos2 - 1] -= 8. * boundary_values.at(1) / (3. * dzdz);
  rhs[pos2 - 0] -= -8. * boundary_values.at(1) / (15. * dzdz);
}
// _____________________________________________________________________________
void PlanarPoissonSolverAny::rearrange_solution(
    double* solution, std::vector<double> boundary_values) {
  // Shift everything by two until first boundary is reached
  size_t pos1{boundary_positions.at(0)};
  size_t pos2{boundary_positions.at(1)};
  // Shifting
  for (size_t i = 0; i < pos1 - 1; ++i) {
    solution[i] = solution[i + 2];
  }
  // Missing value at the boundary is computed via linear extrapolation
  solution[pos1 - 1] = 
      (1. * solution[pos1] +
      8. * boundary_values.at(0) -
      3. * solution[pos1 + 1]) / 6.;
  // Continue shifting
  for (size_t i = pos1; i < pos2 - 1; ++i) {
    solution[i] = solution[i + 1];
  }
  // Missing value at the boundary is computed via linear extrapolation
  solution[pos2 - 1] =
      (-3. * solution[pos2 - 1] +
      8. * boundary_values.at(1) +
      1. * solution[pos2]) / 6.;
}
// _____________________________________________________________________________
