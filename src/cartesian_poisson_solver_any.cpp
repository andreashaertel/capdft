// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file cartesian_poisson_solver_any.cpp
 *  \brief Source file for the CartesianPoissonSolverAny class.
 *
 *  The file contains the definitions of the CartesianPoissonSolverAny class.
 */
#include "cartesian_poisson_solver_any.hpp"  // NOLINT
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include "stl_algorithms.hpp"
// _____________________________________________________________________________
CartesianPoissonSolverAny::CartesianPoissonSolverAny() {
  //
}
// _____________________________________________________________________________
CartesianPoissonSolverAny::CartesianPoissonSolverAny(
    std::vector<size_t> bin_count,
    std::vector<double> bin_size,
    std::vector<bool> periodic_boundaries,
    std::vector<std::vector<std::vector<std::vector<size_t>>>>& boundary_positions,  // NOLINT
    std::vector<std::vector<std::vector<double>>>& boundary_values)
  : SparseMatrix(bin_count.at(0) * bin_count.at(1) * bin_count.at(2)),
    bin_count(bin_count),
    bin_size(bin_size),
    periodic_boundaries(periodic_boundaries) {
  // Calculate the square of the bin sizes
  bin_size_squared.resize(bin_size.size());
  for (size_t i = 0; i < bin_size.size(); ++i) {
    bin_size_squared.at(i) = bin_size.at(i) * bin_size.at(i);
  }
  set_laplacian();
  set_boundary_values(boundary_positions, boundary_values);
}
// _____________________________________________________________________________
CartesianPoissonSolverAny::~CartesianPoissonSolverAny() {
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::solve(
    std::vector<double>& rhs, std::vector<double>& solution) {
  // Termination condition: maximum allowed residual norm
  double accuracy = pow(std::numeric_limits<double>::epsilon(), 1. / 3.);
  // The algorithm becomes inefficient in terms of memory if too many (>50)
  // iterations are performed.
  size_t max_iterations{100};
  size_t iteration_count{0};
  double deviation{0.};
  // Dummy right-hand side with certain elements removed
  add_boundary_values(rhs);
  std::vector<double> rhs_dummy(rhs);
  remove_boundary_values(rhs);
  stl_algorithm::erase(rows_to_remove, rhs_dummy);
  // Dummy solution with reduced size
  std::vector<double> solution_dummy(solution);
  stl_algorithm::erase(rows_to_remove, solution_dummy);
  // Solver loop
  while (!SparseMatrix::solve(
      rhs_dummy, solution_dummy, max_iterations, accuracy, &deviation)) {
    iteration_count += max_iterations;
    std::cout << "CartesianPoissonSolverAny::solve(): \"";
    std::cout << "iterations: " << iteration_count;
    std::cout << "; deviation: " << deviation << " > " << accuracy << "\"";
    std::cout << std::endl << "\033[A\033[K";
  }
  // Add the known potential values at the right place
  for (size_t i = 0; i < cols_to_remove.size(); ++i) {
    solution_dummy.insert(
        solution_dummy.begin() + cols_to_remove.at(i).first,
        cols_to_remove.at(i).second);
  }
  solution = solution_dummy;
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_laplacian() {
  set_bare_laplacian();
  set_boundary_conditions();
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_bare_laplacian() {
  size_t row{0}, col{0};
  double diagonal_value{0.};
  std::vector<double> offdiag_values(0);
  // Calculate all possible matrix values
  diagonal_value = -2. / bin_size_squared.at(0) - 2. / bin_size_squared.at(1) -
      2. / bin_size_squared.at(2);
  for (auto& bss : bin_size_squared) {
    offdiag_values.push_back(1. / bss);
  }
  // Convert the 3D cartesian laplacian into a matrix
  for (size_t i = 0; i < bin_count.at(0); ++i) {
    for (size_t j = 0; j < bin_count.at(1); ++j) {
      for (size_t k = 0; k < bin_count.at(2); ++k) {
        coordinates_to_index(i, j, k, &row);
        SparseMatrix::set(row, row, diagonal_value);
        // The six (spacially) adjacent voxels might contribute to the
        // laplacian. If the voxel does not exist, it will be ignored and
        // handled as boundary condition later.
        if(coordinates_to_index(i - 1, j, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(0));
        }
        if(coordinates_to_index(i + 1, j, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(0));
        }
        if(coordinates_to_index(i, j - 1, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(1));
        }
        if(coordinates_to_index(i, j + 1, k, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(1));
        }
        if(coordinates_to_index(i, j, k - 1, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(2));
        }
        if(coordinates_to_index(i, j, k + 1, &col)) {
          SparseMatrix::set(row, col, offdiag_values.at(2));
        }
      }
    }
  }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_boundary_conditions() {
  size_t row{0}, col{0};
  if (periodic_boundaries.at(0)) {  // if periodic
    for (size_t i = 0; i < bin_count.at(1); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(0, i, j, &row);  // one face of cube
        coordinates_to_index(bin_count.at(0) - 1, i, j, &col);  // opposing face
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(0));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(0));
      }
    }
  } else {  // if not periodic
  }
  if (periodic_boundaries.at(1)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(i, 0, j, &row);
        coordinates_to_index(i, bin_count.at(1) - 1, j, &col);
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(1));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(1));
      }
    }
  } else {
  }
  if (periodic_boundaries.at(2)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(1); ++j) {
        coordinates_to_index(i, j, 0, &row);
        coordinates_to_index(i, j, bin_count.at(2) - 1, &col);
        // PBCs are symmetric:
        // bin on one side is linked to the other side and vice versa
        SparseMatrix::set(row, col, 1. / bin_size_squared.at(2));
        SparseMatrix::set(col, row, 1. / bin_size_squared.at(2));
      }
    }
  } else {
  }
}
// _____________________________________________________________________________
bool CartesianPoissonSolverAny::coordinates_to_index(
    size_t i, size_t j, size_t k, size_t* index) {
  // i, j, k cannot be smaller 0, instead they jump to the largest value.
  if (i >= bin_count.at(0) || j >= bin_count.at(1) || k >= bin_count.at(2)) {
    return false;
  }
  *index = k + bin_count.at(2) * j + bin_count.at(2) * bin_count.at(1) * i;
  return true;
}
// _____________________________________________________________________________
bool CartesianPoissonSolverAny::coordinates_to_index(
    std::vector<size_t> pos, size_t* index) {
  return coordinates_to_index(pos.at(0), pos.at(1), pos.at(2), index);
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::set_boundary_values(
    std::vector<std::vector<std::vector<std::vector<size_t>>>>& boundary_positions,  // NOLINT
    std::vector<std::vector<std::vector<double>>>& boundary_values) {
  std::vector<size_t> position1, position2;
  std::vector<size_t> col_indices;
  std::vector<double> col_values;
  double boundary_value1, boundary_value2;
  size_t index1{0}, index2{0};  // matrix indices
  // Reset global variables
  cols_to_remove.clear();
  rows_to_remove.clear();
  rhs_addition.clear();
  rhs_addition.resize(bin_count.at(0) * bin_count.at(1) * bin_count.at(2), 0.);
  for (size_t dir = 0; dir < boundary_positions.size(); ++dir) {
    for (size_t i = 0; i < boundary_positions.at(dir).size(); ++i) {
      // Convert coordinate pair to indices
      position1 = boundary_positions.at(dir).at(i).at(0);
      position2 = boundary_positions.at(dir).at(i).at(1);
      boundary_value1 = boundary_values.at(dir).at(i).at(0);
      boundary_value2 = boundary_values.at(dir).at(i).at(1);
      coordinates_to_index(position1, &index1);
      coordinates_to_index(position2, &index2);
      cols_to_remove.push_back(std::make_pair(index1, boundary_value1));
      cols_to_remove.push_back(std::make_pair(index2, boundary_value2));
      // The first and last entry of that direction are removed
      if (dir == 0) { 
        coordinates_to_index(0, position1.at(1), position1.at(2), &index1);
        coordinates_to_index(
            bin_count.at(0) - 1, position1.at(1), position1.at(2), &index2);
      } else if (dir == 1) {
        coordinates_to_index(position1.at(0), 0, position1.at(2), &index1);
        coordinates_to_index(
            position1.at(0), bin_count.at(1) - 1, position1.at(2), &index2);
      } else if (dir == 2) {
        coordinates_to_index(position1.at(0), position1.at(1), 0, &index1);
        coordinates_to_index(
            position1.at(0), position1.at(1), bin_count.at(2) - 1, &index2);
      }
      rows_to_remove.push_back(index1);
      rows_to_remove.push_back(index2);
    }
  }
  // The voxels where we know the values are transferred to the rhs and the
  // respective matrix columns are removed (in reverse index order).
  std::sort(cols_to_remove.begin(), cols_to_remove.end());
  std::vector<size_t> cols_to_remove_indices(cols_to_remove.size());
  for (int i = cols_to_remove.size() - 1; i >= 0; --i) {
    // For removal save only the column indices
    cols_to_remove_indices.at(i) = cols_to_remove.at(i).first;
    // Calculate the modified rhs
    col_indices = SparseMatrix::get_col_indices(cols_to_remove.at(i).first);
    col_values = SparseMatrix::get_col_values(cols_to_remove.at(i).first);
    for (size_t j = 0; j < col_indices.size(); ++j) {
      rhs_addition.at(col_indices.at(j)) -=
          col_values.at(j) * cols_to_remove.at(i).second;
    }
  }
  SparseMatrix::remove_columns(cols_to_remove_indices);
  // The first and last position in a direction with boundary condition are now
  // redundant, and are hence removed to keep a square matrix.
  // Rows are removed by removing columns from the transpose.
  std::sort(rows_to_remove.begin(), rows_to_remove.end());
  SparseMatrix::transpose();
  SparseMatrix::remove_columns(rows_to_remove);
  SparseMatrix::transpose();
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::add_boundary_values(
    std::vector<double>& rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) { rhs.at(i) += rhs_addition.at(i); }
}
// _____________________________________________________________________________
void CartesianPoissonSolverAny::remove_boundary_values(
    std::vector<double>& rhs) {
  for (size_t i = 0; i < rhs.size(); ++i) { rhs.at(i) -= rhs_addition.at(i); }
}
// _____________________________________________________________________________
// _____________________________________________________________________________
