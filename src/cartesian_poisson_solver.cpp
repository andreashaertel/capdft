// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file cartesian_poisson_solver.cpp
 *  \brief Source file for the CartesianPoissonSolver class.
 *
 *  The file contains the definitions of the CartesianPoissonSolver class.
 */
#include "cartesian_poisson_solver.hpp"  // NOLINT
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
// _____________________________________________________________________________
CartesianPoissonSolver::CartesianPoissonSolver() {
  //
}
// _____________________________________________________________________________
CartesianPoissonSolver::CartesianPoissonSolver(
    std::vector<size_t> bin_count, std::vector<double> bin_size,
    std::vector<bool> periodic_boundaries)
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
}
// _____________________________________________________________________________
CartesianPoissonSolver::~CartesianPoissonSolver() {
}
// _____________________________________________________________________________
void CartesianPoissonSolver::solve(
    std::vector<double>& rhs,
    std::vector<std::vector<double>> boundary_values,
    std::vector<double>& solution) {
  // Termination condition: maximum allowed residual norm
  double accuracy = pow(std::numeric_limits<double>::epsilon(), 1. / 3.);
  // The algorithm becomes inefficient in terms of memory if too many (>50)
  // iterations are performed.
  size_t max_iterations{25};
  size_t iteration_count{0};
  double deviation{0.};
  add_boundary_values(boundary_values, rhs);
  while (!SparseMatrix::solve(
      rhs, solution, max_iterations, accuracy, &deviation)) {
    iteration_count += max_iterations;
    std::cout << "CartesianPoissonSolver::solve(): \"";
    std::cout << "iterations: " << iteration_count;
    std::cout << "; deviation: " << deviation << " > " << accuracy << "\"";
    std::cout << std::endl << "\033[A\033[K";
  }
  remove_boundary_values(boundary_values, rhs);
}
// _____________________________________________________________________________
void CartesianPoissonSolver::set_laplacian() {
  set_bare_laplacian();
  set_boundary_conditions();
}
// _____________________________________________________________________________
void CartesianPoissonSolver::set_bare_laplacian() {
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
void CartesianPoissonSolver::set_boundary_conditions() {
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
    // TODO: set the boundary to half bin size
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
    // TODO: 
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
    // TODO: 
  }
}
// _____________________________________________________________________________
bool CartesianPoissonSolver::coordinates_to_index(
    size_t i, size_t j, size_t k, size_t* index) {
  // i, j, k cannot be smaller 0, instead they jump to the largest value
  if (i >= bin_count.at(0) || j >= bin_count.at(1) || k >= bin_count.at(2)) {
    return false;
  }
  *index = k + bin_count.at(2) * j + bin_count.at(2) * bin_count.at(1) * i;
  return true;
}
// _____________________________________________________________________________
void CartesianPoissonSolver::add_boundary_values(
    std::vector<std::vector<double>>& boundary_values,
    std::vector<double>& rhs) {
  size_t index{0};
  if (!periodic_boundaries.at(0)) {  // if periodic
    for (size_t i = 0; i < bin_count.at(1); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(0, i, j, &index);
        rhs.at(index) += -boundary_values.at(0).at(0) / bin_size_squared.at(0);
        coordinates_to_index(bin_count.at(0) - 1, i, j, &index);
        rhs.at(index) += -boundary_values.at(0).at(1) / bin_size_squared.at(0);
      }
    }
  }
  if (!periodic_boundaries.at(1)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(i, 0, j, &index);
        rhs.at(index) += -boundary_values.at(1).at(0) / bin_size_squared.at(1);
        coordinates_to_index(i, bin_count.at(1) - 1, j, &index);
        rhs.at(index) += -boundary_values.at(1).at(1) / bin_size_squared.at(1);
      }
    }
  }
  if (!periodic_boundaries.at(2)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(1); ++j) {
        coordinates_to_index(i, j, 0, &index);
        rhs.at(index) += -boundary_values.at(2).at(0) / bin_size_squared.at(2);
        coordinates_to_index(i, j, bin_count.at(2) - 1, &index);
        rhs.at(index) += -boundary_values.at(2).at(1) / bin_size_squared.at(2);
      }
    }
  }
}
// _____________________________________________________________________________
void CartesianPoissonSolver::remove_boundary_values(
    std::vector<std::vector<double>>& boundary_values,
    std::vector<double>& rhs) {
  size_t index{0};
  if (!periodic_boundaries.at(0)) {  // if periodic
    for (size_t i = 0; i < bin_count.at(1); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(0, i, j, &index);
        rhs.at(index) -= -boundary_values.at(0).at(0) / bin_size_squared.at(0);
        coordinates_to_index(bin_count.at(0) - 1, i, j, &index);
        rhs.at(index) -= -boundary_values.at(0).at(1) / bin_size_squared.at(0);
      }
    }
  }
  if (!periodic_boundaries.at(1)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(2); ++j) {
        coordinates_to_index(i, 0, j, &index);
        rhs.at(index) -= -boundary_values.at(1).at(0) / bin_size_squared.at(1);
        coordinates_to_index(i, bin_count.at(1) - 1, j, &index);
        rhs.at(index) -= -boundary_values.at(1).at(1) / bin_size_squared.at(1);
      }
    }
  }
  if (!periodic_boundaries.at(2)) {
    for (size_t i = 0; i < bin_count.at(0); ++i) {
      for (size_t j = 0; j < bin_count.at(1); ++j) {
        coordinates_to_index(i, j, 0, &index);
        rhs.at(index) -= -boundary_values.at(2).at(0) / bin_size_squared.at(2);
        coordinates_to_index(i, j, bin_count.at(2) - 1, &index);
        rhs.at(index) -= -boundary_values.at(2).at(1) / bin_size_squared.at(2);
      }
    }
  }
}
// _____________________________________________________________________________
// _____________________________________________________________________________
