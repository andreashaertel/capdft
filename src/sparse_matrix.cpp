// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file banded_matrix.cpp
 *  \brief Source file for the SparseMatrix class.
 *
 *  The file contains the class definitions of the SparseMatrix class.
 */
#include "sparse_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
// _____________________________________________________________________________
SparseMatrix::SparseMatrix(size_t row_count, size_t col_count)
  : row_count(row_count),
    col_count(col_count),
    size(0) {
  zero();
}
// _____________________________________________________________________________
SparseMatrix::SparseMatrix(size_t row_col_count)
  : SparseMatrix(row_col_count, row_col_count) {
}
// _____________________________________________________________________________
SparseMatrix::SparseMatrix() : SparseMatrix(0, 0) {
}
// _____________________________________________________________________________
SparseMatrix::~SparseMatrix() {
}
// _____________________________________________________________________________
void SparseMatrix::zero() {
  values.clear();
  row_index.clear();
  values.resize(col_count);
  row_index.resize(col_count);
}
// _____________________________________________________________________________
double SparseMatrix::get(size_t row, size_t col) {
  for (size_t i = 0; i < row_index.at(col).size(); ++i) {
    if (row_index.at(col).at(i) == row) {
      return values.at(col).at(i);
    }
  }
  return 0.;
}
// _____________________________________________________________________________
void SparseMatrix::set(size_t row, size_t col, double value) {
  size_t insert_index{0};
  for (size_t i = 0; i < row_index.at(col).size() + 1; ++i) {
    if (i == row_index.at(col).size()) {  // spot comes after last element
      insert_index = i;
      break;
    } else if (row < row_index.at(col).at(i)) {  // we jumped over the spot
      insert_index = i;
      break;
    } else if (row == row_index.at(col).at(i)) {  // we found the spot
      if (fabs(value) > 0.) { values.at(col).at(i) = value; }
      else { erase(i, col); }
      return;
    }
  }
  if (fabs(value) > 0.) {
    insert(insert_index, row, col, value);
  }
}
// _____________________________________________________________________________
void SparseMatrix::insert(
    size_t position, size_t row, size_t col, double value) {
  values.at(col).insert(values.at(col).begin() + position, value);
  row_index.at(col).insert(row_index.at(col).begin() + position, row);
  ++size;
}
// _____________________________________________________________________________
void SparseMatrix::erase(size_t position, size_t col) {
  values.at(col).erase(values.at(col).begin() + position);
  row_index.at(col).erase(row_index.at(col).begin() + position);
  --size;
}
// _____________________________________________________________________________
std::vector<double> SparseMatrix::matrix_times_vector(std::vector<double>& x) {
  std::vector<double> result(row_count, 0.);
  if (x.size() != col_count) {
    std::cerr << "SparseMatrix::matrix_times_vector(): \"ERROR: ";
    std::cerr << "Matrix and vector have incompatible sizes!\"" << std::endl;
    exit(1);
  }
  for (size_t i = 0; i < row_index.size(); ++i) {  // go through columns
    for (size_t j = 0; j < row_index.at(i).size(); ++j) {  // go through rows
      result.at(row_index.at(i).at(j)) +=
          values.at(i).at(j) * x.at(i);
    }
  }
  return result;
}
// _____________________________________________________________________________
std::vector<double> SparseMatrix::vector_times_matrix(std::vector<double>& x) {
  std::vector<double> result(col_count, 0.);
  if (x.size() != row_count) {
    std::cerr << "SparseMatrix::vector_times_matrix(): \"ERROR: ";
    std::cerr << "Matrix and vector have incompatible sizes!\"" << std::endl;
    exit(1);
  }
  for (size_t i = 0; i < row_index.size(); ++i) {  // go through columns
    for (size_t j = 0; j < row_index.at(i).size(); ++j) {  // go through rows
      result.at(i) += values.at(i).at(j) * x.at(row_index.at(i).at(j));
    }
  }
  return result;
}
// _____________________________________________________________________________
double SparseMatrix::vector_times_vector(
    std::vector<double>& x, std::vector<double>& y) {
  if (x.size() != y.size()) {
    std::cerr << "SparseMatrix::vector_times_vector(): \"ERROR: ";
    std::cerr << "The vectors have incompatible sizes!\"" << std::endl;
    exit(1);
  }
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.);
}
// _____________________________________________________________________________
bool SparseMatrix::solve(
    std::vector<double>& rhs, std::vector<double>& solution,
    size_t max_iterations, double epsilon) {
  bool converged{false};
  size_t j{0};
  double dummy1{0.}, dummy2{0.}, beta{0.};
  std::vector<double> residual(rhs.size(), 0.);
  std::vector<std::vector<double>> h(
      max_iterations + 1, std::vector<double>(max_iterations));
  std::vector<std::vector<double>> v(0);  // normalised residuals
  std::vector<double> temp(rhs.size(), 0.);
  std::vector<double> q(rhs.size(), 0.);
  std::vector<double> w(rhs.size(), 0.);
  std::vector<double> s(1, 0.);
  std::vector<double> c(1, 0.);
  std::vector<double> gamma(max_iterations + 1, 0.);
  std::vector<double> y(max_iterations, 0.);
  // Check if input vectors have correct size
  if (rhs.size() != col_count || solution.size() != col_count) {
    std::cerr << "SparseMatrix::solve(): \"ERROR: \"";
    std::cerr << "Provided vectors have incorrect size!\"" << std::endl;
    exit(1);
  }
  // Warn if the number of iterations is large
  if (max_iterations > 50) {
    std::cerr << "SparseMatrix::solve(): \"Warning: \"";
    std::cerr << "The maximum number of iterations is very large. Consider ";
    std::cerr << "repeating the algorithm multiple times instead.\"";
    std::cerr << std::endl;
  }
  // Calculate initial residual
  temp = matrix_times_vector(solution);
  for (size_t i = 0; i < rhs.size(); ++i) {
    residual.at(i) = rhs.at(i) - temp.at(i);
  }
  gamma.at(0) = sqrt(vector_times_vector(residual, residual));
  for (size_t i = 0; i < residual.size(); ++i) {
    residual.at(i) = residual.at(i) / gamma.at(0);
  }
  v.push_back(residual);
  if (gamma.at(0) < epsilon) {
    converged = true;
    return converged;
  }
  // Iteration loop
  for (j = 0; j < max_iterations; ++j) {
    q = matrix_times_vector(v.at(j));
    for (size_t i = 0; i <= j; ++i) {
      h.at(i).at(j) = vector_times_vector(v.at(i), q);
    }
    w = q;
    for (size_t k = 0; k < w.size(); ++k) {  // vector elements
      for (size_t i = 0; i <= j; ++i) {  // sum
        w.at(k) -= h.at(i).at(j) * v.at(i).at(k);
      }
    }
    h.at(j + 1).at(j) = sqrt(vector_times_vector(w, w));
    for (size_t i = 0; i < j; ++i) {
      dummy1 = c.at(i + 1) * h.at(i).at(j) + s.at(i + 1) * h.at(i + 1).at(j);
      dummy2 = -s.at(i + 1) * h.at(i).at(j) + c.at(i + 1) * h.at(i + 1).at(j);
      h.at(i).at(j) = dummy1;
      h.at(i + 1).at(j) = dummy2;
    }
    beta = sqrt(pow(h.at(j).at(j), 2) + pow(h.at(j + 1).at(j), 2));
    s.push_back(h.at(j + 1).at(j) / beta);
    c.push_back(h.at(j).at(j) / beta);
    h.at(j).at(j) = beta;
    gamma.at(j + 1) = -s.at(j + 1) * gamma.at(j);
    gamma.at(j) *= c.at(j + 1);
    if (fabs(gamma.at(j + 1)) >= epsilon) {
      v.push_back(w);
      for (size_t i = 0; i < v.at(j + 1).size(); ++i) {
        v.at(j + 1).at(i) /= h.at(j + 1).at(j);
      }
    } else {
      ++j;
      converged = true;
      break;
    }
  }
  --j;
  for (long int i = j; i >= 0; --i) {
    y.at(i) = gamma.at(i) / h.at(i).at(i);
    for (size_t k = i + 1; k <= j; ++k) {
      y.at(i) -= h.at(i).at(k) * y.at(k) / h.at(i).at(i);
    }
  }
  for (size_t k = 0; k < solution.size(); ++k) {  // vector elements
    for (size_t i = 0; i <= j; ++i) {  // sum
      solution.at(k) += y.at(i) * v.at(i).at(k);
    }
  }
  return converged;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
