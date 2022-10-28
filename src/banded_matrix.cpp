// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file banded_matrix.cpp
 *  \brief Source file for the BandedMatrix class.
 *
 *  The file contains the class definitions of the BandedMatrix class.
 */
#include "banded_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
// _____________________________________________________________________________
BandedMatrix::BandedMatrix() {
  //
}
// _____________________________________________________________________________
BandedMatrix::BandedMatrix(size_t row_count, size_t upper_bw, size_t lower_bw)
  : row_count(row_count),
    upper_bw(upper_bw),
    lower_bw(lower_bw) {
  recreate_banded_matrix();
  recreate_lu_containers();
}
// _____________________________________________________________________________
BandedMatrix::~BandedMatrix() {
}
// _____________________________________________________________________________
void BandedMatrix::zero() {
  for (auto& band : bands) {
    for (auto& element : band) {
      element = 0.;
    }
  }
}
// _____________________________________________________________________________
size_t BandedMatrix::size() {
  return row_count;  // size of diagonal band
}
// _____________________________________________________________________________
double& BandedMatrix::at(size_t band, size_t index) {
  if (is_out_of_bounds(band, index)) {
    std::cerr << "BandedMatrix::at(" << band << ", " << index;
    std::cerr << "): \"ERROR: ";
    std::cerr << "Chosen element is out of bounds!\"";
    std::cerr << std::endl;
    exit(1);
  }
  return bands.at(band).at(index);
}
// _____________________________________________________________________________
double& BandedMatrix::element(size_t row, size_t col) {
  size_t band{0};
  long int diff_int{0};
  long int lower_bw_int{static_cast<long int>(lower_bw)};
  diff_int = static_cast<long int>(col) - static_cast<long int>(row);
  band = static_cast<size_t>(lower_bw_int + diff_int);
  return at(band, row);
}
// _____________________________________________________________________________
void BandedMatrix::lu_decomposition() {
  size_t band_count{lower_bw + upper_bw + 1};
  size_t l{0}, m{0};
  double dummy{0.};
  u_matrix = bands;
  // Move all rows to the left, such that all zeros appear on the right
  l = lower_bw;
  for (size_t i = 0; i < lower_bw; ++i) {
    for (size_t j = lower_bw - i; j < band_count; ++j) {
      u_matrix.at(j - l).at(i) = u_matrix.at(j).at(i);
    }
    --l;
    for (size_t j = band_count - l - 1; j < band_count; ++j) {
      u_matrix.at(j).at(i) = 0.;
    }
  }
  // LU decomposition algorithm
  l = lower_bw;
  for (size_t i = 0; i < row_count; ++i) {
    // Check for pivoting
    dummy = u_matrix.at(0).at(i);
    m = i;
    if (l < row_count) { ++l; }
    for (size_t j = i + 1; j < l; ++j) {
      if (fabs(u_matrix.at(0).at(j)) > fabs(dummy)) {
        dummy = u_matrix.at(0).at(j);
        m = j;
      }
    }
    pivot.at(i) = m + 1;
    if (dummy == 0.) {
      u_matrix.at(0).at(i) = std::numeric_limits<double>::epsilon();
      std::cerr << "BandedMatrix::solve(): \"Warning: ";
      std::cerr << "Banded matrix appears to be singular.\"";
      std::cerr << std::endl;
    }
    // Swap rows
    if (m != i) {
      for (size_t j = 0; j < band_count; ++j) {
        std::swap(u_matrix.at(j).at(i), u_matrix.at(j).at(m));
      }
    }
    // LU decomposition
    for (size_t j = i + 1; j < l; ++j) {
      dummy = u_matrix.at(0).at(j) / u_matrix.at(0).at(i);
      l_matrix.at(j - i - 1).at(i) = dummy;
      for (size_t k = 1; k <  band_count; ++k) {
        u_matrix.at(k - 1).at(j) = 
            u_matrix.at(k).at(j) - dummy * u_matrix.at(k).at(i);
      }
      u_matrix.at(band_count - 1).at(j) = 0.;
    }
  }
}
// _____________________________________________________________________________
void BandedMatrix::lu_solve(double* rhs, double* solution) {
  size_t band_count{lower_bw + upper_bw + 1};
  size_t l{0}, m{0};
  double dummy{0.};
  for (size_t i = 0; i < row_count; ++i) {
    solution[i] = rhs[i];
  }
  l = lower_bw;
  for (size_t i = 0; i < row_count; ++i) {
    m = pivot.at(i) - 1;
    if (i != m) { std::swap(solution[i], solution[m]); }
    if (l < row_count) { ++l; }
    for (size_t j = i + 1; j < l; ++j) {
      solution[j] += -l_matrix.at(j - i - 1).at(i) * solution[i];
    }
  }
  l = 1;
  for (long int i = row_count - 1; i >= 0; --i) {
    dummy = solution[i];
    for (size_t j = 1; j < l; ++j) {
      dummy += -u_matrix.at(j).at(i) * solution[i + j];
    }
    solution[i] = dummy / u_matrix.at(0).at(i);
    if (l < band_count) { ++l; }
  }
}
// _____________________________________________________________________________
bool BandedMatrix::is_out_of_bounds(size_t band, size_t index) {
  bool is_out_of_bounds{false};
  // Check if out of bounds
  if (band < lower_bw) {  // lower
    if (index < (lower_bw - band) || index >= row_count) {
      is_out_of_bounds = true;
    }
  } else if (band >= lower_bw && band <= lower_bw + upper_bw) {  // upper
    if (index >= row_count - (lower_bw - band)) {
      is_out_of_bounds = true;
    }
  } else {
    is_out_of_bounds = true;
  }
  return is_out_of_bounds;
}
// _____________________________________________________________________________
void BandedMatrix::recreate_banded_matrix() {
  bands.clear();
  // Lower bands
  for (size_t i = 0; i < lower_bw; ++i) {
    bands.push_back(std::vector<double>(row_count));
  }
  // Diagonal
  bands.push_back(std::vector<double>(row_count));
  // Upper bands
  for (size_t i = 0; i < upper_bw; ++i) {
    bands.push_back(std::vector<double>(row_count));
  }
  zero();
}
// _____________________________________________________________________________
void BandedMatrix::recreate_lu_containers() {
  u_matrix.clear();
  u_matrix.resize(lower_bw + upper_bw + 1, std::vector<double>(row_count));
  l_matrix.clear();
  l_matrix.resize(lower_bw, std::vector<double>(row_count));
  pivot.clear();
  pivot.resize(row_count);
}
// _____________________________________________________________________________
