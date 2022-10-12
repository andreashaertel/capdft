// SPDX-FileCopyrightText: 2020 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "integration.hpp"
#include <cmath>
#include "data_frame.hpp"
// _____________________________________________________________________________
double integration_1d_radial_open_closed(
    const DataFrame<1, double>& integrand, double bin_size) {
  double integral{0.};
  double r{0.};
  // Trapezoidal rule
  for (size_t i = 0; i != integrand.size(); ++i) {
    r = bin_size * static_cast<double>(i);
    if (i == 0) {  // open
      integral += 1.5 * integrand.element(i) * r * r;
    } else if (i == integrand.size()-1) {  // closed
      integral += .5 * integrand.element(i) * r * r;
    } else {
      integral += 1. * integrand.element(i) * r * r;
    }
  }
  integral *= bin_size;
  integral *= 4. * M_PI;
  return integral;
}
// _____________________________________________________________________________
double integration_1d_closed(
    const DataFrame<1, double>& integrand, double bin_size) {
  double integral{0.};
  // Trapezoidal rule
  for (size_t i = 0; i != integrand.size(); ++i) {
    if (i == 0 || i == integrand.size()) {
      integral += .5 * integrand.element(i);
    } else {
      integral += 1. * integrand.element(i);
    }
  }
  integral *= bin_size;
  return integral;
}
// _____________________________________________________________________________
double integration_2d_closed(
    const DataFrame<2, double>& integrand, std::vector<double> bin_sizes) {
  double integral{0.};
  std::vector<size_t> grid_counts{integrand.size_dim()};
  DataFrame<1, double> integrand_1d(grid_counts.at(1));
  DataFrame<1, double> integrated_once(grid_counts.at(0));
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      integrand_1d.at(j) = integrand.element(i, j);
    }
    integrated_once.at(i) = 
        integration_1d_closed(integrand_1d, bin_sizes.at(1));
  }
  integral = integration_1d_closed(integrated_once, bin_sizes.at(0));
  return integral;
}
// _____________________________________________________________________________
double integration_3d_closed(
    const DataFrame<3, double>& integrand, std::vector<double> bin_sizes) {
  double integral{0.};
  std::vector<size_t> grid_counts{integrand.size_dim()};
  std::vector<size_t> grid_counts_reduced{integrand.size_dim()};
  std::vector<double> bin_sizes_reduced{bin_sizes};
  grid_counts_reduced.pop_back();
  bin_sizes_reduced.pop_back();
  DataFrame<1, double> integrand_1d(grid_counts.at(2));
  DataFrame<2, double> integrated_once(grid_counts_reduced);
  for (size_t i = 0; i != grid_counts.at(0); ++i) {
    for (size_t j = 0; j != grid_counts.at(1); ++j) {
      for (size_t k = 0; k != grid_counts.at(2); ++k) {
        integrand_1d.at(k) = integrand.element(i, j, k);
      }
      integrated_once.at(i, j) =
          integration_1d_closed(integrand_1d, bin_sizes.at(2));
    }
  }
  integral = integration_2d_closed(integrated_once, bin_sizes_reduced);
  return integral;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
