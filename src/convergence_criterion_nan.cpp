// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file src/convergence_criterion_nan.cpp
 *  \brief This file contains the definitions of the ConvergenceCriterionNan
 *         class.
 */
#include "convergence_criterion_nan.hpp"  // NOLINT
#include <cmath>
#include <bitset>
#include <iostream>
#include <string>
// _____________________________________________________________________________
ConvergenceCriterionNan::ConvergenceCriterionNan(
    const std::vector<DataFrame<1, double>>& old_profile,
    const std::vector<DataFrame<1, double>>& new_profile,
    const int& nan_threshold)
  : ConvergenceCriterion(old_profile, new_profile, nan_threshold) {
}
// _____________________________________________________________________________
ConvergenceCriterionNan::~ConvergenceCriterionNan() {
}
// _____________________________________________________________________________
bool ConvergenceCriterionNan::check(
    double* nan_count) {
  int nan_counter{0};
  // Calculate maximum deviations per species
  for (auto& species : new_profile) {
    for (size_t i = 0; i < species.size(); ++i) {
      if (std::isnan(species.element(i))) {
        ++nan_counter;
      }
      if (nan_counter > threshold_int) {
        *nan_count = static_cast<double>(nan_counter);
        return 1;
      }
    }
  }
  *nan_count = static_cast<double>(nan_counter);
  return 0;
}
// _____________________________________________________________________________
std::string ConvergenceCriterionNan::name() {
  return "Nans";
}
// _____________________________________________________________________________
