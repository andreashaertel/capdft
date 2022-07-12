// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file src/convergence_criterion_max_dev.cpp
 *  \brief This file contains the definitions of the ConvergenceCriterionMaxDev
 *         class.
 */
#include "convergence_criterion_max_dev.hpp"  // NOLINT
#include <algorithm>
#include <string>
// _____________________________________________________________________________
ConvergenceCriterionMaxDev::ConvergenceCriterionMaxDev(
    const std::vector<DataFrame<1, double>>& old_profile,
    const std::vector<DataFrame<1, double>>& new_profile,
    const double& target_deviation)
  : ConvergenceCriterion(old_profile, new_profile, target_deviation) {
}
// _____________________________________________________________________________
ConvergenceCriterionMaxDev::~ConvergenceCriterionMaxDev() {
}
// _____________________________________________________________________________
bool ConvergenceCriterionMaxDev::check(
    double* max_deviation) {
  std::vector<double> max_deviations(0);
  // Calculate maximum deviations per species
  for (size_t i = 0; i < old_profile.size(); ++i) {
    max_deviations.push_back(max(abs(
        static_cast<std::vector<DataFrame<1, double>>>(old_profile).at(i) -
        static_cast<std::vector<DataFrame<1, double>>>(new_profile).at(i))));
  }
  // Find the largest deviation among species
  *max_deviation = *std::max_element(
      max_deviations.begin(), max_deviations.end());
  return (*max_deviation < threshold);
}
// _____________________________________________________________________________
std::string ConvergenceCriterionMaxDev::name() {
  return "maximum deviation";
}
// _____________________________________________________________________________
