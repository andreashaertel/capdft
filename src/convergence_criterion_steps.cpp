// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file src/convergence_criterion_steps.cpp
 *  \brief This file contains the definitions of the ConvergenceCriterionSteps
 *         class.
 */
#include "convergence_criterion_steps.hpp"  // NOLINT
#include <string>
// _____________________________________________________________________________
ConvergenceCriterionSteps::ConvergenceCriterionSteps(
    const std::vector<DataFrame<1, double>>& old_profile,
    const std::vector<DataFrame<1, double>>& new_profile,
    const int& steps)
  : ConvergenceCriterion(old_profile, new_profile, steps) {
}
// _____________________________________________________________________________
ConvergenceCriterionSteps::~ConvergenceCriterionSteps() {
}
// _____________________________________________________________________________
bool ConvergenceCriterionSteps::check(double* max_deviation) {
  *max_deviation = static_cast<double>(threshold_int - step_count);
  ++step_count;
  return (step_count > threshold_int);
}
// _____________________________________________________________________________
std::string ConvergenceCriterionSteps::name() {
  return "step counter";
}
// _____________________________________________________________________________
