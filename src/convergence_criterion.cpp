// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file src/convergence_criterion.cpp
 *  \brief This file contains the definitions of the ConvergenceCriterion class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <string>
// _____________________________________________________________________________
ConvergenceCriterion::ConvergenceCriterion(
    const std::vector<DataFrame<1, double>>& old_profile,
    const std::vector<DataFrame<1, double>>& new_profile,
    const double& threshold)
  : old_profile(old_profile),
    new_profile(new_profile),
    threshold(threshold),
    threshold_int(0.) {
}
// _____________________________________________________________________________
ConvergenceCriterion::ConvergenceCriterion(
    const std::vector<DataFrame<1, double>>& old_profile,
    const std::vector<DataFrame<1, double>>& new_profile,
    const int& threshold_int)
  : old_profile(old_profile),
    new_profile(new_profile),
    threshold(0.),
    threshold_int(threshold_int) {
}
// _____________________________________________________________________________
