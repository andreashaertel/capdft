// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file src/convergence_criterion.cpp
 *  \brief This file contains the definitions of the ConvergenceCriterion class
 *
 */
#include "convergence_criterion.hpp"  // NOLINT
// _____________________________________________________________________________
ConvergenceCriterion::ConvergenceCriterion(
    const std::vector<DataFrame<1, double>>& old_profile, 
    const std::vector<DataFrame<1, double>>& new_profile)
  : old_profile(old_profile),
    new_profile(new_profile) {
}
// _____________________________________________________________________________
