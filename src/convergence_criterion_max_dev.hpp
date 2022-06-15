// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
#define SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
/** \file src/convergence_criterion_max_dev.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterion class
 *
 */
#include "convergence_criterion.hpp"
#include <vector>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines convergence criteria based on the old and the new
 *         density profiles.
 *
 */
class ConvergenceCriterionMaxDev : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   *
   */
  ConvergenceCriterionMaxDev(
      const std::vector<DataFrame<1, double>>& old_profile, 
      const std::vector<DataFrame<1, double>>& new_profile);
  /** \brief Destructor
   *
   */
  virtual ~ConvergenceCriterionMaxDev() {};
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  Calculates the largest difference between old and new (proposed)
   *  density profiles and checks if it is smaller than the target_deviation.
   *
   *  \param Uses a pointer to return the current_deviation.
   *
   *  \param The target_deviation determines, when the convergence criterion is
   *         reached.
   */
  virtual bool check(double* current_deviation, double target_deviation);
};
#endif  // SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
