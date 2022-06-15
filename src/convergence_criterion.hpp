// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_HPP_
#define SRC_CONVERGENCE_CRITERION_HPP_
/** \file src/convergence_criterion.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterion class
 *
 */
#include <vector>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines convergence criteria based on the old and the new
 *         density profiles.
 *
 */
class ConvergenceCriterion {
 public:
  /** \brief Constructor
   *
   */
  ConvergenceCriterion(const std::vector<DataFrame<1, double>>& old_profile, 
      const std::vector<DataFrame<1, double>>& new_profile);
  /** \brief Destructor
   *
   */
  virtual ~ConvergenceCriterion() {};
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  The implementation of the criterion should go here.
   *
   *  \param Uses a pointer to return the convergence progress
   *         (e.g. the largest deviation between the profiles).
   *
   *  \param The threshold determines, when the convergence criterion is reached
   *         given the old and the new density profiles.
   */
  virtual bool check(double* progress, double threshold) = 0;
  /** \brief Pointer to the old density profile
   *
   */
  const std::vector<DataFrame<1, double>>& old_profile;
  /** \brief Pointer to the proposed density profile
   *
   */
  const std::vector<DataFrame<1, double>>& new_profile;
};
#endif  // SRC_CONVERGENCE_CRITERION_HPP_
