// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_HPP_
#define SRC_CONVERGENCE_CRITERION_HPP_
/** \file src/convergence_criterion.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterion
 *         class.
 */
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines convergence criteria based on the old and the new
 *         density profiles.
 */
class ConvergenceCriterion {
 public:
  /** \brief Constructor for floating point thresholds
   */
  ConvergenceCriterion(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const double& threshold);
  /** \brief Constructor for integer thresholds
   *
   *  This constructor is used, when counting steps etc.
   */
  ConvergenceCriterion(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& threshold_int);
  /** \brief Destructor
   *
   */
  virtual ~ConvergenceCriterion() {}
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  The implementation of the criterion should go here.
   *
   *  \param Uses a pointer to return the convergence progress
   *         (e.g. the largest deviation between the profiles).
   *
   *  \param The threshold determines, when the convergence criterion is reached
   *         given the old and the new density profiles.
   *
   *  \return true if converged, false otherwise
   */
  virtual bool check(double* progress) = 0;
  /** \brief Return name of the criterion
   */
  virtual std::string name() = 0;
  /** \brief Pointer to the old density profile
   */
  const std::vector<DataFrame<1, double>>& old_profile;
  /** \brief Pointer to the proposed density profile
   */
  const std::vector<DataFrame<1, double>>& new_profile;
  /** \brief Threshold value used in the check function (floating point)
   */
  double threshold;
  /** \brief Threshold value used in the check function (integer)
   */
  int threshold_int;
};
#endif  // SRC_CONVERGENCE_CRITERION_HPP_
