// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
#define SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
/** \file src/convergence_criterion_max_dev.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterionMaxDev
 *         class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines a convergence criterion based on the largest
 *         occuring deviation between old and new density profiles.
 */
class ConvergenceCriterionMaxDev : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   */
  ConvergenceCriterionMaxDev(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const double& target_deviation);
  /** \brief Destructor
   */
  virtual ~ConvergenceCriterionMaxDev();
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  Calculates the largest difference between old and new (proposed)
   *  density profiles and checks if it is smaller than the target_deviation.
   *
   *  \param Uses a pointer to return the current_deviation.
   */
  virtual bool check(double* current_deviation);
  /** \brief Return name of the criterion
   */
  virtual std::string name();
  // The following section is found in the base class
  // /** \brief Pointer to the old density profile
  //  */
  // const std::vector<DataFrame<1, double>>& old_profile;
  // /** \brief Pointer to the proposed density profile
  //  */
  // const std::vector<DataFrame<1, double>>& new_profile;
  // /** \brief Threshold value used in the check function
  //  */
  // double threshold;
  // /** \brief Threshold value used in the check function (integer)
  //  */
  // int threshold_int;
};
#endif  // SRC_CONVERGENCE_CRITERION_MAX_DEV_HPP_
