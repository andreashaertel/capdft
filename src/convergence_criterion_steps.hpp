// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_STEPS_HPP_
#define SRC_CONVERGENCE_CRITERION_STEPS_HPP_
/** \file src/convergence_criterion_steps.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterion class
 *
 */
#include "convergence_criterion.hpp"
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines convergence criteria based on the old and the new
 *         density profiles.
 *
 */
class ConvergenceCriterionSteps : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   *
   */
  ConvergenceCriterionSteps(
      const std::vector<DataFrame<1, double>>& old_profile, 
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& steps);
  /** \brief Destructor
   *
   */
  virtual ~ConvergenceCriterionSteps() {};
  /** \brief Check if the step criterion is fulfilled (fulfilled = true)
   *
   *  Keeps track of the number of steps (number of time check() was called).
   *
   *  \param Uses a pointer to return the number of remaining steps
   *
   */
  virtual bool check(double* current_deviation);
  /** \brief Return name of the criterion
   *
   */
  virtual std::string name();
  // The following section is found in the base class
  ///** \brief Pointer to the old density profile
  // *
  // */
  //const std::vector<DataFrame<1, double>>& old_profile;
  ///** \brief Pointer to the proposed density profile
  // *
  // */
  //const std::vector<DataFrame<1, double>>& new_profile;
  ///** \brief Threshold value used in the check function
  // *
  // */
  //double threshold;
  ///** \brief Threshold value used in the check function (integer)
  // *
  // */
  //int threshold_int;

 private:
  /** \brief Number of steps, i.e. number of times check()
   *
   */
  int step_count;
};
#endif  // SRC_CONVERGENCE_CRITERION_STEPS_HPP_
