// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_CONVERGENCE_CRITERION_NAN_HPP_
#define SRC_CONVERGENCE_CRITERION_NAN_HPP_
/** \file src/convergence_criterion_nan.hpp
 *  \brief This file contains the declarations of the ConvergenceCriterionNan
 *         class.
 */
#include "convergence_criterion.hpp"  // NOLINT
#include <vector>
#include <string>
#include "data_frame.hpp"  // NOLINT
/** \brief This class defines a convergence criterion based on the number of
 *         nans in the density profiles.
 *
 *  This is usually used to interrupt the iteration process if the iterator
 *  "exploded".
 */
class ConvergenceCriterionNan : public ConvergenceCriterion {
 public:
  /** \brief Constructor
   */
  ConvergenceCriterionNan(
      const std::vector<DataFrame<1, double>>& old_profile,
      const std::vector<DataFrame<1, double>>& new_profile,
      const int& nan_threshold);
  /** \brief Destructor
   */
  virtual ~ConvergenceCriterionNan();
  /** \brief Check if the criterion is fulfilled (fulfilled = true)
   *
   *  Checks for nans.
   *
   *  \param Uses a pointer to return the nan_count.
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
#endif  // SRC_CONVERGENCE_CRITERION_NAN_HPP_
