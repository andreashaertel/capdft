// SPDX-FileCopyrightText: 2020 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_INTEGRATION_HPP_
#define SRC_INTEGRATION_HPP_
/** \file integration.hpp
 *  \brief Header file for functions that integrate over DataFrames.
 *
 *  The file contains the function declarations for various integration
 *  algorithms with DataFrames as input.
 *
 */
// _____________________________________________________________________________
// Includes
#include "data_frame.hpp"
/** \brief 1D radial integration
 *
 *  The integrand is integrated with the specified bin_sizes over the entire
 *  space using the algorithms found in "Numerical Recipes 3rd edition"
 *  in equations 4.1.3-4.1.18 (the terms open and closed are also explained
 *  there).
 *
 *  Since this is a radial integration the result is multiplied by 4\pi and the
 *  integrand is multiplied by the radial coordinate squared.
 */
double integration_1d_radial_open_closed(
    const DataFrame<1, double>& integrand, double bin_size);
/** \brief 1D regular integration
 *
 *  The integrand is integrated with the specified bin_sizes over the entire
 *  space using the algorithms found in "Numerical Recipes 3rd edition"
 *  in equations 4.1.3-4.1.18 (the terms open and closed are also explained
 *  there).
 */
double integration_1d_closed(
    const DataFrame<1, double>& integrand, double bin_size);
/** \brief 2D regular integration
 *
 *  The integrand is integrated with the specified bin_sizes over the entire
 *  space using the algorithms found in "Numerical Recipes 3rd edition"
 *  in equations 4.1.3-4.1.18 (the terms open and closed are also explained
 *  there).
 *
 *  For the 2D integration the integrand is successively integrated in one
 *  dimension, thus using the 1d algorithm multiple times.
 */
double integration_2d_closed(
    const DataFrame<2, double>& integrand, std::vector<double> bin_sizes);
/** \brief 3D regular integration
 *
 *  The integrand is integrated with the specified bin_sizes over the entire
 *  space using the algorithms found in "Numerical Recipes 3rd edition"
 *  in equations 4.1.3-4.1.18 (the terms open and closed are also explained
 *  there).
 *
 *  For the 3D integration the integrand is successively integrated in one
 *  dimension, thus using the 1d algorithm multiple times.
 */
double integration_3d_closed(
    const DataFrame<3, double>& integrand, std::vector<double> bin_sizes);
#endif  // SRC_INTEGRATION_HPP_
