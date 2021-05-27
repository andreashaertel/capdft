// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_TEMPLATE_DATA_HPP_
#define SRC_TEMPLATE_DATA_HPP_
/** \file template_data.hpp
 *  \brief Header file for the templated Data class.
 *
 *  The file contains the class declarations of the templated Data class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <typeinfo>
#include "data.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief Data class is a universal type class
 *
 */
template <typename T> class TemplateData : public Data {
 public:
  /** \brief Constructor
   *
   */
  explicit TemplateData(T value) {
    type = &(typeid(T));
    this->value = value;
  }
  /** \brief Stores the value of template type
   *
   */
  T value;
};
#endif  // SRC_TEMPLATE_DATA_HPP_
