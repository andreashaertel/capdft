// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
/** \file functional_es_mf_spherical.hpp
 *  \brief Header file for the FunctionalESMFSpherical class.
 *
 *  The file contains the class declarations of the FunctionalESMFSpherical
 *  class.
 *
 */
// _____________________________________________________________________________
// Includes
// Class forward declarations
// _____________________________________________________________________________
/** \brief This class calculates the elctrostatic mean field functional.
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values for point charges in the mean field approximation.
 *
 */
class FunctionalESMFSpherical {
 public:
  /** \brief Constructor
   *
   */
  FunctionalESMFSpherical();
  /** \brief Destructor
   *
   */
  ~FunctionalESMFSpherical();
 private:
 protected:
};
#endif  // SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
