// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
/** \file functional_fmt_spherical.hpp
 *  \brief Header file for the FunctionalFMTSpherical class.
 *
 *  The file contains the class declarations of the FunctionalFMTSpherical
 *  class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <vector>
#include "src/properties.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief FunctionalFMTSpherical calculates the FMT functional in spherical geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *
 */
class FunctionalFMTSpherical {
 public:
  /** \brief Constructor
   *
   */
  FunctionalFMTSpherical(
      const Properties& system_properties,
      const std::vector<Properties>& species_properties);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTSpherical();

 private:
  /** \brief System length (radius of the sperical geometry)
   *
   */
  double length;
  /** \brief Number of grid points
   *
   */
  size_t grid_count;
  /** \brief Number of species
   *
   */
  size_t species_count;
  /** \brief Hard sphere diameters
   *
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   *
   */
  std::vector<double> bulk_densities;

 protected:
};
#endif  // SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
