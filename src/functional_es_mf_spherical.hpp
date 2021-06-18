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
#include <vector>
#include "functional.hpp"
#include "properties.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief This class calculates the elctrostatic mean field functional.
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values for point charges in the mean field approximation.
 *
 */
class FunctionalESMFSpherical : public Functional {
 public:
  /** \brief Constructor
   *
   */
  FunctionalESMFSpherical(
      const Properties& system_properties,
      const std::vector<Properties>& species_properties,
      double** density_profile);
  /** \brief Destructor
   *
   */
  ~FunctionalESMFSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal array.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   */
  virtual void calc_derivative();
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative();
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the energy value of this functional, which approaches the excess
   *  free energy contribution of the hard sphere interactions.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

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
  /** \brief Bjerrum length
   *
   */
  double bjerrum_length;
  /** \brief System temperature
   *
   */
  double temperature;
  /** \brief Dielectric constant
   *
   */
  double dielectric;
  /** \brief Charge of the center sphere
   *
   */
  double center_charge;
  /** \brief Hard sphere diameters
   *
   */
  std::vector<double> valencies;
  /** \brief Hard sphere diameters
   *
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   *
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   *
   */
  double** density_profile;
  /** \brief Charge density profiles
   *
   */
  double** charge_density_profile;
  /** \brief Species-dependent electrostatic potentials
   *
   */
  double** species_potentials;
  /** \brief Functional derivatives
   *
   */
  double** functional_derivative;
  /** \brief Calculate the weighted densities
   *
   */
  void calc_weighted_densities();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_MF_SPHERICAL_HPP_
