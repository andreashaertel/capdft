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
      const std::vector<Properties>& species_properties,
      double** density_porfile);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal array.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   *  check_weighted_densities(),
   *  calc_partial_derivatives()
   */
  void calc_derivative();
  /** Calculate bulk derivatives
   *
   */
  void calc_bulk_derivative();
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the energy value of this functional, which approaches the excess
   *  free energy contribution of the hard sphere interactions.
   *
   *  \return Returns the functional energy value
   */
  double calc_energy();
  /** \brief Calculate the ideal energy (not part of this functional)
   *
   *  Calculates the ideal free energy of all species that interact via this
   *  functional. The ideal free energy ist not really part of this excess
   *  functional.
   *
   *  \return Returns the ideal free energy value
   */
  double calc_ideal_free_energy();
  /** \brief Calculate the grand potential (not part of this functional)
   *
   * Calculates the grand potential given a certain chemical potential and
   * external potential.
   *
   * \param[in] chemical_potential The chemical potentials of the species.
   * \param[in] external_potential The external potential of the system.
   * \return Returns the grand potential energy
   */
  double calc_grand_potential(
      std::vector<double> chemical_potential,
      std::vector<std::vector<double>> external_potential);

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
  /** \brief Pointer to density profiles
   *
   */
  double** density_profile;
  /** \brief Functional derivatives
   *
   */
  double** functional_derivative;
  /** \brief Calculate the weighted densities
   *
   */
  void calc_weighted_densities();
  /** \brief Check if unphysical values appear in the weighted densities
   *
   */
  void check_weighted_densities();
  /** \brief Calculate the partial derivatives of the free energy densities
   *
   */
   void calc_partial_derivatives();

 protected:
};
#endif  // SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
