// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_FMT_PLANAR_HPP_
#define SRC_FUNCTIONAL_FMT_PLANAR_HPP_
/** \file functional_fmt_planar.hpp
 *  \brief Header file for the FunctionalFMTPlanar class.
 *
 *  The file contains the class declarations of the FunctionalFMTPlanar
 *  class.
 *
 */
// _____________________________________________________________________________
// Includes
#include <fftw3.h>
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "functional.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
// Class forward declarations
// _____________________________________________________________________________
/** \brief FunctionalFMTPlanar calculates the FMT functional in the planar
 *  geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *  It is derived from the abstract Functional base class.
 *
 */
class FunctionalFMTPlanar : public Functional {
 public:
  /** \brief Standard Constructor
   *
   */
  FunctionalFMTPlanar();
  /** \brief Manual Constructor
   *
   * This constructor chooses the particle species that are supplied from the
   * affected_species vector.
   *
   */
  FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      const std::vector<size_t>& affected_species);
  /** \brief Automated Constructor
   *
   * This constructor chooses the particle species that are interacting via this
   * functional by checking if the hard sphere diameter exists.
   *
   */
  FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTPlanar();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  supplied DataFrame.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   *  check_weighted_densities(),
   *  calc_partial_derivatives() (which uses calc_local_partial_derivatives())
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  void calc_derivative_warnings(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the (excess free) energy value of this functional
   *
   *  Calculate the energy value of this functional, which corresponds to the
   *  excess free energy contribution of the hard sphere interactions in
   *  equilibrium.
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
  /** \brief Bin sizes in real and Fourier space
   *
   */
  double dz;
  double dkz;
  /** \brief Number of species
   *
   */
  size_t species_count;
  /** \brief Vector that remembers the species, that are affected by this
   *  functional
   */
  std::vector<size_t> affected_species;
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
  const std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief From the system object extract the system properties
   *
   */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From the system object extract the species properties
   *
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Initialize all data frame vectors
   *
   */
  void initialize_all_data_frames();
  /** \brief Calculate the weight functions in Fourier space
   *
   */
  void calc_weights();

 protected:
};
#endif  // SRC_FUNCTIONAL_FMT_PLANAR_HPP_
