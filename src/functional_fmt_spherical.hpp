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
#include <fftw3.h>
#include "functional.hpp"
#include "properties.hpp"
#include "system.hpp"
// Class forward declarations
// _____________________________________________________________________________
/** \brief FunctionalFMTSpherical calculates the FMT functional in spherical
 *  geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *
 */
class FunctionalFMTSpherical : public Functional {
 public:
  /** \brief Constructors
   *
   */
  FunctionalFMTSpherical();
  explicit FunctionalFMTSpherical(System* system);
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
  virtual void calc_derivative();
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative();
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
  double dr;  // TODO: think about renaming
  double dkr;
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
  DataField<double>* density_profile_pointer;
  /** \brief Density profiles times the radial position
   *
   *  The sine transform of fftw does not require the point r=0, because it is
   *  zero anyways. The cosine transform, however has a non-zero value there.
   *  For this reason the the array is of length grid_count+1 and a sine
   *  transform starts at the second element.
   *
   */
  DataField<double>* density_profile_times_r;
  /** \brief Fourier transformed density profile
   *
   *  Since it is the Sine transformed density_profile_times_r it also contains
   *  arrays of length grid_count+1.
   *
   */
  DataField<double>* density_profile_four;
  /** \brief Functional derivatives
   *
   */
  DataField<double>* functional_derivative;
  /** \brief Weighted densities
   *
   * There are three weighted density types: scalar, vectorial, tensorial.
   * There are four arrays containing the four scalar weighted densities.
   * There are two arrays containing the z-components of the two vectorial
   * weighted densities.
   * There are two arrays containing the first (= second) and third element of
   * the tensorial weighted density.
   *
   */
  DataField<double>* scalar_weighted_dens_real;
  DataField<double>* vector_weighted_dens_real;
  DataField<double>* tensor_weighted_dens_real;
  DataField<double>* scalar_weighted_dens_four;
  DataField<double>* vector_weighted_dens_four;
  DataField<double>* tensor_weighted_dens_four;
  /** \brief Weight functions
   *
   *  In the radially symmetric case we only need a few weight functions.
   *  Every vector element contains all weight function of another species.
   *
   */
  std::vector<DataField<double>> weights_four;
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
  void calc_partial_derivatives();  // TODO(Moritz): referenz rückgabe
  /** \brief Calculate the energy density
   *  
   *  \param The index of the position of which the energy density value is
   *  sought.
   *
   *  \return Returns the functional energy density at the specified position.
   *
   */
   double calc_local_energy_density(size_t position);

 protected:  // TODO(Moritz): swap protected and private
  /** \brief From the system object extract the system properties
   *
   */
  void extract_system_properties(System* sys);
  /** \brief From the system object extract the species properties
   *
   */
  void extract_species_properties(System* sys);
  /** \brief The density profile is changed outside of this functional, thus
   *  we need to update the internal arrays it before calculating something.
   */
  void update_density_times_r();
  /** \brief Initialize the weight functions in Fourier space
   *  For radially symmetric systems only a few weights are actually needed.
   *  The different weighted densities are obtained by introducing the
   *  prefactors in the convolution.
   *
   */
  void initialize_weights();
  /** \brief Radial integration  // TODO(Moritz): toolbox?
   *
   */
  double radial_integration(double* data, int n, double delta);
  /** \brief Flags for the sine and cosine transforms
   *
   *  The first one preserves the input, while the second one might destroy it.
   *
   */
  static const unsigned int flags_keep = FFTW_MEASURE | FFTW_PRESERVE_INPUT;
  static const unsigned int flags_destroy = FFTW_MEASURE; 
};
#endif  // SRC_FUNCTIONAL_FMT_SPHERICAL_HPP_
