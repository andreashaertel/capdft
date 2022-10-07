// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
#define SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
/** \file functional_fmt_cartesian.hpp
 *  \brief Header file for the FunctionalFMTCartesian class.
 *
 *  The file contains the class declarations of the FunctionalFMTCartesian
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
/** \brief FunctionalFMTCartesian calculates the FMT functional in the 3D
 *         cartesian geometry
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values.
 *  It is derived from the abstract Functional base class.
 *
 */
class FunctionalFMTCartesian : public Functional {
 public:
  /** \brief Standard Constructor
   *
   */
  FunctionalFMTCartesian();
  /** \brief Manual Constructor
   *
   * This constructor chooses the particle species that are supplied from the
   * affected_species vector.
   *
   */
  FunctionalFMTCartesian(
      const std::vector<DataFrame<3, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      const std::vector<size_t>& affected_species);
  /** \brief Automated Constructor
   *
   * This constructor chooses the particle species that are interacting via this
   * functional by checking if the hard sphere diameter exists.
   *
   */
  FunctionalFMTCartesian(
      const std::vector<DataFrame<3, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   *
   */
  ~FunctionalFMTCartesian();
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
  std::vector<double> lengths;
  /** \brief Number of grid points in real and Fourier space
   *
   * voxel_count is the total ammount of voxels in the system, while
   * grid_counts holds the number of grid points in every dimension.
   */
  size_t voxel_count;
  std::vector<size_t> grid_counts;
  std::vector<size_t> grid_counts_four;
  /** \brief Bin sizes in real and Fourier space
   *
   */
  std::vector<double> bin_sizes;
  std::vector<double> bin_sizes_four;
  /** \brief Number of species
   */
  size_t species_count;
  /** \brief Vector that remembers the species, that are affected by this
   *  functional
   */
  std::vector<size_t> affected_species;
  /** \brief Hard sphere diameters
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   */
  const std::vector<DataFrame<3, double>>* density_profiles_pointer;
  /** \brief Internal density profiles
   *  
   *  The fftw package resets the arrays it is supposed to transform, thus
   *  we need an internal copy of the density profiles.
   */
  std::vector<DataFrame<3, double>> density_profiles;
  /** \brief Fourier transformed density profiles
   */
  std::vector<DataFrame<3, fftw_complex>> density_profiles_four;
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
  std::vector<DataFrame<3, double>> scalar_weighted_dens_real;
  std::vector<DataFrame<3, double>> vector_weighted_dens_real;
  std::vector<DataFrame<3, double>> tensor_weighted_dens_real;
  std::vector<DataFrame<3, fftw_complex>> scalar_weighted_dens_four;
  std::vector<DataFrame<3, fftw_complex>> vector_weighted_dens_four;
  std::vector<DataFrame<3, fftw_complex>> tensor_weighted_dens_four;
  /** \brief Weight functions
   *
   *  Every vector element contains all weight functions of another species.
   *
   */
  std::vector<std::vector<DataFrame<3, fftw_complex>>> scalar_weights_four;
  std::vector<std::vector<DataFrame<3, fftw_complex>>> vector_weights_four;
  std::vector<std::vector<DataFrame<3, fftw_complex>>> tensor_weights_four;
  /** \brief Partial derivatives of the free energy density w.r.t. the
   *  weighted densities
   *
   * There are three partial derivative types: scalar, vectorial, tensorial.
   * There are four arrays containing the four scalar partial derivatives.
   * There are two arrays containing the z-components of the two vectorial
   * partial derivatives.
   * There are two arrays containing the first (= second) and third element of
   * the tensorial partial derivative.
   *
   */
  std::vector<DataFrame<1, double>> scalar_partial_derivative_real;
  std::vector<DataFrame<1, double>> vector_partial_derivative_real;
  std::vector<DataFrame<1, double>> tensor_partial_derivative_real;
  std::vector<DataFrame<1, fftw_complex>> scalar_partial_derivative_four;
  std::vector<DataFrame<1, fftw_complex>> vector_partial_derivative_four;
  std::vector<DataFrame<1, fftw_complex>> tensor_partial_derivative_four;
  /** Fourier transform of the functional derivative
   */
  std::vector<DataFrame<1, fftw_complex>> functional_derivative_four;
  /** \brief Flags for the Fourier transforms
   *
   *  The first one preserves the input, while the second one might destroy it.
   */
  static const unsigned int flags_keep = FFTW_MEASURE | FFTW_PRESERVE_INPUT;
  static const unsigned int flags_destroy = FFTW_MEASURE;
  /** \brief Calculate the weighted densities
   */
  void calc_weighted_densities();
  /** \brief Check if unphysical values appear in the weighted densities
   */
  void check_weighted_densities();
  /** \brief Calculate the partial derivatives of the free energy densities
   */
  void calc_partial_derivatives();
  /** \brief Calculate the partial derivatives of the free energy densities at
   *         one position
   */
  void calc_local_partial_derivatives(size_t i);
  /** \brief Calculate the weighted partial derivatives of the free energy
   *  densities
   */
  void calc_weighted_partial_derivatives(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate the energy density
   *  
   *  \param The index of the position of which the energy density value is
   *  sought.
   *
   *  \return Returns the functional energy density at the specified position.
   */
  double calc_local_energy_density(size_t position);
  /** \brief From the system object extract the system properties
   */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From the system object extract the species properties
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Initialize all data frame vectors
   */
  void initialize_all_data_frames();
  /** \brief Updates the internal density profile arrays
   *         (before calculating something)
   */
  void update_density_profiles();
  /** \brief Calculate the weight functions in Fourier space
   */
  void calc_weights();
  /** \brief Set all weights to zero
   */
  void set_weights_to_zero();
  /** \brief Radial integration
   */
  double integration(double* data, int n, double delta);

 protected:
};
#endif  // SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
