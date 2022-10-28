// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_DELTA_PLANAR_HPP_
#define SRC_FUNCTIONAL_ES_DELTA_PLANAR_HPP_
/** \file functional_es_delta_planar.hpp
 *  \brief Header file for the FunctionalESDeltaPlanar class.
 *
 *  The file contains the class declarations of the FunctionalESDeltaPlanar
 *  class.
 */
// Includes
#include "functional.hpp"  // NOLINT
#include <fftw3.h>
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
#include "planar_poisson_solver.hpp"  // NOLINT
/** \brief This class calculates the elctrostatic delta functional in the
 *  planar geometry.
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values of point charges using the mean field approximation in the
 *  planar geometry.
 *
 *  The theory for this functional was published in
 *  [M. Bültmann and A. Härtel 2022 J. Phys.: Condens. Matter 34 235101].
 */
class FunctionalESDeltaPlanar : public Functional {
 public:
  /** \brief Standard Constructor */
  FunctionalESDeltaPlanar();
  /** \brief Manual Constructor */
  FunctionalESDeltaPlanar(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor */
  FunctionalESDeltaPlanar(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor */
  ~FunctionalESDeltaPlanar();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal arrays.
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   *  They are always zero due to the required charge neutrality condition.
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
  /** \brief Calculate the energy value of this functional
   *
   *  Calculate the excess free energy of the mean-field electrostatic
   *  Functional. For that, one only multiplies the electrostatic potential with
   *  the total charge density and integrates the result over the entire space.
   *
   *  \return Returns the functional energy value
   */
  virtual double calc_energy();

 private:
  /** \brief System length */
  double length;
  /** \brief Number of grid points */
  size_t grid_count;
  /** \brief Number of grid points in Fourier space (extended system) */
  size_t grid_count_four;
  /** \brief Number of grid points of the extended grid.
   *
   * Since the weighted densities can penetrate the physically relevant system
   * bounaries, an extended grid with more points is established.
   */
  size_t extended_grid_count;
  /** \brief Offset between extended and regular system in terms of grid points
   *
   *  The extended system starts extended_system_offset grid points to the left
   *  of the regular system. To be more precise, the regular system starts at
   *  z=0, while the extended system starts at z=-extended_system_offset*dz.
   */
  size_t extended_system_offset;
  /** \brief Bin width */
  double dz;
  /** \brief Bin width in Fourier space */
  double dkz;
  /** \brief Number of species */
  size_t species_count;
  /** \brief Bjerrum length */
  double bjerrum;
  /** \brief System temperature */
  double temperature;
  /** \brief Dielectric constant */
  double dielectric;
  /** \brief A list with indices of the affected species */
  std::vector<size_t> affected_species;
  /** \brief Hard-sphere diameters of every species */
  std::vector<double> diameters;
  /** \brief Bulk densities of every species */
  std::vector<double> bulk_densities;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Total charge density profile */
  DataFrame<1, double> total_charge_density_profile;
  /** \brief Total charge density profile */
  DataFrame<1, double> total_poisson_rhs;
  /** \brief Total charge density profile */
  DataFrame<1, double> total_potential;
  /** \brief Pointer to density profiles */
  std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Charge density profile of every species */
  std::vector<DataFrame<1, double>> charge_density_profiles;
  /** \brief Right hand side of the Poisson equation for every charge density */
  std::vector<DataFrame<1, double>> poisson_rhs;
  /** \brief Electrostatic potentials of every species */
  std::vector<DataFrame<1, double>> potentials;
  /** \brief Weight functions (delta-functions)
   *
   *  This is a std::vector-matrix, because the weights depend on all
   *  combinations of diameter pairs.
   */
  std::vector<std::vector<DataFrame<1, fftw_complex>>> weights_four;
  /** \brief Weighted (delta-functions) densities in Fourier space
   *
   *  This is a std::vector-matrix, because the weights depend on all
   *  combinations of diameter pairs.
   */
  std::vector<std::vector<DataFrame<1, fftw_complex>>> weighted_densities_four;
  /** \brief Weighted (delta-functions) densities
   *
   *  This is not a std::vector-matrix anymore, because a summation over one
   *  index occured.
   */
  std::vector<DataFrame<1, double>> weighted_densities_real;
  /** \brief Poisson solver for total charge densities
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation.
   */
  PlanarPoissonSolver* total_poisson_solver;
  /** \brief Poisson solver for species potentials
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation.
   */
  PlanarPoissonSolver* poisson_solver;
  /** \brief Functional derivatives */
  double** functional_derivative;
  /** \brief Extract the system Properties required for this functional */
  void extract_system_properties(const Properties& system_properties);
  /** \brief From two of the three electrical properties, the third on can be
   *         calculated.
   */
  void extract_electrical_properties(const Properties& system_properties);
  /** \brief Extract the species specific Properties required for this
   *         functional
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  /** \brief Calculate the properties of the extended system
   *  
   *  The extended system is one hard-sphere diameter larger on both ends
   *  than the regular system. For that the larges hard-sphere diameter is used.
   */
  void calc_extended_system_bounds();
  /** \brief Allocate memory for all DataFrame objects */
  void initialize_all_data_frames();
  /** \brief Initialize the PlanarPoissonSolver object */
  void initialize_poisson_solver();
  /** \brief Initialize the weight functions in Fourier space
   *  
   *  The weight functions are delta-functions with a radius of the average of
   *  two hard-sphere diameters.
   */
  void initialize_weights();
  /** \brief Calculate charge density profiles
   *
   *  calculates the total_charge_density_profile and the charge density
   *  profiles of the species.
   */
  void calc_charge_densities();
  /** \brief Calculate the weighted densities from the charge densities. */
  void calc_weighted_densities();
  /** \brief Calculate the rhs of the Poisson equation for every species */
  void calc_poisson_rhs();
  /** \brief From the charge densities calculate the electrostatic potential */
  void calc_potential();
  /** \brief Integration over charge densities to obtain net charge */
  double calc_net_charge();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_DELTA_PLANAR_HPP_
