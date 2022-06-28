// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
#define SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
/** \file functional_es_delta_spherical.hpp
 *  \brief Header file for the FunctionalESDeltaSpherical class.
 *
 *  The file contains the declarations of the FunctionalESDeltaSpherical class.
 */
// Includes
#include "functional.hpp"  // NOLINT
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
#include "radial_poisson_solver.hpp"  // NOLINT
/** \brief This class calculates the electrostatic interactions via the delta
 *         functional.
 *  
 *  The theory for this functional was published in
 *  [M. Bültmann and A. Härtel 2022 J. Phys.: Condens. Matter 34 235101].
 */
class FunctionalESDeltaSpherical : public Functional {
 public:
  /** \brief Standard Constructor
   */
  FunctionalESDeltaSpherical();
  /** \brief Manual Constructor
   */
  FunctionalESDeltaSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor
   */
  FunctionalESDeltaSpherical(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   */
  ~FunctionalESDeltaSpherical();
  /** \brief Calculate the functional derivatives
   *
   *  This function calculates the functional derivatives and updates the
   *  corresponding internal arrays.
   *
   *  It uses the functions:
   *  calc_weighted_densities(),
   */
  virtual void calc_derivative(
      std::vector<DataFrame<1, double>>* functional_derivative);
  /** \brief Calculate bulk derivatives
   *
   */
  virtual void calc_bulk_derivative(std::vector<double>* bulk_derivative);
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
   */
  double length;
  /** \brief Number of grid points
   */
  size_t grid_count;
  /** \brief Bin width
   */
  double dr;
  /** \brief Bin width in Fourier space
   */
  double dkr;
  /** \brief Number of species
   */
  size_t species_count;
  /** \brief Bjerrum length
   */
  double bjerrum;
  /** \brief System temperature
   */
  double temperature;
  /** \brief Dielectric constant
   */
  double dielectric;
  /** \brief A list with indices of the affected species
   */
  std::vector<size_t> affected_species;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Hard-sphere diameters
   */
  std::vector<double> diameters;
  /** \brief Bulk densities
   */
  std::vector<double> bulk_densities;
  /** \brief Pointer to density profiles
   */
  std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Charge density profile for every species
   */
  std::vector<DataFrame<1, double>> charge_density_profiles;
  /** \brief Right hand side of the Poisson equation for each species
   */
  std::vector<DataFrame<1, double>> poisson_rhs;
  /** \brief Electrostatic potentials (numerical solution of Poisson equation)
   *         for each species.
   */
  std::vector<DataFrame<1, double>> potentials;
  /** \brief Fourier transformed delta weight functions
   */
  std::vector<std::vector<DataFrame<1, double>>> weights_delta;
  /** \brief Weighted densities (delta-functions)
   */
  std::vector<std::vector<DataFrame<1, double>>> weighted_densities;
  /** \brief Poisson solver
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation.
   */
  RadialPoissonSolver* poisson_solver;
  /** \brief Functional derivatives
   */
  double** functional_derivative;
  /** \brief Extract the system Properties required for this functional
   */
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
  /** \brief Allocate memory for all DataFrame objects
   */
  void initialize_all_data_frames();
  /** \brief Initialize the RadialPoissonSolver object
   */
  void initialize_poisson_solver();
  /** \brief Initialize Fourier transformed weight functions
   */
  void initialize_weights();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   */
  void calc_charge_densities();
  /** \brief Calculate the right-hand side of the Poisson equation for every
   *         species.
   */
  void calc_poisson_rhs();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   *
   *  See equation (32a).
   */
  void calc_weighted_densities();
  /** \brief From the charge densities calculate the electrostatic potential
   *
   *  Unlike the mean-field electrostatic functional, the potential of every
   *  species is calculated separately. Moreover, this functional uses weighted
   *  densities instead of the regular charge densities.
   *  The boundary conditions for the solution of the Poisson equations are:
   *  - inner one equals 0 (Neumann) due to the radial symmetry
   *  - outer one equals net charge divided by radial position (Dirichlet)
   *    due to Gauss' theorem.
   *  Note, that there is no external charge at the center.
   *
   */
  void calc_potentials();
  /** \brief Integration over weighted charge densities to obtain net charge
   */
  std::vector<double> calc_charge_weight_dens();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_DELTA_SPHERICAL_HPP_
