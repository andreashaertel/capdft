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
  /** \brief Standard Constructor
   */
  FunctionalESDeltaPlanar();
  /** \brief Manual Constructor
   */
  FunctionalESDeltaPlanar(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor
   */
  FunctionalESDeltaPlanar(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);
  /** \brief Destructor
   */
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
  /** \brief System length
   */
  double length;
  /** \brief Number of grid points
   */
  size_t grid_count;
  /** \brief Bin width
   */
  double dz;
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
  /** \brief Hard-sphere diameters of every species
   */
  std::vector<double> diameters;
  /** \brief Bulk densities of every species
   */
  std::vector<double> bulk_densities;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Pointer to density profiles
   */
  std::vector<DataFrame<1, double>>* density_profiles_pointer;
  /** \brief Charge density profile of every species
   */
  std::vector<DataFrame<1, double>> charge_density_profiles;
  /** \brief Right hand side of the Poisson equation for every charge density
   */
  std::vector<DataFrame<1, double>> poisson_rhs;
  /** \brief Electrostatic potentials of every species
   */
  std::vector<DataFrame<1, double>> potentials;
  /** \brief Weighted (delta-functions) densities
   */
  std::vector<DataFrame<1, double>> weighted_densities;
  /** \brief Poisson solver
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
  /** \brief Allocate memory for all DataFrame objects */
  void initialize_all_data_frames();
  /** \brief Initialize the PlanarPoissonSolver object */
  void initialize_poisson_solver();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   */
  void calc_charge_densities();
  /** \brief From the charge densities calculate the electrostatic potential */
  void calc_potential();
  /** \brief Integration over charge densities to obtain net charge */
  double calc_net_charge();

 protected:
};
#endif  // SRC_FUNCTIONAL_ES_DELTA_PLANAR_HPP_
