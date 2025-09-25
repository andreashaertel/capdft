// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_FUNCTIONAL_ES_MF_CARTESIAN_HPP_
#define SRC_FUNCTIONAL_ES_MF_CARTESIAN_HPP_
/** \file functional_es_mf_cartesian.hpp
 *  \brief Header file for the FunctionalESMFCartesian class.
 *
 *  The file contains the class declarations of the FunctionalESMFCartesian
 *  class.
 */
// Includes
#include "functional.hpp"  // NOLINT
#include <vector>
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
#include "poisson_solver_cartesian.hpp"  // NOLINT
#include "system.hpp"
/** \brief This class calculates the electrostatic mean field functional in the
 *         in the 3D cartesian geometry.
 *  
 *  This class contains the tools to calculate functional and functional
 *  derivative values of point charges using the mean field approximation in the
 *  3D cartesian geometry. The boundary conditions are defined within the
 *  PoissonSolver which is given as constructor argument.
 */
class FunctionalESMFCartesian : public Functional {
 public:
  /** \brief Standard Constructor
   */
  FunctionalESMFCartesian();
  /** \brief Constructor
   *
   * \param system: defines system and species properties
   * \param poisson_solver: Poisson solver, already initialized with required
   * 		boundary conditions (boundary value zero all over the boundary)
   */
  FunctionalESMFCartesian(
      std::vector<DataFrame<3, double>>* density_profiles,
      System<3>& system, PoissonSolverCartesian* poisson_solver);
  /** \brief Constructor (without System, for planar boundaries)
   *
   * Deprecated version, used in the example
   * ../example/cartesian_functionals_all_tools
   */
  FunctionalESMFCartesian(
      std::vector<DataFrame<3, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties,
      std::vector<size_t> affected_species);
  /** \brief Automated Constructor (without System, for planar boundaries)
   *
   * Deprecated version, used in none of the examples here?
   */
  FunctionalESMFCartesian(
      std::vector<DataFrame<3, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties);

  /** \brief Destructor
   */
  ~FunctionalESMFCartesian();
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
  /** \brief Return the electrostatic potential profile as DataFrame
   */
  void get_potential(System<3>& system, DataFrame<3, double>* data) const;
 private:
  /** \brief Vector that remembers the species, that are affected by this
   *  functional
   */
  std::vector<size_t> affected_species;
  /** \brief Pointer to density profiles
   */
  const std::vector<DataFrame<3, double>>* density_profiles_pointer;
  /** \brief System length (radius of the sperical geometry)
   */
  std::vector<double> lengths;
  /** \brief Number of grid points
   *
   * voxel_count is the total ammount of voxels in the system, while
   * grid_counts holds the number of grid points in every dimension.
   */
  size_t voxel_count;
  std::vector<size_t> grid_counts;
  /** \brief Periodic boundary conditions
   *
   *  Contains tree booleans, that are "true" if the corresponding dimension
   *  has periodic boundary conditions and false otherwise.
   *  Note, that at least one has to be non-periodic, since this required by
   *  the Poisson solver.
   */
  std::vector<bool> periodic_boundaries;
  /** \brief Bin sizes
   */
  std::vector<double> bin_sizes;
  /** \brief Bjerrum length
   */
  double bjerrum;
  /** \brief System temperature
   */
  double temperature;
  /** \brief Dielectric constant
   */
  double dielectric;
  /** \brief Valency of every species
   *
   *  To be more general, valencies must be given as "double".
   */
  std::vector<double> valencies;
  /** \brief Number of species
   */
  size_t species_count;
  /** \brief Total charge density profile
   */
  std::vector<double> charge_density_profile;
  /** \brief Right hand side of the Poisson equation
   */
  std::vector<double> poisson_rhs;
  /** \brief Electrostatic potential (numerical solution of Poisson equation)
   */
  std::vector<double> potential;
  /** \brief Poisson solver
   *
   *  This object contains the matrix representation of the numerical Poisson
   *  equation in the 3D cartesian geometry. This can be either the
   *  CartesianPoissonSolver (planar boundaries) or the
   *  CartesianPoissonSolverAny (structured boundaries).
   */
  PoissonSolverCartesian* poisson_solver;
  /** \brief Extract the system Properties required for this functional
   *
   * (first version is deprecated, needed for the deprecated constructors)
   */
  void extract_system_properties(const Properties& system_properties);
  void extract_system_properties(System<3>& system);
  /** \brief From two of the three electrical properties, the third on can be
   *         calculated.
   */
  void extract_electrical_properties(const Properties& system_properties);
  /** \brief Extract the species specific Properties required for this
   *         functional
   *
   * (first version is deprecated, needed for the deprecated constructors)
   */
  void extract_species_properties(
      const std::vector<Properties>& species_properties);
  void extract_species_properties(System<3>& system);
  /** \brief Allocate memory for all DataFrame objects */
  void initialize_all_data_frames();
  /** \brief Initialize the PlanarPoissonSolver object if no other PoissonSolver
   * is given in the constructor
   *
   * (deprecated, needed for the deprecated constructors)
   */
  void initialize_poisson_solver();
  /** \brief Calculate the net charge density profile and the rhs of the Poisson
   *         equation.
   */
  void calc_charge_densities();
  /** \brief From the charge densities calculate the electrostatic potential */
  void calc_potential();
};
#endif  // SRC_FUNCTIONAL_ES_MF_CARTESIAN_HPP_
