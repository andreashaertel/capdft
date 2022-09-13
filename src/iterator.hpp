// SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_ITERATOR_HPP_
#define SRC_ITERATOR_HPP_
/** \file src/iterator.hpp
 *  \brief This file contains the declarations of the Iterator class
 *
 *  The class holds the system and all functionals defined on the system. 
 *  It provides methods to iterate the system's density profiles and to perform
 *  calculations in the framework of DFT. 
 *
 */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <vector>
#include <map>
#include "convergence_criterion.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "functional.hpp" // NOLINT
#include "properties.hpp"  // NOLINT
#include "system.hpp"  // NOLINT
/** \brief Class that holds the system and all functionals and provides DFT 
 *         methods
 *
 *  The class is a container for the physical system which have to be provided 
 *  during declaration. The system particulalry holds all properties and stores
 *  the density profiles. 
 *  The class further provides methods for calculations in the framework of DFT.
 *
 */
class Iterator {
 public:
  /** \brief Constructor
   * 
   *  \param density_profiles The density profiles of different species given as
   *         DataFrame
   *  \param external_potentials The external potentials for each species
   *
   *  Get the density profiles pointer and the external potentials.
   *  The system data is provided to the functionals and is not needed here.
   *
   */
  // TODO(Andreas): add a functional_ideal object
  Iterator(
      std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<DataFrame<1, double>>& exp_external_potentials,
      const std::vector<Properties>& species_properties);
  /** \brief Destructor
   *
   */
  ~Iterator();
  /** \brief Add an excess functional
   *
   *  This function adds the specified functional to a list (std::vector).
   *  This works for any capdft functional, because they are derived from the
   *  abstract Functional base class.
   *
   *  For example, the function can be used as follows: 
   *  Iterator my_iterator(my_densities, my_external_potentials);
   *  FunctionalFMTSpherical my_functional(...);
   *  my_iterator.add_excess_functional(my_functional);
   *
   */
  void add_excess_functional(Functional* functional);
  /** \brief Remove a functional according to its index. 
   * 
   *  Removes the specified functional from the DFT framework. 
   *  If the index specifies no functional, nothing will happen. 
   * 
   *  \param index Index of the functional that shall be removed.
   *
   *  \return Indicates whether a functional has been removed or not. 
   *
   */
  void remove_excess_functional(size_t index);
  /** \brief Remove all functionals
   *
   */
  void clear_functionals();
  /** \brief Iterates the system densities according to the Picard iteration
   *         method.
   *
   *  \param Mixing factor of the Picard iteration scheme.
   */
  void run_picard(double mixing);
  /** \brief Iterates the system densities according to the Andersen mixing
   *         method.
             https://doi.org/10.1016/j.fluid.2017.03.023
   *
   *  \param Mixing factor of the Picard iteration scheme.
   *  \param memory is the number of density profiles from the past iterations
   *         the scheme will remember.
   */
  void run_anderson(double mixing, size_t memory);
  /** \brief Calculate the excess free energy.
   *
   *  The grand canonical energy of the system is calculated for its current
   *  state without further iteration. 
   *
   *  \return The grand canonical energy in kT per system volume.
   *
   */
  double calculate_excess_free_energy();
  /** \brief Calculate the grand canonical energy of the system. 
   *
   *  The grand canonical energy of the system is calculated for its current
   *  state without further iteration. 
   *
   *  \return The grand canonical energy in kT per system volume.
   *
   */
  double calculate_gc_energy();
  /** \brief Add a convergence criterion
   *
   *  The standard criterion (ConvergenceCriterionMaxDev) is added when a new
   *  Iterator is created.
   *
   */
  template <typename T>
  void add_convergence_criterion(double threshold) {
    ConvergenceCriterion* criterion = 
        new T(*density_profiles, proposed_densities, threshold);
    convergence_criteria.push_back(criterion);
  }
  template <typename U>
  void add_convergence_criterion(int threshold_int) {
    ConvergenceCriterion* criterion =
        new U(*density_profiles, proposed_densities, threshold_int);
    convergence_criteria.push_back(criterion);
  }
  /** \brief Clear convergence criteria
   *
   *  Clears the convergence criteria. This is especially important when
   *  overwriting the standard criteria. If run_*() is executed without a
   *  ConvergenceCriterion, it will run forever.
   *
   */
  void clear_convergence_criteria();

 private:
  /** \brief Density profiles
   *
   */
  std::vector<DataFrame<1, double>>* density_profiles;
  /** \brief Proposed density profiles by Picard iteration
   *
   */
  std::vector<DataFrame<1, double>> proposed_densities;
  /** \brief External potentials
   *
   */
  const std::vector<DataFrame<1, double>>* exp_external_potentials;
  /** \brief Species properties pointer
   *
   *  Species properties are required for the "bulk density" property.
   *
   */
  const std::vector<Properties>* species_properties;
  /** \brief Container holding all functionals that have benn added
   *
   */
  std::vector<Functional*> excess_functionals;
  /** \brief Container holding all convergence criteria to be checked
   *
   */
  std::vector<ConvergenceCriterion*> convergence_criteria;
  /** \brief Functional derivatives
   *
   */
  std::vector<std::vector<DataFrame<1, double>>> functional_derivatives;
  /** \brief Bulk functional derivatives
   *
   */
  std::vector<std::vector<double>> bulk_derivatives;
  /** \brief Check all convergence criteria and return true if one is fulfilled
   *
   */
  bool check_convergence_criteria();
  /** \brief Calculates all possible scalar product combinations of a profile
   *         history
   *
   *  \param profile_history is a vector of density profiles from which all
   *         scalar product combinations are calculated
   *
   *  \return A matrix that contains all scalar products.
   */
  std::vector<std::vector<double>> calc_scalar_products(
      std::vector<std::vector<DataFrame<1, double>>>& profile_history);
  /** \brief Calculates the linear combination of vectors with scalar_products
   *         that minimizes the 2-norm under the condition that the sum of all
   *         coefficients equals 1.
   *
   *  \param scalar_products is a matrix of all scalar products of the vectors
   *         which are used in the linear combination.
   *
   *  \return Linear combination coefficients
   */
  std::vector<double> shortest_linear_combination(
      std::vector<std::vector<double>>& scalar_products);

 protected:
};
/** \brief Function that defines the root search problem that arises from the
 *         Lagrange multiplier formalism in a way that GSL can process it.
 */
int anderson_f(const gsl_vector* x, void* params, gsl_vector* f);
#endif  // SRC_ITERATOR_HPP_
