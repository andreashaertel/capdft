// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_ITERATOR_HPP_
#define SRC_ITERATOR_HPP_
/** \file src/dft.hpp
 *  \brief This file contains the declarations of the Iterator class
 *
 *  The class holds the system and all functionals defined on the system. 
 *  It provides methods to iterate the system's density profiles and to perform
 *  calculations in the framework of DFT. 
 *
 */
#include <vector>
#include <map>
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
  void add_excess_functional(Functional& functional);
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
  /** \brief Iterates the system densities according to the set iteration 
   *         method. 
   *
   *  \return Returns the deviation between the new and old density profiles. 
   *
   */
  double run();
  /** \brief Calculate the grand canonical energy of the system. 
   *
   *  The grand canonical energy of the system is calculated for its current
   *  state without further iteration. 
   *
   *  \return The grand canonical energy in kT per system volume.
   *
   */
  double calculate_gc_energy();

 private:
  /** \brief Density profiles
   *
   */
  std::vector<DataFrame<1, double>>* density_profiles;
  /** \brief External potentials
   *
   */
  const std::vector<DataFrame<1, double>>* exp_external_potentials;  // TODO(Moritz): maybe not pointer? // NOLINT
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
  /** \brief Functional derivatives
   *
   */
  std::vector<std::vector<DataFrame<1, double>>> functional_derivatives;
  /** \brief Bulk functional derivatives
   *
   */
  std::vector<std::vector<double>> bulk_derivatives;

 protected:
};
#endif  // SRC_ITERATOR_HPP_
