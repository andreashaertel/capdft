// SPDX-FileCopyrightText: 2022 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_DFT_HPP_
#define SRC_DFT_HPP_
/** \file src/dft.hpp
 *  \brief This file contains the declarations of the Dft class
 *
 *  The class holds the system and all functionals defined on the system. 
 *  It provides methods to iterate the system's density profiles and to perform
 *  calculations in the framework of DFT. 
 *
 */
#include <vector>
#include <map>
#include "data_field.hpp"
#include "properties.hpp"
#include "system.hpp"
/** \brief Class that holds the system and all functionals and provides DFT 
 *         methods
 *
 *  The class is a container for the physical system which have to be provided 
 *  during declaration. The system particulalry holds all properties and stores
 *  the density profiles. 
 *  The class further provides methods for calculations in the framework of DFT.
 *
 */
class Dft {
 public:
  /** \brief Constructor
   * 
   *  \param system The system on which the Dft class shall work. 
   *
   *  The physical system must be provided during declaration and cannot be 
   *  changed at any later time. 
   *
   */
  Dft(System* system);
  /** \brief Destructor
   *
   */
  ~Dft();
  /** \brief Add a functional to work on the system
   *
   *  This template function creates a new functional of the specified type 
   *  AnyFunctional and stores it to be applied to the system. 
   *
   *  For example, the function can be used as follows: 
   *  Dft* my_dft = new Dft(system);
   *  my_dft->add_excess_functional<FunctionalFMTSpherical>();
   *
   *  Note that each functional that shall be available must be explicitly 
   *  defined in the Dft.cpp code file. 
   *
   *  \return Index of the functional that has been add. 
   *
   */
  template <typename AnyFunctional>
  size_t add_excess_functional();
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
  bool remove_excess_functional(size_t index);
  /** \brief Sets the chemical potentials of the system according to the 
   *         specified bulk densities. 
   *
   *  The chemical potentials for all species are determined such that the 
   *  given bulk densities minimize the functional. 
   *  Then the species properties are updated accordingly with the new chemical
   *  potentials. This method should be called after adding all excess 
   *  functionals to the Dft class. 
   *
   *  \param bulk_densities A vector that defines the bulk densities of all 
   *         species. 
   *
   */
  void set_chempots_from_bulk_densities(std::vector<double>* bulk_densities);
  /** \brief Sets the chemical potentials of the system according to the bulk
   *         densities specified in the species properties. 
   *
   *  The chemical potentials for all species are determined such that the bulk
   *  densities minimize the functional which are given in the species 
   *  properties. The field bulk_density is expected in the properties vector
   *  of each species of the system. 
   *
   *  Then the species properties are updated accordingly with the new chemical
   *  potentials. This method should be called after adding all excess 
   *  functionals to the Dft class. 
   *
   */
  void set_chempots_from_bulk_densities();
  /** \brief Iterates the system densities according to the set iteration 
   *         method. 
   *
   *  Performs one iteration step according to the set iteration method. 
   *  The iteration method also defines the deviation measure between the 
   *  old density profiles and the newly calculated ones. 
   * 
   *  \return Returns the deviation between the new and old density profiles. 
   *
   */
  double iterate_densities();
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
  /** \brief Physical system holding properties and density DataFields. */
  System* system;
  /** \brief Functionals that have been add. */
  std::map<size_t,Functional*> functional;
  /** \brief Functional index of the last add functional, initial value is 0. */
  size_t functional_index;

 protected:
};
#endif  // SRC_DFT_HPP_
