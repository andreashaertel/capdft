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
#include "data_field.hpp"  // NOLINT
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
  // TODO(Andreas): input as reference(&), output as pointers(*) (from Moritz)
  explicit Dft(System* system);
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
  /** \brief Sets the fugacity of each species in the system according to the 
   *         specified bulk densities. 
   *
   *  The fugacities for all species are determined from the Euler equation of 
   *  DFT, which in three dimensions read 
   *  \f[\rho_\nu(\vec{r}) = \frac{e^{\beta\mu_\nu}}{\Lambda_\nu^3}
        \exp(c^{(1)}_\nu(\vec{r})-\beta v_\nu(\vec{r})) . \f] 
   *  The fugacities \f$z_\nu\f$ are defined such that the Euler equation 
   *  follows as
   *  \f[\rho_\nu(\vec{r}) = z_\nu 
        \exp(c^{(1)}_\nu(\vec{r})-\beta v_\nu(\vec{r})) . \f]
   *  Thus, in bulk (vanishing external potentials) we have 
   *  \f[z_\nu=\rho_\nu/\exp(c^{(1)}_\nu) .\f]
   *  Here, the bulk densities and functional derivatives (via \f$c_\nu^{(1)}\f$
   *  enter. 
   *
   *  Then the species properties are updated accordingly with the new 
   *  fugacities. This method should be called after adding all excess 
   *  functionals to the Dft class. 
   *
   *  \param bulk_densities A vector that defines the bulk densities of all 
   *         species. 
   *
   */
  void set_fugacities_from_bulk_densities(std::vector<double>* bulk_densities);
  /** \brief Sets the fugacity of each species in the system according to the
   *         bulk densities specified in the species properties. 
   *
   *  The fugacities for all species are determined from the Euler equation of 
   *  DFT, which in three dimensions read 
   *  \f[\rho_\nu(\vec{r}) = \frac{e^{\beta\mu_\nu}}{\Lambda_\nu^3}
        \exp(c^{(1)}_\nu(\vec{r})-\beta v_\nu(\vec{r})) . \f] 
   *  The fugacities \f$z_\nu\f$ are defined such that the Euler equation 
   *  follows as
   *  \f[\rho_\nu(\vec{r}) = z_\nu 
        \exp(c^{(1)}_\nu(\vec{r})-\beta v_\nu(\vec{r})) . \f]
   *  Thus, in bulk (vanishing external potentials) we have 
   *  \f[z_\nu=\rho_\nu/\exp(c^{(1)}_\nu) .\f]
   *  Here, the bulk densities and functional derivatives (via \f$c_\nu^{(1)}\f$
   *  enter. The bulk densities are red for each species from the property 
   *  bulk_density. 
   *
   *  Then the species properties are updated accordingly with the new
   *  fugacities. This method should be called after adding all excess 
   *  functionals to the Dft class. 
   *
   */
  void set_fugacities_from_bulk_densities();
  /** \brief Iterates the system densities according to the set iteration 
   *         method. 
   *
   *  Performs one iteration step according to the set iteration method. 
   *  The iteration method also defines the deviation measure between the 
   *  old density profiles and the newly calculated ones. <br>
   *  The fugacity of each species (property fugacity) is used during the 
   *  iteration according to the Euler equation. 
   *
   *  TODO: Picard mixing parameter, Check for numerical correct solution, and
   *  more must still be incorporated. 
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
  std::map<size_t, Functional*> functional;
  /** \brief Functional index of the last add functional, initial value is 0. */
  size_t functional_index;

 protected:
};
#endif  // SRC_DFT_HPP_
