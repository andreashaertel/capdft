// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_SYSTEM_HPP_
#define SRC_SYSTEM_HPP_
/** \file src/system.hpp
 *  \brief This file contains the declarations of the System class
 *
 */
#include <vector>
#include "data_field.hpp"  // NOLINT
#include "df_spherical.hpp"  // NOLINT
#include "data_frame.hpp"  // NOLINT
#include "properties.hpp"  // NOLINT
/** \brief Template Container class that contains all system Properties and
 *         density profiles stored in data fields (DF*)
 *
 */
template <typename T>  // Template for different data frames DF*
class System {
 public:
  /** \brief Constructor
   *
   */
  System();
  System(
      Properties system_properties,
      std::vector<Properties> species_properties);
  /** \brief Destructor
   *
   */
  ~System();
  /** \brief Access system properties (no modification)
   *
   */
  const Properties& get_system_properties() const;
  /** \brief Access species properties (no modification)
   *
   */
  const std::vector<Properties>& get_species_properties() const;
  /** \brief Write bulk values to the density profile
   *
   */
  void set_bulk_densities();
  /** \brief Set/Update fugacities in all species properties. 
   *
   *  Updates the property fugacity in each species_properties 
   *  element in System. The fugacity \f$z_\nu\f$ of a species \f$\nu\f$ is 
   *  defined via the Euler equation 
   *  \f[\rho_\nu(\vec{r}) = z_\nu 
        \exp(c^{(1)}_\nu(\vec{r})-\beta v_\nu(\vec{r})) . \f]
   *  Accordingly, in three dimensions the fugacity of a species is defined by 
   *  \f[z_nu=e^{\beta\mu_\nu}/\Lambda_\nu^3. <br>
   *
   *  \param fugacities Vector of fugacities for all species. 
   *
   */
  void set_fugacities(std::vector<double>& fugacities);
  /** \brief Obtain pointer of the density_profiles
   *
   */
  std::vector<T>* get_density_profiles_pointer();
  /** \brief Get current density profiles by reference. 
   *
   *  The density profiles are constant and cannot be changed. In order to 
   *  update them, use update_density_profiles(). 
   *
   *  \return Constant vector of DF* by reference. The vector is over 
   *          species in the system. 
   *
   */
  std::vector<T>& get_density_profiles();
  /** \brief Update density profiles to new values. 
   * 
   *  Updates the density profiles to new values. 
   *
   *  \param density_profiles The new density profiles as DF*s in a 
   *         vector over the species of the system. The data is not changed,
   *         values are just copied. 
   *
   */
  void update_density_profiles(
      const std::vector<T>& other_density_profiles);

 private:
  /** \brief Supplied system properties
   *
   */
  Properties system_properties;
  /** \brief Supplied species properties
   *
   */
  std::vector<Properties> species_properties;
  /** \brief Data field pointer that contains the density profile if initialized
   *
   */
  std::vector<T> density_profiles;

 protected:
};
#endif  // SRC_SYSTEM_HPP_
