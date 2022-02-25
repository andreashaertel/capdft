// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_SYSTEM_HPP_
#define SRC_SYSTEM_HPP_
/** \file src/system.hpp
 *  \brief This file contains the declarations of the System class
 *
 */
#include <vector>
#include "data_field.hpp"
#include "data_frame.hpp"
#include "properties.hpp"
/** \brief Container class that contains all Properties and DataFields
 *
 */
class System {
 public:
  /** \brief Constructor
   *
   */
  //System();
  template <typename AnyDataFrame>
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
  void bulk(); //VORSCHLAG: Aendern in: set_bulk_density()
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
  void set_fugacities(std::vector<double>* fugacities);
  /** \brief Obtain pointer of the density_profile
   *
   */
  DataField<double>* get_density_profile_pointer();

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
  DataField<double>* density_profile;
  std::vector<DataFrame> density_profiles;

 protected:
};
#endif  // SRC_SYSTEM_HPP_
