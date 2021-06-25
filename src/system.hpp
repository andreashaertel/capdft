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
#include "properties.hpp"
/** \brief Container class that contains all Properties and DataFields
 *
 */
class System {
 public:
  /** \brief Constructor
   *
   */
  System();
  System(
      Properties system_properties,
      std::vector<Properties> species_properties,
      DataField density_profile);
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

 private:
  /** \brief Supplied system properties
   *
   */
  Properties system_properties;
  /** \brief Supplied species properties
   *
   */
  std::vector<Properties> species_properties;
  /** \brief Data field that contains the density profile
   *
   */
  DataField density_profile;

 protected:
};
#endif  // SRC_SYSTEM_HPP_
