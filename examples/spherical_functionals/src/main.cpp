// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
 *  \brief Main file of the example of the spherical functionals.
 *  
 *  This main file contains examples to show how the FunctionalFMTSpherical
 *  class works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/data_field.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/system.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_spherical.hpp"
#include "../../../src/functional_es_mf_spherical.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  // Set the desired system properties
  /* We start by defining the necessary geometric and physical properties.
   * The general properties of a system are put into one Properties conainer.
   * The properties of every species are put in separate Properties container.
   * The containers containing the species properties are encapsulated by an
   * std::vector.
   *
   * We will define three particle species.
   * Two of them will interact via this functional. The other species has no
   * diameter and will be ignored.
   * Which means, that having no diameter is not the same as having a diameter
   * of 0.
   *
   * In general programs should obtain their parameters via the command line
   * or a parameter file (e.g. via the ParameterHandler).
   * However, for this example we hard-coded the numbers to reduce complexity.
   */
// _____________________________________________________________________________
  size_t grid_count = 1001;
  double system_length = 19.821782178217823;  // external pot. between two bins
  double bjerrum_length = 1.;  // not really needed
  // Create objects of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length", system_length);
  system_properties.add_property<double>("bjerrum length", bjerrum_length);
  system_properties.add_property<size_t>("grid count", grid_count);
  // First species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .1);
  properties.add_property<double>("valency", -1.);  // not really needed
  species_properties.push_back(properties);
  properties.clear();
  // Second species
  properties.add_property<double>("bulk density", .2);
  properties.add_property<double>("valency", +1.);  // not really needed
  species_properties.push_back(properties);
  properties.clear();
  // Third species
  properties.add_property<double>("diameter", 1.);
  properties.add_property<double>("bulk density", .3);
  species_properties.push_back(properties);
  properties.clear();
// _____________________________________________________________________________
  /* All the supplied data is now brought together in the System class.
   * The system class creates a density profile (DataField) from the given
   * information and does not allow modifictation of the properties.
   * The System object is then used as input for our functional.
   * In general the functional only copies the parameters it requires to
   * calculate the quantities of interest.
   * It also detects which species interact via this functional and calculates
   * the relevant quantities exclusively for them.
   */
// _____________________________________________________________________________
  System my_system(system_properties, species_properties);
  // Obtain the density profile pointer for modification purposes (e.g. Picard)
  DataField<double>* density_profile = my_system.get_density_profile_pointer();
  // Create FMT Functional object
  // Either specify the species or let the constructor decide which species can
  // interact via this functional
  // FunctionalFMTSpherical my_fmt_functional(&my_system);
  std::vector<size_t> affected_species{0, 2};  // selected species 0 and 2
  FunctionalFMTSpherical my_fmt_functional(&my_system, affected_species);
// _____________________________________________________________________________
  // Picard iteration with the functional derivatives
  /* Now that we processed all the necessary parameters, it is time to calculate
   * a density profile by iteratively calculating the functional derivatives.
   *
   * For that, we first have to specify an external potential, since otherwise
   * our profiles would be completely flat (bulk case).
   *
   * We decided to mimic the shape of a hard spherical capacitor
   * (without any charges). This external potential can also be used to
   * calculate pair correlation functions via the Percus trick, if the inner
   * part and the outer part are far away from each other.
   * From a numerical point of view it is also more stable to limit the system
   * in such a way.
   */
// _____________________________________________________________________________
  DataField<double> fmt_derivatives(species_properties.size(), grid_count);
  double deviation{std::numeric_limits<double>::max()};
  double maximum_deviation{1.0e-6};
  double r{0.};
  double dr{system_length / static_cast<double>(grid_count)};
  double diameter{0.};
  double bulk_density{0.};
  double proposed_density{0.};
  double mixing{1.0e-1};
  size_t step{0};
  std::vector<double> bulk_derivatives(species_properties.size());
  // Set external (hard) potential
  DataField<double> exp_ext_potential(species_properties.size(), grid_count);
  for (size_t i = 0; i != species_properties.size(); ++i) {
    species_properties.at(i).get_property("diameter", &diameter);
    species_properties.at(i).get_property("bulk density", &bulk_density);
    for (size_t j = 0; j != grid_count; ++j) {
      r = dr * static_cast<double>(j + 1);
      if (r < diameter) {
        exp_ext_potential.at(i, j) = 0.;
      } else if (system_length - r < 1.5 * diameter) {
        exp_ext_potential.at(i, j) = 0.;
      } else {
        exp_ext_potential.at(i, j) = 1.;
      }
    }
  }
  // Initial guess for the density profiles
  *(density_profile) = *(density_profile) * exp_ext_potential;
  // Get the bulk derivatives
  my_fmt_functional.calc_bulk_derivative(&bulk_derivatives);
  // Picard iterations
  while (deviation > maximum_deviation) {
    ++step;
    deviation = 0.;
    // Calculate the functional derivative; Suppress warnings about unphysical
    // values
    my_fmt_functional.calc_derivative_no_warnings(&fmt_derivatives);
    for (size_t i = 0; i != species_properties.size(); ++i) {
      species_properties.at(i).get_property("bulk density", &bulk_density);
      for (size_t j = 0; j != grid_count; ++j) {
        // Calculate the right hand side of the update formula of cDFT
        proposed_density = bulk_density * exp_ext_potential.at(i, j) *
            exp(bulk_derivatives.at(i) - fmt_derivatives.at(i, j));
        // Mix the new solution with the old one
        density_profile->at(i, j) = density_profile->at(i, j) * (1. - mixing) +
            mixing * proposed_density;
        // Calculate how much the new density deviates from the old one at most
        if (fabs(density_profile->at(i, j) - proposed_density) > deviation) {
          deviation = fabs(density_profile->at(i, j) - proposed_density);
        }
      }
    }
    // Print the deviation and the number of iteration steps taken
    std::cout << "Picard iteration step " << step;
    std::cout << ". Deviation: " << deviation << std::endl;
  }
  // Write density profile to file
  std::fstream out_stream;
  out_stream.open("spherical_profile.dat", std::ios::out);
  density_profile->set_bin_width(dr);
  density_profile->print(out_stream);
  out_stream.close();
  // Obtain grand potential of the system
  double energy;
  energy = my_fmt_functional.calc_energy();
  std::cout << "Excess free energy of FMT functional: ";
  std::cout << energy << std::endl;
  return 0;
}
