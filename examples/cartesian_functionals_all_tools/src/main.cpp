// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/cartesian_functionals_all_tools/src/main.cpp
 *  \brief Main file of the example of the 3D cartesian functionals.
 *  
 *  This main file contains examples to show how Functional classes in the
 *  3D cartesian geometry work.
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../../src/convergence_criterion.hpp"
#include "../../../src/convergence_criterion_max_dev.hpp"
#include "../../../src/convergence_criterion_steps.hpp"
#include "../../../src/convergence_criterion_nan.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/functional.hpp"
#include "../../../src/functional_fmt_cartesian.hpp"
#include "../../../src/iterator.hpp"
#include "../../../src/properties.hpp"
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
   *
   * In general programs should obtain their parameters via the command line
   * or a parameter file (e.g. via our ParameterHandler).
   * However, for this example we hard-coded the numbers to reduce complexity.
   *
   * Note, that the system_length is chosen such that the external potential
   * step lies exactly between two bins.
   */
// _____________________________________________________________________________
<<<<<<< HEAD
  std::vector<size_t> grid_counts{128, 128, 128};  // (x, y, z)
=======
  std::vector<size_t> grid_counts{63, 63, 63};  // (x, y, z)
  //std::vector<size_t> grid_counts{1000, 1, 1};  // (x, y, z)
  //std::vector<size_t> grid_counts{1, 1000, 1};  // (x, y, z)
  //std::vector<size_t> grid_counts{1, 1, 1000};  // (x, y, z)
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
  std::vector<double> system_lengths{2., 2., 2.};  // in nm (x, y, z)
  // Create objects of Properties class
  Properties properties;
  Properties system_properties;
  std::vector<Properties> species_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length x", system_lengths.at(0));
  system_properties.add_property<double>("length y", system_lengths.at(1));
  system_properties.add_property<double>("length z", system_lengths.at(2));
  system_properties.add_property<size_t>("grid count x", grid_counts.at(0));
  system_properties.add_property<size_t>("grid count y", grid_counts.at(1));
  system_properties.add_property<size_t>("grid count z", grid_counts.at(2));
  // One hard-sphere species
  properties.add_property<double>("diameter", .3);
  properties.add_property<double>("bulk density", 7.);  // 0.1 volume fraction
  species_properties.push_back(properties);
  properties.clear();
// _____________________________________________________________________________
  /* The functional objects are created. As input they require the Properties of
   * the system and the std::vector<Properties> from the species. Moreover
   * the pointer of the three-dimensional density profiles needs to be passed
   * on, as well as a std::vector<size_t> with the species numbers, that are
   * affected by the functional.
   * 
   * Remark: Even though the DataFrames (e.g. density_profiles) are
   * three-dimensional, the underlying data arrays are still one-dimensional.
   * The only difference are the 3-dim. helper functions like
   * at(size_t,size_t,size_t).
   */
// _____________________________________________________________________________
  // Create density profiles
  std::vector<DataFrame<3, double>> density_profiles(0);
  for (size_t i = 0; i < species_properties.size(); ++i) {
    density_profiles.push_back(DataFrame<3, double>(grid_counts));
  }
  // Create an FMT Functional object. For this we specify the
  // species which are interacting via this functional in affected_species.
  std::vector<size_t> affected_species_fmt{0};  // selected species 0
  FunctionalFMTCartesian my_fmt_functional(&density_profiles,
      species_properties, system_properties, affected_species_fmt);
// _____________________________________________________________________________
  // Picard iterations
  /* For the Picard iterations the Iterator class is used. For that we define
   * the external potential, add our functionals "my_fmt_functional" and
   * "my_es_functional" to the functional list and use the run_picard() function
   * to carry out the Picard iterations. In this case we use two planar hard
   * walls at distance system_length as external potential.
   *
   * If you want to save loads of time use the Anderson mixing algorithm
   * run_andersen() instead of the Picard iterations. They usually are faster
   * by a factor of 20.
   *
   * The convergence criteria are also added to the Iterator. They are added via
   * a template argument that specifies the type of criterion and a function
   * argument that specifies the threshold of the specified criterion.
   *
   * For example: 
   * my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-4);
   * This will terminate the iterations after the largest difference
   * (ConvergenceCriterionMaxDev) of the old density profile and the new one is
   * smaller than 1.0e-4.
   */
// _____________________________________________________________________________
  // Create external potential DataFrames
  double diameter{0.}, bulk_density{0.};
  double x{0.}, y{0.}, z{0.};
<<<<<<< HEAD
  double frequency{4. * M_PI / system_lengths.at(0)};
=======
  double frequency{2. * M_PI / system_lengths.at(0)};
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
  double wave_left{0.}, wave_right{0.};;
  double dx{system_lengths.at(0) / static_cast<double>(grid_counts.at(0))};
  double dy{system_lengths.at(1) / static_cast<double>(grid_counts.at(1))};
  double dz{system_lengths.at(2) / static_cast<double>(grid_counts.at(2))};
  std::vector<DataFrame<3, double>> exp_ext_potential(
      species_properties.size(), DataFrame<3, double>(grid_counts));
  // Set external hard potential: two plates with sine wave shape
<<<<<<< HEAD
  // The left (z=0) plate is wave shaped in the x direction, while the right
  // (z=L) plate is wave shaped in the y direction.
=======
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
  for (auto& species : affected_species_fmt) {
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t i = 0; i != grid_counts.at(0); ++i) {
      x = dx * static_cast<double>(i);
      for (size_t j = 0; j != grid_counts.at(1); ++j) {
        y = dy * static_cast<double>(j);
        for (size_t k = 0; k != grid_counts.at(2); ++k) {
          z = dz * static_cast<double>(k);
<<<<<<< HEAD
          wave_left = .5 + 2. * pow(sin(frequency * x), 2);
          wave_right =.5 + 2. * pow(sin(frequency * y), 2);
=======
          wave_left = 1 + .5 * sin(frequency * x);
          wave_right = 1 + .5 * sin(frequency * y);
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
          if (z < diameter * wave_left) {
            exp_ext_potential.at(species).at(i, j, k) = 0.;
          } else if ((system_lengths.at(2) - z) < diameter * wave_right) {
            exp_ext_potential.at(species).at(i, j, k) = 0.;
          } else {
            exp_ext_potential.at(species).at(i, j, k) = 1.;
          }
        }
      }
    }
  }
  // Initial guess for the density profiles
  for (size_t i = 0; i < species_properties.size(); ++i) {
    species_properties.at(i).get_property("bulk density", &bulk_density);
    density_profiles.at(i).set_all_elements_to(bulk_density);
    density_profiles.at(i) *= exp_ext_potential.at(i);
  }
  // Create iterator and run iterations
  Iterator my_iterator(&density_profiles, exp_ext_potential,
      species_properties);
  my_iterator.add_excess_functional(&my_fmt_functional);
  my_iterator.clear_convergence_criteria();
  my_iterator.add_convergence_criterion<ConvergenceCriterionSteps>(2e3);
<<<<<<< HEAD
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(1.0e-6);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  //my_iterator.run_picard(1e-1);
  my_iterator.run_anderson(1e-1, 20);
=======
  my_iterator.add_convergence_criterion<ConvergenceCriterionMaxDev>(-1.0e-4);
  my_iterator.add_convergence_criterion<ConvergenceCriterionNan>(0);
  my_iterator.run_picard(1e-1);
  //my_iterator.run_anderson(1e-1, 10);
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
// _____________________________________________________________________________
  /* All done!
   * Now we produce some output and view it in gnuplot.
   * We also supplied this example with a pdf that shows the plot in case you do
   * not use gnuplot.
   *
   * You can also try to reproduce the density profiles of
   * [M. Bültmann and A. Härtel 2022 J. Phys.: Condens. Matter 34 235101].
   */
// _____________________________________________________________________________
<<<<<<< HEAD
  // Write density profile to file
  std::fstream out_stream;
  out_stream.open("3d_profile.dat", std::ios::out);
  out_stream << "# [x] [y] [z] [density profile]" << std::endl;
  for (size_t i = 0; i < grid_counts.at(0); ++i) {
    x = dx * static_cast<double>(i);
    for (size_t j = 0; j < grid_counts.at(1); ++j) {
      y = dy * static_cast<double>(j);
      for (size_t k = 0; k < grid_counts.at(2); ++k) {
        z = dz * static_cast<double>(k);
        out_stream << x << " " << y << " " << z << " ";
        out_stream << density_profiles.at(0).at(i, j, k);
        out_stream << std::endl;
      }
    }
  }
  out_stream.close();
  // Obtain grand potential of the system
  double energy;
  energy = my_fmt_functional.calc_energy();
  std::cout << "Excess free energy of the FMT functional: ";
  std::cout << energy << std::endl;
=======
  //// Write density profile to file
  //std::fstream out_stream;
  //out_stream.open("planar_profile.dat", std::ios::out);
  //for (size_t i = 0; i < grid_count; ++i) {
  //  z = dz * static_cast<double>(i);
  //  out_stream << z << " ";
  //  for (size_t j = 0; j < species_properties.size(); ++j) {
  //    out_stream << density_profiles.at(j).at(i) << " ";
  //  }
  //  out_stream << std::endl;
  //}
  //out_stream.close();
  //// Obtain grand potential of the system
  //double energy;
  //energy = my_fmt_functional.calc_energy();
  //std::cout << "Excess free energy per square nanometer of FMT functional: ";
  //std::cout << energy << std::endl;
>>>>>>> 7277f8181f7df11e7866e0736c890ae3e0b65b8a
  //energy = my_es_functional.calc_energy();
  //std::cout << "Excess free energy of mean-field electrostatic functional: ";
  //std::cout << energy << std::endl;
  return 0;
}
