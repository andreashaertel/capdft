// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef EXAMPLES_CARTESIAN_FUNCTIONALS_ALL_TOOLS_SRC_EXTPOT_HPP_
#define EXAMPLES_CARTESIAN_FUNCTIONALS_ALL_TOOLS_SRC_EXTPOT_HPP_
/** \file examples/cartesian_functionals_all_tools/src/extpot.hpp
 *  \brief Header file that contains functions to calculate external potentials.
 */
// Includes
#include <chrono>
#include <cmath>
#include <vector>
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
// Own namespace for external potential tool kit
namespace extpot {
// _____________________________________________________________________________
/** \brief Set external hard potential: two plates with sine wave shape
 *
 *  The left plate is wave shaped in the x direction, while the right plate
 *  is wave shaped in the y direction.
 */
void hard_wavelike(
    Properties& system_properties,
    std::vector<Properties>& species_properties,
    std::vector<size_t>& affected_species,
    std::vector<DataFrame<3, double>>* exp_ext_potential) {
  std::vector<double> system_lengths{0., 0., 0.};
  std::vector<size_t> grid_counts{0, 0, 0};
  system_properties.get_property("length x", &system_lengths.at(0));
  system_properties.get_property("length y", &system_lengths.at(1));
  system_properties.get_property("length z", &system_lengths.at(2));
  system_properties.get_property("grid count x", &grid_counts.at(0));
  system_properties.get_property("grid count y", &grid_counts.at(1));
  system_properties.get_property("grid count z", &grid_counts.at(2));
  double dx{system_lengths.at(0) / static_cast<double>(grid_counts.at(0))};
  double dy{system_lengths.at(1) / static_cast<double>(grid_counts.at(1))};
  double dz{system_lengths.at(2) / static_cast<double>(grid_counts.at(2))};
  double frequency{4. * M_PI / system_lengths.at(0)};
  double wave_left{0.}, wave_right{0.};;
  double x{0.}, y{0.}, z{0.};
  double diameter{0.};
  // exp(Hard wall external potential)
  for (auto& species : affected_species) {
    species_properties.at(species).get_property("diameter", &diameter);
    for (size_t i = 0; i != grid_counts.at(0); ++i) {
      x = dx * static_cast<double>(i);
      for (size_t j = 0; j != grid_counts.at(1); ++j) {
        y = dy * static_cast<double>(j);
        for (size_t k = 0; k != grid_counts.at(2); ++k) {
          z = dz * static_cast<double>(k);
          wave_left = .5 + 2. * pow(sin(frequency * x), 2);
          wave_right =.5 + 2. * pow(sin(frequency * y), 2);
          if (z < diameter * wave_left) {
            exp_ext_potential->at(species).at(i, j, k) *= 0.;
          } else if ((system_lengths.at(2) - z) < diameter * wave_right) {
            exp_ext_potential->at(species).at(i, j, k) *= 0.;
          } else {
            exp_ext_potential->at(species).at(i, j, k) *= 1.;
          }
        }
      }
    }
  }
}
// _____________________________________________________________________________
/** \brief Set external electrostatic potential: linear potenial in z-direction
 *  
 *  This potential works well with planar hard walls.
 */
void electrostatic_planar(
    Properties& system_properties,
    std::vector<Properties>& species_properties,
    std::vector<size_t>& affected_species,
    std::vector<DataFrame<3, double>>* exp_ext_potential) {
  // Set external electrostatic potential (linear between two charged plates)
  double system_length_z{0.};
  double potential{0.};
  std::vector<size_t> grid_counts{0, 0, 0};
  system_properties.get_property("length z", &system_length_z);
  system_properties.get_property("grid count x", &grid_counts.at(0));
  system_properties.get_property("grid count y", &grid_counts.at(1));
  system_properties.get_property("grid count z", &grid_counts.at(2));
  system_properties.get_property("potential", &potential);
  double dz{system_length_z / static_cast<double>(grid_counts.at(2))};
  double z{0.};
  double valency{0.};
  for (auto& species : affected_species) {
    species_properties.at(species).get_property("valency", &valency);
    for (size_t i = 0; i != grid_counts.at(0); ++i) {
      for (size_t j = 0; j != grid_counts.at(1); ++j) {
        for (size_t k = 0; k != grid_counts.at(2); ++k) {
          z = dz * (static_cast<double>(k) + 0.5);
          exp_ext_potential->at(species).at(i, j, k) *=
              exp(valency * potential * (1. / 2. - z / system_length_z));
        }
      }
    }
  }
}
// _____________________________________________________________________________
///** \brief Set external electrostatic potential: wavelike
// *  
// *  Set the external electrostatic potential in the following way.
// *  Due to the wave-like hard walls we also have wave like boundary conditions.
// *  We set these boundary conditions for a system without free charges and
// *  calculate the potential, which we then impose as external potential.
// */
//void electrostatic_wavelike(
//    Properties& system_properties,
//    std::vector<Properties>& species_properties,
//    std::vector<size_t>& affected_species,
//    std::vector<DataFrame<3, double>>* exp_ext_potential) {
//  // Extract the needed system properties
//  double potential{0.};
//  std::vector<double> system_lengths{0., 0., 0.};
//  std::vector<size_t> grid_counts{0, 0, 0};
//  system_properties.get_property("length x", &system_lengths.at(0));
//  system_properties.get_property("length y", &system_lengths.at(1));
//  system_properties.get_property("length z", &system_lengths.at(2));
//  system_properties.get_property("grid count x", &grid_counts.at(0));
//  system_properties.get_property("grid count y", &grid_counts.at(1));
//  system_properties.get_property("grid count z", &grid_counts.at(2));
//  system_properties.get_property("potential", &potential);
//  std::vector<double> bin_sizes{
//      system_lengths.at(0) / static_cast<double>(grid_counts.at(0)),
//      system_lengths.at(1) / static_cast<double>(grid_counts.at(1)),
//      system_lengths.at(2) / static_cast<double>(grid_counts.at(2))};
//  std::vector<bool> periodic_boundaries{true, true, false};
//  size_t voxel_count{1};
//  for (auto& gc : grid_counts) { voxel_count *= gc; }
//  // Create the wavelike boundary conditions
//  std::vector<std::vector<std::vector<std::vector<size_t>>>> boundary_pos(3);
//  std::vector<std::vector<std::vector<double>>> boundary_vals(3);
//  std::vector<size_t> pos_left, pos_right;
//  size_t k_left{0}, k_right{0};;
//  double x{0.}, y{0.};
//  double wave_left{0.}, wave_right{0.};;
//  double frequency{4. * M_PI / system_lengths.at(0)};
//  size_t zig = 8;
//  double k_left_base =
//      static_cast<size_t>(.5 * .3 / bin_sizes.at(2) + .5);
//  double k_right_base =
//      static_cast<size_t>((system_lengths.at(2) - .5 * .3) / bin_sizes.at(2) + .5);
//  for (size_t i = 0; i < grid_counts.at(0); ++i) {
//    for (size_t j = 0; j < grid_counts.at(1); ++j) {
//      // Produces zig-zag boundary conditions
//      if ((j / zig) % 2 == 0) { k_left = k_left_base + (j%zig) + 1; }
//      else { k_left = k_left_base + zig - j%zig - 1; }
//      if ((i / zig) % 2 == 0) { k_right = k_right_base - ((i%zig) + 1); }
//      else { k_right = k_right_base - (zig - i%zig - 1); }
//      pos_left = {i, j, k_left};
//      pos_right = {i, j, k_right};
//      boundary_pos.at(2).push_back({pos_left, pos_right});
//      boundary_vals.at(2).push_back({potential, -potential});
//    }
//  }
//  // Initial guess
//  std::vector<double> solution(voxel_count, 0.);
//  size_t index{0};
//  double slope{0.};
//  slope = -2. * potential / static_cast<double>(k_right_base-k_left_base);
//  for (size_t i = 0; i < grid_counts.at(0); ++i) {
//    for (size_t j = 0; j < grid_counts.at(1); ++j) {
//      for (size_t k = 0; k < grid_counts.at(2); ++k) {
//        index = i * grid_counts.at(2) * grid_counts.at(1) +
//            j * grid_counts.at(2) + k;
//        solution.at(index) = potential + slope * static_cast<double>(k-k_left_base);
//      }
//    }
//  }
//  // Solve the Poisson equation for every point in the system
//  CartesianPoissonSolverAny poisson_solver(
//      grid_counts, bin_sizes, periodic_boundaries, boundary_pos, boundary_vals);
//  std::vector<double> rhs(voxel_count, 0.);
//  poisson_solver.solve(rhs, solution);
//}
// _____________________________________________________________________________
}
#endif  // EXAMPLES_CARTESIAN_FUNCTIONALS_ALL_TOOLS_SRC_EXTPOT_HPP_
