// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_cartesian.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "data_frame.hpp"  // NOLINT
#include "integration.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalFMTCartesian::FunctionalFMTCartesian() {
  //
}
// _____________________________________________________________________________
FunctionalFMTCartesian::FunctionalFMTCartesian(
    const std::vector<DataFrame<3, double>>* density_profiles,
    const std::vector<Properties>& species_properties,
    const Properties& system_properties,
    const std::vector<size_t>& affected_species)
  : affected_species(affected_species),
    density_profiles_pointer(density_profiles) {
  // Clear all std::vectors
  diameters.clear();
  bulk_densities.clear();
  // Get system properties
  extract_system_properties(system_properties);
  // Get species properties; excludes all species without diameter property
  extract_species_properties(species_properties);
  // Initialize all data frames
  initialize_all_data_frames();
  update_density_profiles();
  // Calculate weights
  calc_weights();
}
// _____________________________________________________________________________
FunctionalFMTCartesian::FunctionalFMTCartesian(
      const std::vector<DataFrame<3, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties)
  : FunctionalFMTCartesian(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalFMTCartesian::~FunctionalFMTCartesian() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::extract_system_properties(
    const Properties& system_properties) {
  // Extract properties
  lengths = std::vector<double>(3);
  grid_counts = std::vector<size_t>(3);
  system_properties.get_property("length x", &lengths.at(0));
  system_properties.get_property("length y", &lengths.at(1));
  system_properties.get_property("length z", &lengths.at(2));
  system_properties.get_property("grid count x", &grid_counts.at(0));
  system_properties.get_property("grid count y", &grid_counts.at(1));
  system_properties.get_property("grid count z", &grid_counts.at(2));
  voxel_count = grid_counts.at(0) * grid_counts.at(1) * grid_counts.at(2);
  // Calculate bin sizes
  for (size_t i = 0; i < lengths.size(); ++i) {
    bin_sizes.push_back(lengths.at(i) / static_cast<double>(grid_counts.at(i)));
    bin_sizes_four.push_back(2. * M_PI / lengths.at(i));
  }
  // Since we only do real-->complex Fourier transforms (or vice versa),
  // we only need half of the Fourier space (see FFTW docs).
  // In 3D that means that the first two dimensions have the regular number of
  // grid points and the last one has half the grid points.
  grid_counts_four.push_back(grid_counts.at(0));
  grid_counts_four.push_back(grid_counts.at(1));
  grid_counts_four.push_back(grid_counts.at(2) / 2 + 1);
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::extract_species_properties(
    const std::vector<Properties>& species_properties) {
  // Sort the affected species numbers
  std::sort(affected_species.begin(), affected_species.end());
  // Remove duplicates
  affected_species.erase(
      unique(affected_species.begin(), affected_species.end()),
      affected_species.end());
  // If no affected species were specified, find them automatically
  const std::vector<Properties>& spec_prop = species_properties;
  double diameter{0.};
  double bulk_density{0.};
  if (affected_species.empty()) {
    for (auto it = spec_prop.begin(); it != spec_prop.end(); ++it) {
      if (it->get_property("diameter", &diameter)) {  // only species with diam.
        if (!it->get_property("bulk density", &bulk_density)) {
          std::cerr << "FunctionalFMTCartesian::extract_species_properties(): ";
          std::cerr << "\"Error: A species with a diameter but no density ";
          std::cerr << "was detected.\"" << std::endl;
          exit(1);
        }
        affected_species.push_back(it - spec_prop.begin());
      }
    }
  }
  // Extract properties
  for (auto& species : affected_species) {
    if (!spec_prop.at(species).get_property("diameter", &diameter) ||
        !spec_prop.at(species).get_property("bulk density", &bulk_density)) {
      std::cerr << "FunctionalFMTCartesian::extract_species_properties(): ";
      std::cerr << "\"Error: One species is lacking a required parameter.";
      std::cerr << std::endl;
      exit(1);
    }
    diameters.push_back(diameter);
    bulk_densities.push_back(bulk_density);
  }
  // Count species that interact via the hard sphere potential
  species_count = affected_species.size();
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::initialize_all_data_frames() {
  // Initialize internal density profiles
  for (size_t i = 0; i < species_count; ++i) {
    density_profiles.push_back(DataFrame<3, double>(grid_counts));
    density_profiles_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
  // Initialize weights
  for (size_t i = 0; i < species_count; ++i) {
    scalar_weights_four.push_back(std::vector<DataFrame<3, fftw_complex>>(0));
    vector_weights_four.push_back(std::vector<DataFrame<3, fftw_complex>>(0));
    tensor_weights_four.push_back(std::vector<DataFrame<3, fftw_complex>>(0));
    for (size_t j = 0; j < 4; ++j) {
      scalar_weights_four.at(i).push_back(
          DataFrame<3, fftw_complex>(grid_counts_four));
    }
    for (size_t j = 0; j < 6; ++j) {
      vector_weights_four.at(i).push_back(
          DataFrame<3, fftw_complex>(grid_counts_four));
      tensor_weights_four.at(i).push_back(
          DataFrame<3, fftw_complex>(grid_counts_four));
    }
  }
  // Initialize weighted densities
  for (size_t i = 0; i < 4; ++i) {
    scalar_weighted_dens_real.push_back(DataFrame<3, double>(grid_counts));
    scalar_weighted_dens_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
  for (size_t i = 0; i < 6; ++i) {
    vector_weighted_dens_real.push_back(DataFrame<3, double>(grid_counts));
    tensor_weighted_dens_real.push_back(DataFrame<3, double>(grid_counts));
    vector_weighted_dens_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
    tensor_weighted_dens_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
  // Initialize partial derivatives
  for (size_t i = 0; i < 4; ++i) {
    scalar_partial_derivative_real.push_back(DataFrame<3, double>(grid_counts));
    scalar_partial_derivative_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
  for (size_t i = 0; i < 6; ++i) {
    vector_partial_derivative_real.push_back(DataFrame<3, double>(grid_counts));
    tensor_partial_derivative_real.push_back(DataFrame<3, double>(grid_counts));
    vector_partial_derivative_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
    tensor_partial_derivative_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
  // Initialize internal functional derivatives (Fourier transform)
  for (size_t i = 0; i < species_count; ++i) {
    functional_derivative_four.push_back(
        DataFrame<3, fftw_complex>(grid_counts_four));
  }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::update_density_profiles() {
  size_t index{0};
  for (size_t i = 0; i < species_count; ++i) {
    index = affected_species.at(i);
    density_profiles.at(i) = density_profiles_pointer->at(index);
  }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_weights() {
  // Declare auxiliary variables. R means radius of the hard sphere.
  // k is the absolute value of the position in Fourier space.
  // RR is the square of R, kk is the square of k, and so on.
  double kx{0.}, ky{0.}, kz{0.};
  double k{0.}, kk{0.}, kkk{0.};
  double R{0.}, RR{0.}, RRR{0.};
  double Rk{0.};
  set_weights_to_zero();
  for (size_t i = 0; i != species_count; ++i) {
    R = diameters.at(i) / 2.;
    RR = R * R;
    RRR = RR * R;
    for (size_t l = 0; l != grid_counts_four.at(0); ++l) {
      for (size_t m = 0; m != grid_counts_four.at(1); ++m) {
        for (size_t n = 0; n != grid_counts_four.at(2); ++n) {
          if (l < grid_counts_four.at(0) / 2 + 1) {  // stupid convention
            kx = bin_sizes_four.at(0) * static_cast<double>(l);
          } else {
            kx = -bin_sizes_four.at(0) *
                static_cast<double>(grid_counts_four.at(0)-l);
          }
          if (m < grid_counts_four.at(1) / 2 + 1) {  // stupid convention
            ky = bin_sizes_four.at(1) * static_cast<double>(m);
          } else {
            ky = -bin_sizes_four.at(1) *
                static_cast<double>(grid_counts_four.at(1)-m);
          }
          kz = bin_sizes_four.at(2) * static_cast<double>(n);
          kk = kx * kx + ky * ky + kz * kz;
          k = sqrt(kk);
          kkk = k * kk;
          Rk = R * k;
          // Calculate limits k --> 0
          if (l == 0 && m == 0 && n == 0) {
            // scalar weight 3 (theta)
            scalar_weights_four.at(i).at(0).at(l, m, n)[0] = 4.*M_PI * RRR / 3.;
            // scalar weight 2 (delta)
            scalar_weights_four.at(i).at(1).at(l, m, n)[0] = 4. * M_PI * RR;
            // scalar weight 1 (delta)
            scalar_weights_four.at(i).at(2).at(l, m, n)[0] = R;
            // scalar weight 0 (delta)
            scalar_weights_four.at(i).at(3).at(l, m, n)[0] = 1.;
          } else {
            // scalar weight 3 (theta)
            scalar_weights_four.at(i).at(0).at(l, m, n)[0] = 4. * M_PI *
                (sin(Rk) - Rk * cos(Rk)) / kkk;  // real
            // scalar weight 2 (delta)
            scalar_weights_four.at(i).at(1).at(l, m, n)[0] = 4. * M_PI *
                R * sin(Rk) / k;
            // scalar weight 1 (delta)
            scalar_weights_four.at(i).at(2).at(l, m, n)[0] = sin(Rk) / k;
            // scalar weight 0 (delta)
            scalar_weights_four.at(i).at(3).at(l, m, n)[0] = sin(Rk) / Rk;
            // vectorial weight 2 1st element
            vector_weights_four.at(i).at(0).at(l, m, n)[1] =
                -kx * scalar_weights_four.at(i).at(0).at(l, m, n)[0];
            // vectorial weight 2 2nd element
            vector_weights_four.at(i).at(1).at(l, m, n)[1] =
                -ky * scalar_weights_four.at(i).at(0).at(l, m, n)[0];
            // vectorial weight 2 3rd element
            vector_weights_four.at(i).at(2).at(l, m, n)[1] =
                -kz * scalar_weights_four.at(i).at(0).at(l, m, n)[0];
            // vectorial weight 1 1st element
            vector_weights_four.at(i).at(3).at(l, m, n)[1] =
                vector_weights_four.at(i).at(0).at(l, m, n)[1] / (4.*M_PI * R);
            // vectorial weight 1 2nd element
            vector_weights_four.at(i).at(4).at(l, m, n)[1] =
                vector_weights_four.at(i).at(1).at(l, m, n)[1] / (4.*M_PI * R);
            // vectorial weight 1 3rd element
            vector_weights_four.at(i).at(5).at(l, m, n)[1] =
                vector_weights_four.at(i).at(2).at(l, m, n)[1] / (4.*M_PI * R);
            // tensorial weight (1,1)-element
            tensor_weights_four.at(i).at(0).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (kx * kx / kk - 1. / 3.);
            // tensorial weight (1,2)-element
            tensor_weights_four.at(i).at(1).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (kx * ky / kk);
            // tensorial weight (1,3)-element
            tensor_weights_four.at(i).at(2).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (kx * kz / kk);
            // tensorial weight (2,2)-element
            tensor_weights_four.at(i).at(3).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (ky * ky / kk - 1. / 3.);
            // tensorial weight (2,3)-element
            tensor_weights_four.at(i).at(4).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (ky * kz / kk);
            // tensorial weight (3,3)-element
            tensor_weights_four.at(i).at(5).at(l, m, n)[0] =
                (scalar_weights_four.at(i).at(1).at(l, m, n)[0] -
                3. * scalar_weights_four.at(i).at(0).at(l, m, n)[0] / R) *
                (kz * kz / kk - 1. / 3.);
          }  // else
        }  // for n (kz)
      }  // for m (ky)
    }  // for l (kx)
  }  // for i (species)
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::set_weights_to_zero() {
  for (size_t i = 0; i != species_count; ++i) {
    for (size_t j = 0; j != scalar_weights_four.at(i).size(); ++j) {
      scalar_weights_four.at(i).at(j).zero();
    }
    for (size_t j = 0; j != vector_weights_four.at(i).size(); ++j) {
      vector_weights_four.at(i).at(j).zero();
    }
    for (size_t j = 0; j != tensor_weights_four.at(i).size(); ++j) {
      tensor_weights_four.at(i).at(j).zero();
    }
  }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != voxel_count) {
      std::cerr << "FunctionalFMTCartesian::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  // Calculate the weighted densities
  calc_weighted_densities();
  // From the weighted densities calculate the partial derivatives of the
  // excess free energy
  calc_partial_derivatives();
  // Calculate the weighted partial derivatives
  calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_derivative_warnings(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != voxel_count) {
      std::cerr << "FunctionalFMTCartesian::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  // Calculate the weighted densities
  calc_weighted_densities();
  // Check whether there are unphysical values in the  weighted densities
  check_weighted_densities();
  // From the weighted densities calculate the partial derivatives of the
  // excess free energy
  calc_partial_derivatives();
  // Calculate the weighted partial derivatives
  calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  // Auxiliary variables
  double diameter{0.}, diameter2{0.}, diameter3{0.};
  // Calculate the fluid weighted densities and save them to the first position
  for (auto& weighted_density : scalar_weighted_dens_real) {
    weighted_density.at(0) = 0.;
  }
  for (size_t i = 0; i != species_count; ++i) {
    // Define auxiliary values
    diameter = diameters.at(i);
    diameter2 = diameter * diameter;
    diameter3 = diameter2 * diameter;
    // All non-scalar weighted densities are 0 the others are obviously constant
    scalar_weighted_dens_real.at(0).at(0) +=
        bulk_densities.at(i) * diameter3 * M_PI / 6.0;
    scalar_weighted_dens_real.at(1).at(0) +=
        bulk_densities.at(i) * diameter2 * M_PI;
    scalar_weighted_dens_real.at(2).at(0) +=
        bulk_densities.at(i) * diameter / 2.0;
    scalar_weighted_dens_real.at(3).at(0) += bulk_densities.at(i);
  }
  // Evaluate the partial derivatives at the first position
  calc_local_partial_derivatives(0);
  std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
  size_t spec_i{0};
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    spec_i = it - affected_species.begin();
    // Define auxiliary values
    diameter = diameters.at(spec_i);
    diameter2 = diameter * diameter;
    diameter3 = diameter2 * diameter;
    // Calculate the derivatives of the excess free energy functional w.r.t. the
    // weighted densities.
    bulk_derivative->at(*it) = scalar_partial_derivative_real.at(0).at(0) *
        diameter3 * M_PI / 6.0;
    bulk_derivative->at(*it) += scalar_partial_derivative_real.at(1).at(0) *
        diameter2 * M_PI;
    bulk_derivative->at(*it) += scalar_partial_derivative_real.at(2).at(0) *
        diameter / 2.0;
    bulk_derivative->at(*it) += scalar_partial_derivative_real.at(3).at(0);
  }
}
// _____________________________________________________________________________
double FunctionalFMTCartesian::calc_energy() {
  double integral{0.};
  DataFrame<3, double> free_energy_density(grid_counts);
  // Calculate the weighted densities for the current density profile
  calc_weighted_densities();
  check_weighted_densities();
  // Store free energy density
  for (size_t i = 0; i < voxel_count; ++i) {
    free_energy_density.at(i) = calc_local_energy_density(i);
  }
  integral = integration_3d_closed(free_energy_density, bin_sizes);
  return integral;
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_weighted_densities() {
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
  // Specify the plans
  // (doing this locally produces just little numerical overhead)
  for (size_t i = 0; i != species_count; ++i) {
    forward_plans.push_back(
        fftw_plan_dft_r2c_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            density_profiles.at(i).array(),  // input
            density_profiles_four.at(i).array(),  // output
            flags_keep));
  }
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            scalar_weighted_dens_four.at(i).array(),
            scalar_weighted_dens_real.at(i).array(),
            flags_destroy));
  }
  for (size_t i = 0; i != vector_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            vector_weighted_dens_four.at(i).array(),
            vector_weighted_dens_real.at(i).array(),
            flags_destroy));
  }
  for (size_t i = 0; i != tensor_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            tensor_weighted_dens_four.at(i).array(),
            tensor_weighted_dens_real.at(i).array(),
            flags_destroy));
  }
  // Creating plans clears the input array hence we need to update it
  update_density_profiles();
  // Forward Fourier transform of density profiles
  for (auto& plan : forward_plans) {
    fftw_execute(plan);
  }
  // Normalize
  for(auto& profile_four : density_profiles_four) {  // every species
    profile_four *= bin_sizes.at(0) * bin_sizes.at(1) * bin_sizes.at(2);
  }
  // Convolution
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    scalar_weighted_dens_four.at(i).zero();
    for (size_t j = 0; j != species_count; ++j) {
      scalar_weighted_dens_four.at(i) +=
          density_profiles_four.at(j) * scalar_weights_four.at(j).at(i);
    }
  }
  for (size_t i = 0; i != vector_weighted_dens_four.size(); ++i) {
    vector_weighted_dens_four.at(i).zero();
    for (size_t j = 0; j != species_count; ++j) {
      vector_weighted_dens_four.at(i) +=
          (density_profiles_four.at(j) * vector_weights_four.at(j).at(i));
    }
  }
  for (size_t i = 0; i != tensor_weighted_dens_four.size(); ++i) {
    tensor_weighted_dens_four.at(i).zero();
    for (size_t j = 0; j != species_count; ++j) {
      tensor_weighted_dens_four.at(i) +=
          density_profiles_four.at(j) * tensor_weights_four.at(j).at(i);
    }
  }
  // Forward Fourier transform of density profiles
  for (auto& plan : backward_plans) {
    fftw_execute(plan);
  }
  // Normalize
  for(auto& profile_real : scalar_weighted_dens_real) {
    profile_real *= bin_sizes_four.at(0) * bin_sizes_four.at(1) *
        bin_sizes_four.at(2) / pow(2. * M_PI, 3);
  }
  for(auto& profile_real : vector_weighted_dens_real) {
    profile_real *= bin_sizes_four.at(0) * bin_sizes_four.at(1) *
        bin_sizes_four.at(2) / pow(2. * M_PI, 3);
  }
  for(auto& profile_real : tensor_weighted_dens_real) {
    profile_real *= bin_sizes_four.at(0) * bin_sizes_four.at(1) *
        bin_sizes_four.at(2) / pow(2. * M_PI, 3);
  }
  // The convolution may have numerical at both ends of the array.
  // One way of avoiding this is by using an outer external potential
  // forcing the density profile to vanish (i.e. two planar walls).
  for (auto& plan : forward_plans) { fftw_destroy_plan(plan); }
  for (auto& plan : backward_plans) { fftw_destroy_plan(plan); }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::check_weighted_densities() {
  // initialize counters
  size_t n3bigger1_counter = 0;
  size_t n3smaller0_counter = 0;
  size_t n2smaller0_counter = 0;
  // Since calc_weighted_densities() is executed before this function,
  // we can check if the weighted densities are in the legal range.
  for (size_t i = 0; i < voxel_count; ++i) {
    // Check scalar weighted densities for values smaller 0
    if (scalar_weighted_dens_real.at(0).at(i) < 0.) {
      ++n3smaller0_counter;
      scalar_weighted_dens_real.at(0).at(i) = 0.;
    }
    if (scalar_weighted_dens_real.at(1).at(i) < 0.) {
      ++n2smaller0_counter;
      scalar_weighted_dens_real.at(1).at(i) = 0.;
      scalar_weighted_dens_real.at(2).at(i) = 0.;
      scalar_weighted_dens_real.at(3).at(i) = 0.;
    }
    // Check if n3 is bigger than 1
    if (scalar_weighted_dens_real.at(0).at(i) > 1.) {
      ++n3bigger1_counter;
      scalar_weighted_dens_real.at(0).at(i) = 1.;
    }
  }
  // Show warning
  if (n3smaller0_counter > 0) {
    std::cerr << "FunctionalFMTCartesian::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n3 < 0.0 at " << n3smaller0_counter;
    std::cerr << " positions.\"" << std::endl;
  }
  if (n2smaller0_counter > 0) {
    std::cerr << "FunctionalFMTCartesian::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n2 < 0.0 at " << n2smaller0_counter;
    std::cerr << " positions.\"" << std::endl;
  }
  if (n3bigger1_counter > 0) {
    std::cerr << "FunctionalFMTCartesian::check_weighted_densities(): \"";
    std::cerr << "Warning: Local n3 > 1.0 at " << n3bigger1_counter;
    std::cerr << " positions.\"" << std::endl;
  }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_partial_derivatives() {
  // Calculate the partial deriavative at every position
  for (size_t i = 0; i != voxel_count; ++i) {
    calc_local_partial_derivatives(i);
  }
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_local_partial_derivatives(size_t i) {
  // Auxiliary variables for the weight functions
  // (i am truly sorry ...)
  // Skalar weighted densities
  double n3{scalar_weighted_dens_real.at(0).at(i)};
  double n2{scalar_weighted_dens_real.at(1).at(i)};
  double n1{scalar_weighted_dens_real.at(2).at(i)};
  double n0{scalar_weighted_dens_real.at(3).at(i)};
  double n3n3{n3 * n3};
  double n3n3n3{n3 * n3n3};
  double oneMn3{1. - n3};  // ONEMinusN3
  double oneOn3{1. / n3};  // ONEOverN3
  double oneOn3n3{1. / n3n3};
  double oneOoneMn3{1. / oneMn3};
  double n2n2{n2 * n2};
  double n2n2n2{n2 * n2n2};
  // Vector weighted densities
  double nvec2x{vector_weighted_dens_real.at(0).at(i)};
  double nvec2y{vector_weighted_dens_real.at(1).at(i)};
  double nvec2z{vector_weighted_dens_real.at(2).at(i)};
  double nvec1x{vector_weighted_dens_real.at(3).at(i)};
  double nvec1y{vector_weighted_dens_real.at(4).at(i)};
  double nvec1z{vector_weighted_dens_real.at(5).at(i)};
  double nvec1nvec2{nvec1x * nvec2x + nvec1y * nvec2y + nvec1z * nvec2z};
  double nvec2nvec2{nvec2x * nvec2x + nvec2y * nvec2y + nvec2z * nvec2z};
  // Tensor weighted densities (symmetric)
  double ntensorxx{tensor_weighted_dens_real.at(0).at(i)};
  double ntensorxy{tensor_weighted_dens_real.at(1).at(i)};
  double ntensorxz{tensor_weighted_dens_real.at(2).at(i)};
  double ntensoryy{tensor_weighted_dens_real.at(3).at(i)};
  double ntensoryz{tensor_weighted_dens_real.at(4).at(i)};
  double ntensorzz{tensor_weighted_dens_real.at(5).at(i)};
  double nvec2ntensornvec2 =  // nvec2 * ntensor * nvec2
      nvec2x * nvec2x * ntensorxx + nvec2y * nvec2y * ntensoryy +
      nvec2z * nvec2z * ntensorzz + 2. * nvec2x * nvec2y * ntensorxy +
      2. * nvec2x * nvec2z * ntensorxz + 2. * nvec2y * nvec2z * ntensoryz;
  double trace3 =  // Tr(ntensor * ntensor * ntensor)
      pow(ntensorxx, 3) + pow(ntensoryy, 3) + pow(ntensorzz, 3) +
      3. * pow(ntensorxy, 2) * ntensorxx + 3. * pow(ntensorxy, 2) * ntensoryy +
      3. * pow(ntensorxz, 2) * ntensorxx + 3. * pow(ntensorxz, 2) * ntensorzz +
      3. * pow(ntensoryz, 2) * ntensoryy + 3. * pow(ntensoryz, 2) * ntensorzz +
      6. * ntensorxy * ntensorxz * ntensoryz;
  double ntensornvec2x =  // ntensor * nvec2
      ntensorxx * nvec2x + ntensorxy * nvec2y + ntensorxz * nvec2z;
  double ntensornvec2y =
      ntensorxy * nvec2x + ntensoryy * nvec2y + ntensoryz * nvec2z;
  double ntensornvec2z =
      ntensorxz * nvec2x + ntensoryz * nvec2y + ntensorzz * nvec2z;

  double ntensorntensorxx =
      pow(ntensorxx, 2) + pow(ntensorxy, 2) + pow(ntensorxz, 2);
  double ntensorntensorxy =
      ntensorxx * ntensorxy + ntensorxy * ntensoryy + ntensorxz * ntensoryz;
  double ntensorntensorxz =
      ntensorxx * ntensorxz + ntensorxy * ntensoryz + ntensorxz * ntensorzz;
  double ntensorntensoryy =
      pow(ntensorxy, 2) + pow(ntensoryy, 2) + pow(ntensoryz, 2);
  double ntensorntensoryz =
      ntensorxy * ntensorxz + ntensoryy * ntensoryz + ntensoryz * ntensorzz;
  double ntensorntensorzz =
      pow(ntensorxz, 2) + pow(ntensoryz, 2) + pow(ntensorzz, 2);
  // Auxiliary variables \phi^{num}_2 and \phi^{num}_3 (B.18)-(B.21)
  double phi2, phi3;
  // Auxiliary variables \partial\phi^{num}_j/partial n_3 (B.30)-(B.33)
  double dphi2dn3, dphi3dn3;
  // Auxiliary variables \partial\Phi_j/\partial n_k (B.37)-(B.46)
  double dPhi1dn0, dPhi1dn3;
  double dPhi2dn3, dPhi2dn2, dPhi2dn1;
  double dPhi2dnvec2x, dPhi2dnvec2y, dPhi2dnvec2z;
  double dPhi2dnvec1x, dPhi2dnvec1y, dPhi2dnvec1z;
  double dPhi3dn2, dPhi3dn3;
  double dPhi3dnvec2x, dPhi3dnvec2y, dPhi3dnvec2z;
  double dPhi3dntensorxx, dPhi3dntensoryy, dPhi3dntensorzz;
  double dPhi3dntensorxy, dPhi3dntensorxz, dPhi3dntensoryz;
  // Auxiliary constants
  double oneO24pi{1. / (24. * M_PI)};
  // Prefactors that depend on n3. They contain logarithm terms, which need to
  // be approximated (Taylor) in certain situations.
  if (n3 >= sqrt(std::numeric_limits<double>::epsilon())) {
    phi2 = (5./3.) + (2./3.) * (oneMn3*oneOn3) * log(oneMn3) - (n3/3.);
    phi3 = 2. - (2./3.) *
        (oneOn3 + n3 + (oneMn3 * oneMn3 * oneOn3n3) * log(oneMn3));
    dphi2dn3 = -(2./3.) * (oneOn3n3 * log(oneMn3) + oneOn3 + .5);
    dphi3dn3 = (2./3.) * oneOn3n3 *
        (2. - n3 - n3n3 + 2. * (oneMn3*oneOn3) * log(oneMn3));
  } else {
    phi2 = 1. + (1./9.) * n3n3 + (1./18.) * n3n3n3 + (1./30.) * n3n3n3 * n3 +
        (1./45.) * n3n3n3 * n3n3;
    phi3 = 1. - (4./9.) * n3 + (1./18.) * n3n3 + (1./45.) * n3n3n3 +
        (1./90.) * n3n3n3 * n3;
    dphi2dn3 = (2./9.) * n3 + (1./6.) * n3n3 + (2./15) * n3n3n3 +
        (1./9.) * n3n3n3 * n3;
    dphi3dn3 = -(4./9.) + (1./9.) * n3 + (1./15.) * n3n3 + (2./45.) * n3n3n3;
  }
  // Calculate the partial derivatives \partial\Phi_j/\partial n_k (B.37)-(B.46)
  // Phi1 derivatives
  dPhi1dn3 = n0 * oneOoneMn3;
  dPhi1dn0 = -log(oneMn3);
  // Phi2 derivatives
  dPhi2dn3 = (phi2 * oneOoneMn3 + dphi2dn3) * (n2 * n1 - nvec1nvec2) *
      oneOoneMn3;
  dPhi2dn2 = phi2 * n1 * oneOoneMn3;
  dPhi2dn1 = phi2 * n2 * oneOoneMn3;
  dPhi2dnvec2x = -phi2 * nvec1x * oneOoneMn3;
  dPhi2dnvec2y = -phi2 * nvec1y * oneOoneMn3;
  dPhi2dnvec2z = -phi2 * nvec1z * oneOoneMn3;
  dPhi2dnvec1x = -phi2 * nvec2x * oneOoneMn3;
  dPhi2dnvec1y = -phi2 * nvec2y * oneOoneMn3;
  dPhi2dnvec1z = -phi2 * nvec2z * oneOoneMn3;
  // Phi3 derivatives
  dPhi3dn2 = phi3 * (3. * n2n2 - 3. * nvec2nvec2) * oneO24pi * oneOoneMn3 *
      oneOoneMn3;
  dPhi3dn3 = (2. * phi3 * oneOoneMn3 + dphi3dn3) *
      (n2n2n2 - 3. * n2 * nvec2nvec2 +
      4.5 * (nvec2ntensornvec2 - trace3)) * oneO24pi * oneOoneMn3 *
      oneOoneMn3;
  dPhi3dnvec2x = phi3 * (9. * ntensornvec2x - 6. * n2 * nvec2x) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dnvec2y = phi3 * (9. * ntensornvec2y - 6. * n2 * nvec2y) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dnvec2z = phi3 * (9. * ntensornvec2z - 6. * n2 * nvec2z) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensorxx = 4.5 * phi3 * (nvec2x * nvec2x - 3. * ntensorntensorxx) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensorxy = 4.5 * phi3 * (nvec2x * nvec2y - 3. * ntensorntensorxy) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensorxz = 4.5 * phi3 * (nvec2x * nvec2z - 3. * ntensorntensorxz) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensoryy = 4.5 * phi3 * (nvec2y * nvec2y - 3. * ntensorntensoryy) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensoryz = 4.5 * phi3 * (nvec2y * nvec2z - 3. * ntensorntensoryz) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  dPhi3dntensorzz = 4.5 * phi3 * (nvec2z * nvec2z - 3. * ntensorntensorzz) *
      oneO24pi * oneOoneMn3 * oneOoneMn3;
  // Sum partial derivatives \sum_{i=1}^3\partial\Phi_i/\partial n_k
  // Scalar n3
  scalar_partial_derivative_real.at(0).at(i) = dPhi1dn3 + dPhi2dn3 + dPhi3dn3;
  // Scalar n2
  scalar_partial_derivative_real.at(1).at(i) = dPhi2dn2 + dPhi3dn2;
  // Scalar n1
  scalar_partial_derivative_real.at(2).at(i) = dPhi2dn1;
  // Scalar n0
  scalar_partial_derivative_real.at(3).at(i) = dPhi1dn0;
  // Vector n2
  vector_partial_derivative_real.at(0).at(i) = dPhi2dnvec2x + dPhi3dnvec2x;
  vector_partial_derivative_real.at(1).at(i) = dPhi2dnvec2y + dPhi3dnvec2y;
  vector_partial_derivative_real.at(2).at(i) = dPhi2dnvec2z + dPhi3dnvec2z;
  // Vector n1
  vector_partial_derivative_real.at(3).at(i) = dPhi2dnvec1x;
  vector_partial_derivative_real.at(4).at(i) = dPhi2dnvec1y;
  vector_partial_derivative_real.at(5).at(i) = dPhi2dnvec1z;
  // Tensor nm2 (first and second diagonal element are the same)
  tensor_partial_derivative_real.at(0).at(i) = dPhi3dntensorxx;
  tensor_partial_derivative_real.at(1).at(i) = dPhi3dntensorxy;
  tensor_partial_derivative_real.at(2).at(i) = dPhi3dntensorxz;
  tensor_partial_derivative_real.at(3).at(i) = dPhi3dntensoryy;
  tensor_partial_derivative_real.at(4).at(i) = dPhi3dntensoryz;
  tensor_partial_derivative_real.at(5).at(i) = dPhi3dntensorzz;
}
// _____________________________________________________________________________
void FunctionalFMTCartesian::calc_weighted_partial_derivatives(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  size_t spec_i{0};
  std::vector<double> sign_convolution_scalar{1., 1., 1., 1.};
  std::vector<double> sign_convolution_vector{-1., -1., -1., -1., -1., -1.};
  std::vector<double> sign_convolution_tensor{1., 2., 2., 1., 2., 1.};
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
  // Specify the plans
  for (size_t i = 0; i != scalar_partial_derivative_real.size(); ++i) {
    forward_plans.push_back(
        fftw_plan_dft_r2c_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            scalar_partial_derivative_real.at(i).array(),
            scalar_partial_derivative_four.at(i).array(),
            flags_keep));
  }
  for (size_t i = 0; i != vector_partial_derivative_real.size(); ++i) {
    forward_plans.push_back(
        fftw_plan_dft_r2c_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            vector_partial_derivative_real.at(i).array(),
            vector_partial_derivative_four.at(i).array(),
            flags_keep));
  }
  for (size_t i = 0; i != tensor_partial_derivative_real.size(); ++i) {
    forward_plans.push_back(
        fftw_plan_dft_r2c_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            tensor_partial_derivative_real.at(i).array(),
            tensor_partial_derivative_four.at(i).array(),
            flags_keep));
  }
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    spec_i = it - affected_species.begin();  // internal species number
    backward_plans.push_back(
        fftw_plan_dft_c2r_3d(
            grid_counts.at(0), grid_counts.at(1), grid_counts.at(2),
            functional_derivative_four.at(spec_i).array(),
            functional_derivative->at(*it).array(),
            flags_destroy));
  }
  // Transform the partial derivatives into Fourier space.
  for (auto& plan : forward_plans) {
    fftw_execute(plan);
  }
  // Normalize
  for (auto& partial_derivative : scalar_partial_derivative_four) {
      partial_derivative *= bin_sizes.at(0) * bin_sizes.at(1) * bin_sizes.at(2);
  }
  for (auto& partial_derivative : vector_partial_derivative_four) {
      partial_derivative *= bin_sizes.at(0) * bin_sizes.at(1) * bin_sizes.at(2);
  }
  for (auto& partial_derivative : tensor_partial_derivative_four) {
      partial_derivative *= bin_sizes.at(0) * bin_sizes.at(1) * bin_sizes.at(2);
  }
  // Convolution
  for (size_t i = 0; i < species_count; ++i) {
    functional_derivative_four.at(i).zero();
    // Scalar
    for (size_t j = 0; j < scalar_weights_four.at(0).size(); ++j) {
      functional_derivative_four.at(i) += sign_convolution_scalar.at(j) *
          scalar_partial_derivative_four.at(j) *
          scalar_weights_four.at(i).at(j);
    }
    // Vector (sign is negative, because of the assymetry of the vector weights)
    for (size_t j = 0; j < vector_weights_four.at(0).size(); ++j) {
      functional_derivative_four.at(i) += sign_convolution_vector.at(j) *
          vector_partial_derivative_four.at(j) *
          vector_weights_four.at(i).at(j);
    }
    // Tensor
    for (size_t j = 0; j < tensor_weights_four.at(0).size(); ++j) {
      functional_derivative_four.at(i) += sign_convolution_tensor.at(j) *
          tensor_partial_derivative_four.at(j) *
          tensor_weights_four.at(i).at(j);
    }
  }
  // Back transform
  for (auto& plan : backward_plans) {
    fftw_execute(plan);
  }
  // Normalize
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    functional_derivative->at(*it) *= bin_sizes_four.at(0) *
        bin_sizes_four.at(1) * bin_sizes_four.at(2) / pow(2. * M_PI, 3);
  }
  for (auto& plan : forward_plans) { fftw_destroy_plan(plan); }
  for (auto& plan : backward_plans) { fftw_destroy_plan(plan); }
}
// _____________________________________________________________________________
double FunctionalFMTCartesian::calc_local_energy_density(size_t i) {
  // Auxiliary variables for the weight functions
  // (i am truly sorry ...)
  // Skalar weighted densities
  double n3{scalar_weighted_dens_real.at(0).at(i)};
  double n2{scalar_weighted_dens_real.at(1).at(i)};
  double n1{scalar_weighted_dens_real.at(2).at(i)};
  double n0{scalar_weighted_dens_real.at(3).at(i)};
  double n3n3{n3 * n3};
  double n3n3n3{n3 * n3n3};
  double oneMn3{1. - n3};  // ONEMinusN3
  double oneMn3Squared{oneMn3 * oneMn3};
  double logoneMn3{log(oneMn3)};
  double n2n2n2{n2 * n2 * n2};
  // Vector weighted densities
  double nvec2x{vector_weighted_dens_real.at(0).at(i)};
  double nvec2y{vector_weighted_dens_real.at(1).at(i)};
  double nvec2z{vector_weighted_dens_real.at(2).at(i)};
  double nvec1x{vector_weighted_dens_real.at(3).at(i)};
  double nvec1y{vector_weighted_dens_real.at(4).at(i)};
  double nvec1z{vector_weighted_dens_real.at(5).at(i)};
  double nvec1nvec2{nvec1x * nvec2x + nvec1y * nvec2y + nvec1z * nvec2z};
  double nvec2nvec2{nvec2x * nvec2x + nvec2y * nvec2y + nvec2z * nvec2z};
  //double nvec2nvec2{nvec2x * nvec2x + nvec2y * nvec2y + nvec2z * nvec2z};
  // Tensor weighted densities (symmetric)
  double ntensorxx{tensor_weighted_dens_real.at(0).at(i)};
  double ntensorxy{tensor_weighted_dens_real.at(1).at(i)};
  double ntensorxz{tensor_weighted_dens_real.at(2).at(i)};
  double ntensoryy{tensor_weighted_dens_real.at(3).at(i)};
  double ntensoryz{tensor_weighted_dens_real.at(4).at(i)};
  double ntensorzz{tensor_weighted_dens_real.at(5).at(i)};
  double nvec2ntensornvec2 =  // nvec2 * ntensor * nvec2
      nvec2x * nvec2x * ntensorxx + nvec2y * nvec2y * ntensoryy +
      nvec2z * nvec2z * ntensorzz + 2. * nvec2x * nvec2y * ntensorxy +
      2. * nvec2x * nvec2z * ntensorxz + 2. * nvec2y * nvec2z * ntensoryz;
  double trace3 =  // Tr(ntensor * ntensor * ntensor)
      pow(ntensorxx, 3) + pow(ntensoryy, 3) + pow(ntensorzz, 3) +
      3. * pow(ntensorxy, 2) * ntensorxx + 3. * pow(ntensorxy, 2) * ntensoryy +
      3. * pow(ntensorxz, 2) * ntensorxx + 3. * pow(ntensorxz, 2) * ntensorzz +
      3. * pow(ntensoryz, 2) * ntensoryy + 3. * pow(ntensoryz, 2) * ntensorzz +
      6. * ntensorxy * ntensorxz * ntensoryz;
  // Auxiliary variables \phi^{num}_2 and \phi^{num}_3 (3.60),(3.61)
  double phi2, phi3;
  // Auxiliary variables \phi^{num}_2 and \phi^{num}_3 (B.18)-(B.21)
  double phi2num, phi3num;
  // Auxiliary variables \Phi_j (3.57),(3.58), (B.17)
  double Phi1, Phi2, Phi3;
  // Calculate factors in the energy density terms
  if (n3 < sqrt(std::numeric_limits<double>::epsilon())) {
    phi2 = 1. + .5 * n3 + .3 * n3n3 + .2 * n3n3n3;  // + O(n^4)
    phi3 = 1. - .125 * n3 - .05 * n3n3 - .025 * n3n3n3;  // + O(n^4)
  } else {
    phi2 = (6. * n3 - 3. * n3n3 + 6. * oneMn3 * logoneMn3) / n3n3n3;
    phi3 = (6. * n3 - 9. * n3n3 + 6. * n3n3n3 + 6. * oneMn3Squared * logoneMn3)/
        (4. * n3n3n3);
  }
  phi2num = 1. + n3n3 * phi2 / 9.;
  phi3num = 1. - 4. * n3 * phi3 / 9.;
  // Calculate free energy density terms
  Phi1 = -n0 * logoneMn3;
  Phi2 = phi2num *  (n1 * n2 - nvec1nvec2) / oneMn3;
  Phi3 = phi3num * (n2n2n2 - 3. * n2 * nvec2nvec2 + 4.5 * (
      nvec2ntensornvec2 - trace3)) / (24. * M_PI * oneMn3Squared);
  return Phi1 + Phi2 + Phi3;
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
