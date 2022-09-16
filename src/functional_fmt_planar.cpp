// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>  // NOLINT
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_planar.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "data_frame.hpp"  // NOLINT
// _____________________________________________________________________________
FunctionalFMTPlanar::FunctionalFMTPlanar() {
  //
}
// _____________________________________________________________________________
FunctionalFMTPlanar::FunctionalFMTPlanar(
    const std::vector<DataFrame<1, double>>* density_profiles,
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
FunctionalFMTPlanar::FunctionalFMTPlanar(
      const std::vector<DataFrame<1, double>>* density_profiles,
      const std::vector<Properties>& species_properties,
      const Properties& system_properties)
  : FunctionalFMTPlanar(
      density_profiles, species_properties, system_properties,
      std::vector<size_t>(0)) {
}
// _____________________________________________________________________________
FunctionalFMTPlanar::~FunctionalFMTPlanar() {
  //
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::extract_system_properties(
    const Properties& system_properties) {
  system_properties.get_property("length", &length);
  system_properties.get_property("grid count", &grid_count);
  // Calculate bin sizes
  dz = length / static_cast<double>(grid_count);
  dkz = 2. * M_PI / length;
  // Since we only do real-->complex Fourier transforms (or vice versa),
  // we only need half of the Fourier space (see FFTW docs).
  grid_count_fourier = grid_count / 2 + 1;
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::extract_species_properties(
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
          std::cerr << "FunctionalFMTPlanar::extract_species_properties(): ";
          std::cerr << "\"Error: A species with a diameter but no density ";
          std::cerr << "was detected.\"" << std::endl;
          exit(1);
        }
        affected_species.push_back(it - spec_prop.begin());
      }
    }
  }
  // Extract properties
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (!spec_prop.at(*it).get_property("diameter", &diameter) ||
        !spec_prop.at(*it).get_property("bulk density", &bulk_density)) {
      std::cerr << "FunctionalFMTPlanar::extract_species_properties(): ";
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
void FunctionalFMTPlanar::initialize_all_data_frames() {
  // Initialize internal density profiles
  for (size_t i = 0; i < species_count; ++i) {
    density_profiles.push_back(DataFrame<1, double>(grid_count));
    density_profiles_four.push_back(
        DataFrame<1, fftw_complex>(grid_count_fourier));
  }
  // Initialize weights
  for (size_t i = 0; i < species_count; ++i) {
    weights_four.push_back(std::vector<DataFrame<1, fftw_complex>>(0));
    for (size_t j = 0; j < 8; ++j) {
      weights_four.at(i).push_back(
          DataFrame<1, fftw_complex>(grid_count_fourier));
    }
  }
  // Initialize weighted densities
  for (size_t i = 0; i < 4; ++i) {
    scalar_weighted_dens_real.push_back(DataFrame<1, double>(grid_count));
    scalar_weighted_dens_four.push_back(
        DataFrame<1, fftw_complex>(grid_count_fourier));
  }
  for (size_t i = 0; i < 2; ++i) {
    vector_weighted_dens_real.push_back(DataFrame<1, double>(grid_count));
    tensor_weighted_dens_real.push_back(DataFrame<1, double>(grid_count));
    vector_weighted_dens_four.push_back(
        DataFrame<1, fftw_complex>(grid_count_fourier));
    tensor_weighted_dens_four.push_back(
        DataFrame<1, fftw_complex>(grid_count_fourier));
  }
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::update_density_profiles() {
  size_t index{0};
  for (size_t i = 0; i < species_count; ++i) {
    index = affected_species.at(i);
    density_profiles.at(i) = density_profiles_pointer->at(index);
  }
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_weights() {
  // Declare auxiliary variables. R means radius of the hard sphere.
  // kz is the z-component of the position in Fourier space.
  // RR is the square of R, kzkz is the square of kz, and so on.
  double kz{0.}, kzkz{0.}, kzkzkz{0.};
  double R{0.}, RR{0.}, RRR{0.};
  double Rkz{0.};
  set_weights_to_zero();
  for (size_t i = 0; i != species_count; ++i) {
    R = diameters.at(i) / 2.;
    RR = R * R;
    RRR = RR * R;
    // Calculate limits kz --> 0
    // scalar weight 3 (theta)
    weights_four.at(i).at(0).at(0)[0] = 4. * M_PI * RRR / 3.;  // real
    // scalar weight 2 (delta)
    weights_four.at(i).at(1).at(0)[0] = 4. * M_PI * RR;
    // scalar weight 1 (delta)
    weights_four.at(i).at(2).at(0)[0] = R;
    // scalar weight 0 (delta)
    weights_four.at(i).at(3).at(0)[0] = 1;
    // Only calculate the weights for the positive frequencies in Fourier space
    for (size_t j = 1; j != grid_count_fourier; ++j) {
      // Absolute value of Fourier space coordinate and kz are equal, since we
      // can choose kx=ky=0 due to the symmetry and we only have to consider
      // positive kz.
      kz = dkz * static_cast<double>(j);
      Rkz = R * kz;
      kzkz = kz * kz;
      kzkzkz = kzkz * kz;
      // scalar weight 3 (theta)
      weights_four.at(i).at(0).at(j)[0] = 4. * M_PI *
          (sin(Rkz) - Rkz * cos(Rkz)) / kzkzkz;  // real
      weights_four.at(i).at(0).at(j)[1] = 0.;  // imaginary
      // scalar weight 2 (delta)
      weights_four.at(i).at(1).at(j)[0] = 4. * M_PI * R * sin(Rkz) / kz;
      weights_four.at(i).at(1).at(j)[1] = 0.;
      // scalar weight 1 (delta)
      weights_four.at(i).at(2).at(j)[0] = sin(Rkz) / kz;
      weights_four.at(i).at(2).at(j)[1] = 0.;
      // scalar weight 0 (delta)
      weights_four.at(i).at(3).at(j)[0] = sin(Rkz) / Rkz;
      weights_four.at(i).at(3).at(j)[1] = 0.;
      // vectorial weight 2 3rd element (1st and 2nd element are zero)
      weights_four.at(i).at(4).at(j)[0] = 0.;
      weights_four.at(i).at(4).at(j)[1] = 4. * M_PI *
          (Rkz * cos(Rkz) - sin(Rkz)) / kzkz;
      // vectorial weight 1 3rd element (1st and 2nd element are zero)
      weights_four.at(i).at(5).at(j)[0] = 0.;
      weights_four.at(i).at(5).at(j)[1] = (kz * cos(Rkz) - sin(Rkz) / R) / kzkz;
      // tensorial weight (1,1)-element identical to (2,2)-element
      weights_four.at(i).at(6).at(j)[0] =
          -4. * M_PI * R * sin(Rkz) / (3 * kz) +
          4. * M_PI * (sin(Rkz) - Rkz * cos(Rkz)) / (R * kzkzkz);
      weights_four.at(i).at(6).at(j)[1] = 0.;
      // tensorial weight (3,3)-element
      weights_four.at(i).at(7).at(j)[0] =
          8. * M_PI * R * sin(Rkz) / (3 * kz) +
          -8. * M_PI * (sin(Rkz) - Rkz * cos(Rkz)) / (R * kzkzkz);
      weights_four.at(i).at(7).at(j)[1] = 0.;
    }
  }
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::set_weights_to_zero() {
  fftw_complex zero{0., 0.};
  for (size_t i = 0; i != species_count; ++i) {
    for (size_t j = 0; j != weights_four.at(i).size(); ++j) {  // different weights
      weights_four.at(i).at(j).set_all_elements_to(zero);
    }
  }
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_derivative(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != grid_count) {
      std::cerr << "FunctionalFMTPlanar::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  // Calculate the weighted densities
  calc_weighted_densities();
  //// From the weighted densities calculate the partial derivatives of the
  //// excess free energy
  //calc_partial_derivatives();
  //// Calculate the weighted partial derivatives
  //calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_derivative_warnings(
    std::vector<DataFrame<1, double>>* functional_derivative) {
  // Check if the given DataFrame has correct array length
  for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
    if (functional_derivative->at(*it).size() != grid_count) {
      std::cerr << "FunctionalFMTPlanar::calc_derivative(): \"";
      std::cerr << "Error: Supplied DataFrame has incorrect array size.\"";
      std::cerr << std::endl;
      exit(1);
    }
  }
  //// Calculate the weighted densities
  //calc_weighted_densities();
  //// Check whether there are unphysical values in the  weighted densities
  //check_weighted_densities();
  //// From the weighted densities calculate the partial derivatives of the
  //// excess free energy
  //calc_partial_derivatives();
  //// Calculate the weighted partial derivatives
  //calc_weighted_partial_derivatives(functional_derivative);
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_bulk_derivative(
    std::vector<double>* bulk_derivative) {
  //// Auxiliary variables
  //double diameter{0.}, diameter2{0.}, diameter3{0.};
  //// Calculate the fluid weighted densities and save them to the first position
  //for (auto it = scalar_weighted_dens_real.begin();
  //    it != scalar_weighted_dens_real.end(); ++it) {
  //  it->at(0) = 0.;
  //}
  //for (size_t i = 0; i != species_count; ++i) {
  //  // Define auxiliary values
  //  diameter = diameters.at(i);
  //  diameter2 = diameter * diameter;
  //  diameter3 = diameter2 * diameter;
  //  // All non-scalar weighted densities are 0 the others are obviously constant
  //  scalar_weighted_dens_real.at(0).at(0) +=
  //      bulk_densities.at(i) * diameter3 * M_PI / 6.0;
  //  scalar_weighted_dens_real.at(1).at(0) +=
  //      bulk_densities.at(i) * diameter2 * M_PI;
  //  scalar_weighted_dens_real.at(2).at(0) +=
  //      bulk_densities.at(i) * diameter / 2.0;
  //  scalar_weighted_dens_real.at(3).at(0) += bulk_densities.at(i);
  //}
  //// Evaluate the partial derivatives at the first position
  //calc_local_partial_derivatives(0);
  //std::fill(bulk_derivative->begin(), bulk_derivative->end(), 0.);
  //size_t spec_i{0};
  //for (auto it = affected_species.begin(); it != affected_species.end(); ++it) {
  //  spec_i = it - affected_species.begin();
  //  // Define auxiliary values
  //  diameter = diameters.at(spec_i);
  //  diameter2 = diameter * diameter;
  //  diameter3 = diameter2 * diameter;
  //  // Calculate the derivatives of the excess free energy functional w.r.t. the
  //  // weighted densities.
  //  bulk_derivative->at(*it) = scalar_partial_derivative_real.at(0).at(0) *
  //      diameter3 * M_PI / 6.0;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(1).at(0) *
  //      diameter2 * M_PI;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(2).at(0) *
  //      diameter / 2.0;
  //  bulk_derivative->at(*it) += scalar_partial_derivative_real.at(3).at(0);
  //}
}
// _____________________________________________________________________________
double FunctionalFMTPlanar::calc_energy() {
  double integral{0.};
  //DataFrame<1, double> free_energy_density(grid_count);
  //// Calculate the weighted densities for the current density profile
  //calc_weighted_densities();
  //check_weighted_densities();
  //// Store free energy density
  //for (size_t i = 0; i < grid_count; ++i) {
  //  free_energy_density.at(i) = calc_local_energy_density(i);
  //}
  //integral = radial_integration(free_energy_density.array(), grid_count, dr);
  return integral;
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::calc_weighted_densities() {
  std::vector<fftw_plan> forward_plans;
  std::vector<fftw_plan> backward_plans;
  // Specify the plans
  // (doing this locally produces just little numerical overhead)
  for (size_t i = 0; i != species_count; ++i) {
    forward_plans.push_back(
        fftw_plan_dft_r2c_1d(grid_count,
            density_profiles.at(i).array(),  // input
            density_profiles_four.at(i).array(),  // output
            flags_keep));
  }
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_1d(grid_count,
            scalar_weighted_dens_four.at(i).array(),
            scalar_weighted_dens_real.at(i).array(),
            flags_destroy));
  }
  for (size_t i = 0; i != vector_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_1d(grid_count,
            vector_weighted_dens_four.at(i).array(),
            vector_weighted_dens_real.at(i).array(),
            flags_destroy));
  }
  for (size_t i = 0; i != tensor_weighted_dens_four.size(); ++i) {
    backward_plans.push_back(
        fftw_plan_dft_c2r_1d(grid_count,
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
    profile_four *= dz;
  }
  // Convolution
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    for (size_t j = 0; j != species_count; ++j) {
      scalar_weighted_dens_four.at(i) +=
          density_profiles_four.at(j) * weights_four.at(j).at(i);
    }
  }
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    for (size_t j = 0; j != species_count; ++j) {
      vector_weighted_dens_four.at(i) +=
          density_profiles_four.at(j) * weights_four.at(j).at(i);
    }
  }
  for (size_t i = 0; i != scalar_weighted_dens_four.size(); ++i) {
    for (size_t j = 0; j != species_count; ++j) {
      tensor_weighted_dens_four.at(i) +=
          density_profiles_four.at(j) * weights_four.at(j).at(i);
    }
  }
  // Forward Fourier transform of density profiles
  for (auto& plan : backward_plans) {
    fftw_execute(plan);
  }
  // Normalize
  for(auto& profile_four : scalar_weighted_dens_four) {  // every weighted dens
    profile_four *= dkz / (2. * M_PI);
  }
  for(auto& profile_four : vector_weighted_dens_four) {  // every weighted dens
    profile_four *= dkz / (2. * M_PI);
  }
  for(auto& profile_four : tensor_weighted_dens_four) {  // every weighted dens
    profile_four *= dkz / (2. * M_PI);
  }
  // The convolution may have numerical at both ends of the array.
  // One way of avoiding this is by using an outer external potential
  // forcing the density profile to vanish (i.e. two planar walls).
  for (auto& plan : forward_plans) { fftw_destroy_plan(plan); }
  for (auto& plan : backward_plans) { fftw_destroy_plan(plan); }
}
// _____________________________________________________________________________
void FunctionalFMTPlanar::check_weighted_densities() {
  //// initialize counters
  //size_t n3bigger1_counter = 0;
  //size_t n3smaller0_counter = 0;
  //size_t n2smaller0_counter = 0;
  //// Since calc_weighted_densities() is executed before this function,
  //// we can check if the weighted densities are in the legal range.
  //for (size_t i = 0; i < grid_count; ++i) {
  //  // Check scalar weighted densities for values smaller 0
  //  if (scalar_weighted_dens_real.at(0).at(i) < 0.) {
  //    ++n3smaller0_counter;
  //    scalar_weighted_dens_real.at(0).at(i) = 0.;
  //  }
  //  if (scalar_weighted_dens_real.at(1).at(i) < 0.) {
  //    ++n2smaller0_counter;
  //    scalar_weighted_dens_real.at(1).at(i) = 0.;
  //    scalar_weighted_dens_real.at(2).at(i) = 0.;
  //    scalar_weighted_dens_real.at(3).at(i) = 0.;
  //  }
  //  // Check if n3 is bigger than 1
  //  if (scalar_weighted_dens_real.at(0).at(i) > 1.) {
  //    ++n3bigger1_counter;
  //    scalar_weighted_dens_real.at(0).at(i) = 1.;
  //  }
  //}
  //// Show warning
  //if (n3smaller0_counter > 0) {
  //  std::cerr << "FunctionalFMTPlanar::check_weighted_densities(): \"";
  //  std::cerr << "Warning: Local n3 < 0.0 at " << n3smaller0_counter;
  //  std::cerr << " positions.\"" << std::endl;
  //}
  //if (n2smaller0_counter > 0) {
  //  std::cerr << "FunctionalFMTPlanar::check_weighted_densities(): \"";
  //  std::cerr << "Warning: Local n2 < 0.0 at " << n2smaller0_counter;
  //  std::cerr << " positions.\"" << std::endl;
  //}
  //if (n3bigger1_counter > 0) {
  //  std::cerr << "FunctionalFMTPlanar::check_weighted_densities(): \"";
  //  std::cerr << "Warning: Local n3 > 1.0 at " << n3bigger1_counter;
  //  std::cerr << " positions.\"" << std::endl;
  //}
}
// _____________________________________________________________________________
// _____________________________________________________________________________
// _____________________________________________________________________________
