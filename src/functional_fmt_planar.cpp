// SPDX-FileCopyrightText: 2019 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-FileCopyrightText: 2019 Andreas Härtel <http://andreashaertel.anno1982.de/>
// SPDX-License-Identifier: LGPL-3.0-or-later
#include "functional_fmt_planar.hpp"  // NOLINT
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "data_frame.hpp"  // NOLINT
//// _____________________________________________________________________________
//FunctionalFMTPlanar::FunctionalFMTPlanar() {
//  //
//}
//// _____________________________________________________________________________
//FunctionalFMTPlanar::FunctionalFMTPlanar(
//    const std::vector<DataFrame<1, double>>* density_profiles,
//    const std::vector<Properties>& species_properties,
//    const Properties& system_properties,
//    const std::vector<size_t>& affected_species)
//  : affected_species(affected_species),
//    density_profiles_pointer(density_profiles) {
//  // Clear all std::vectors
//  diameters.clear();
//  bulk_densities.clear();
//  // Get system properties
//  extract_system_properties(system_properties);
//  // Get species properties; excludes all species without diameter property
//  extract_species_properties(species_properties);
//  // Initialize all data frames
//  initialize_all_data_frames();
//  // Calculate weights
//  calc_weights();
//}
//// _____________________________________________________________________________
//FunctionalFMTPlanar::FunctionalFMTPlanar(
//      const std::vector<DataFrame<1, double>>* density_profiles,
//      const std::vector<Properties>& species_properties,
//      const Properties& system_properties)
//  : FunctionalFMTPlanar(
//      density_profiles, species_properties, system_properties,
//      std::vector<size_t>(0)) {
//}
//// _____________________________________________________________________________
//FunctionalFMTPlanar::~FunctionalFMTPlanar() {
//  //
//}
