// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
 *  \brief Main file for DataFrame examples.
 *  
 *  This main file contains examples to show how the class DFSpherical works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <cmath>
#include <fftw3.h>
#include "../../../src/properties.hpp"
#include "../../../src/df_spherical.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  /* We start by defining some geometric and physical properties. Some of them
   * will be used to create any DF* object, but others will be neglected.
   * To see how the Properties class works see the "properties" example.
   */
// _____________________________________________________________________________
  size_t grid_count = 101;
  double system_length = 19.821782178217823;
  // Create object of Properties class
  Properties system_properties;
  // Put the system properties into a Properties container
  system_properties.add_property<double>("length", system_length);
  system_properties.add_property<size_t>("grid count", grid_count);
// _____________________________________________________________________________
  /* Now, we create one DFSpherical object from the system_properties.
   * The DFSpherical class is used for handling data in the spherical
   * geometry. Its constructor only uses the "grid count" property.
   * A second DFSpherical object is created with the copy constructor.
   */
// _____________________________________________________________________________
  DFSpherical<double> my_data_frame1(system_properties);
  DFSpherical<double> my_data_frame2(my_data_frame1);
  DFSpherical<fftw_complex> my_data_frame3(system_properties);
  DFSpherical<fftw_complex> my_data_frame4(my_data_frame3);
// _____________________________________________________________________________
  /* In the spherical geometry the data in these objects is accessed with at()
   * and element(), where at() can read and write, but element() can only read.
   * Because the objects are 1D arrays, the size() function just returns the
   * grid_count value. The access functions check for out-of-bounds indices.
   */
// _____________________________________________________________________________
  for (size_t i = 0; i < my_data_frame1.size(); ++i) {
    my_data_frame1.at(i) = static_cast<double>(i+1);
    my_data_frame2.at(i) = -static_cast<double>(i+1);
  }
  for (size_t i = 0; i < my_data_frame1.size(); ++i) {
    my_data_frame3.at(i)[0] = static_cast<double>(i);
    my_data_frame3.at(i)[1] = static_cast<double>(i);
    my_data_frame4.at(i)[0] = -static_cast<double>(i);
    my_data_frame4.at(i)[1] = -static_cast<double>(i);
  }
  // my_data_frame1.at(2001);  // error: out of bounds (2001>1001)
// _____________________________________________________________________________
  /* To show the values in these containers, a print() function is available.
   * Every DF* class provides the operators =,+,-,*,/,+=,-=,*=,/= and the math
   * functions exp() and log_natural().
   */
// _____________________________________________________________________________
  DFSpherical<double> result(my_data_frame1);
  my_data_frame1 *= my_data_frame2;
  my_data_frame1 /= my_data_frame2;
  my_data_frame1 += my_data_frame2;
  my_data_frame1 -= my_data_frame2;
  my_data_frame1.print();
  result = 4.2 * exp(my_data_frame1 + my_data_frame2);
  result.print();
  (my_data_frame3-my_data_frame4).print();
// _____________________________________________________________________________
  return 0;
}
