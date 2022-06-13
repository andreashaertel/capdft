// SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
/** \file examples/spherical_functional/src/main.cpp
 *  \brief Main file for DataFrame examples.
 *  
 *  This main file contains examples to show how the class DataFrame works.
 *
 */
// _____________________________________________________________________________
// Includes
#include <fftw3.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "../../../src/data_frame.hpp"
// _____________________________________________________________________________
// Main function
int main(int argc, char** args) {
// _____________________________________________________________________________
  /* We start by defining some geometric and physical properties.
   */
// _____________________________________________________________________________
  size_t grid_count = 21;
  size_t gc_x = 3;
  size_t gc_y = 3;
  size_t gc_z = 5;
// _____________________________________________________________________________
  /* Now, we create one DataFrame object with the 1D constructor, which only
   * requires the number of grid points.
   * A second DataFrame object is created with the copy constructor.
   * The first template argument determines the dimensionality (i.e. 1, 2, 3)
   * and the second one the data type (i.e. "double", "fftw_complex").
   */
// _____________________________________________________________________________
  DataFrame<1, double> my_data_frame1(grid_count);
  DataFrame<1, double> my_data_frame2(my_data_frame1);
  DataFrame<1, fftw_complex> my_data_frame3(grid_count);
  DataFrame<1, fftw_complex> my_data_frame4(my_data_frame3);
  DataFrame<3, double> my_data_frame5(std::vector<size_t>{gc_x, gc_y, gc_z});
// _____________________________________________________________________________
  /* The data in these objects is accessed with at() and element(), where at()
   * can read and write, but element() can only read.
   * The size() function just returns the grid point count, even if the
   * DataFrame has two or three dimensions.
   * The access functions check for out-of-bounds indices.
   * A quick way to set the same value to all DataFrame elements is
   * set_all_elements_to().
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
  for (size_t i = 0; i < gc_z; ++i) {
    for (size_t j = 0; j < gc_y; ++j) {
      for (size_t k = 0; k < gc_x; ++k) {
        my_data_frame5.at(k, j, i) = static_cast<double>((i+1)*(j+1)*(k+1));
      }
    }
  }
  std::cout << "3D DataFrame size: " << my_data_frame5.size() << std::endl;
  // my_data_frame1.at(2001);  // error: out of bounds (2001>1001)
// _____________________________________________________________________________
  /* To show the values in these containers, a print() function is available.
   * The DataFrame class provides the operators =,+,-,*,/,+=,-=,*=,/= and the math
   * functions exp(), log_natural(), abs() and max().
   */
// _____________________________________________________________________________
  DataFrame<1, double> result(my_data_frame1);
  my_data_frame1 *= my_data_frame2;
  my_data_frame1 /= my_data_frame2;
  my_data_frame1 += my_data_frame2;
  my_data_frame1 -= my_data_frame2;
  my_data_frame1.print();
  result = 4.2 * exp(my_data_frame1 + my_data_frame2);
  result.print();
  (my_data_frame3-my_data_frame4).print();
  my_data_frame5.print();
// _____________________________________________________________________________
  return 0;
}
