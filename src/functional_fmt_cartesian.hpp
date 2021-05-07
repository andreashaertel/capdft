//  Project CapDFT
//  functional_fmt_cartesian.cpp
//  WB-Mark II functional for hard spheres in cartesian coordinates
//  Copyright 2021 Philipp Pelagejcev

#ifndef SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
#define SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_

class FunctionalFMTCartesian {
 public:
  FunctionalFMTCartesian(void);
  ~FunctionalFMTCartesian(void);
  void calc_derivative(void);
  void calc_fluid_derivative(void);
  double **get_derivative_pointer(void);
  double **get_fluid_derivative_pointer(void);
};

#endif  // SRC_FUNCTIONAL_FMT_CARTESIAN_HPP_
