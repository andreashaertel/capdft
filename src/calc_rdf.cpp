// Copyright 2021 Andreas Haertel
#include <vector>
#include <regex>  // NOLINT
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include "parameter_handler.hpp"
int main(int argc, char** args) {
  // Activate the command line tool
  ParameterHandler cmdtool(argc, args);
  cmdtool.add_usage(" This program calculates the pair-distribution function ");
  cmdtool.add_usage(" from given particle positions. ");
  // register parameters of interest
  cmdtool.process_flag_help();
  cmdtool.process_parameters();
  // double boxx = cmdtool.get_double("box_len_x");
  // std::cout << "Boxx: " << boxx << std::endl;
  // std::vector<std::string> ert = cmdtool.get_remaining_cmdline_arguments();
  // for (auto it = ert.begin(); it != ert.end(); ++it) {
  //  std::cout << *it << std::endl;
  // }
  //
  // int flag_pos = 0;
  // while (cmdtool.find_flag("-t", flag_pos)) {
  //  std::cout << cmdtool.get_option(flag_pos, 1) << std::endl;
  //  std::cout << cmdtool.get_option(flag_pos, 2) << std::endl;
  // }
  // if (cmdtool.find_flag("-t")) {
  //  std::cout << cmdtool.find_option("-t", 1) << std::endl;
  // }
  //
  // return 0;
  // read in parameters
  double box_len_x = 0.0;
  double box_len_y = 0.0;
  double box_len_z = 0.0;
  int g_bins = 0;
  cmdtool.add_usage_description(" Parameter box_len_x: ... ");
  try {
    box_len_x = cmdtool.get_double("box_len_x");
    box_len_y = cmdtool.get_double("box_len_y");
    box_len_z = cmdtool.get_double("box_len_z");
    g_bins = cmdtool.get_int("rdf_bins");
  } catch (std::exception *e) {
    cmdtool.show_usage(std::cout);
    return 0;
  }
  double min_box_len = box_len_x;
  if (box_len_y < min_box_len) { min_box_len = box_len_y; }
  if (box_len_z < min_box_len) { min_box_len = box_len_z; }
  double g_len = cmdtool.get_double("rdf_length", 0.5*min_box_len);
  // check for errors
  // cmdtool.check_errors(std::cout);
  if (argc < 3) {
    std::cout << " please specify an input file ... " << std::endl;
    return 0;
  }
  // get first non-parameter argument:
  // CHANGE TO -P ...
  // Dummy string
  std::string line;
  // Regular expression to filter comments beforehand
  std::regex comment_regex("^([^#]*)");
  // Regular expression to get name-value pairs:
  std::regex param_regex("^\\s*([\\w\\.]+)\\s+([\\w\\.]+)\\s+([\\w\\.]+)\\s+([\\w\\.]+)\\s+([\\w\\.]+)");  // NOLINT
  // Make smatch object containing the filtered words
  std::smatch param_match;
  std::smatch comment_match;
  // Try to open the parameter file
  std::fstream parameterFile(args[1], std::ios::in);
  if (!parameterFile.is_open()) {
    // if the file couldnt be opened, throw error
    std::cout << " problem with file ... " << std::endl;
    return 0;
  }
  int N1 = 0;
  int N2 = 0;
  std::vector<double*> part;
  double *entry;
  while (getline(parameterFile, line)) {  // read file line by line
    std::regex_search(line, comment_match, comment_regex);  // filter comment
    line = comment_match.str(1);
    if (std::regex_search(line, param_match, param_regex)) {  // filter param
      entry = new double[5];
      entry[0] = std::stod(param_match.str(1), NULL);
      entry[1] = std::stod(param_match.str(2), NULL);
      entry[2] = std::stod(param_match.str(3), NULL);
      entry[3] = std::stod(param_match.str(4), NULL);
      entry[4] = std::stod(param_match.str(5), NULL);
      part.push_back(entry);
      if (entry[4] == 3) { N1++; }
      if (entry[4] == 4) { N2++; }
    }
  }
  double g_d = g_len / static_cast<double>(g_bins);
  double *gr = new double[g_bins];
  double *gvol = new double[g_bins];
  uint64_t *g11long = new uint64_t[g_bins];
  uint64_t *g12long = new uint64_t[g_bins];
  uint64_t *g22long = new uint64_t[g_bins];
  double *g11 = new double[g_bins];
  double *g12 = new double[g_bins];
  double *g22 = new double[g_bins];
  double r1, r2;
  for (int i = 0; i < g_bins; i++) {
    gr[i] = g_d * (static_cast<double>(i) + 0.5);
    r2 = g_d * static_cast<double>(i+1);
    r1 = g_d * static_cast<double>(i);
    gvol[i] = 4.0 / 3.0 * M_PI * (r2*r2*r2-r1*r1*r1);
    g11long[i] = 0;
    g12long[i] = 0;
    g22long[i] = 0;
    g11[i] = 0.0;
    g12[i] = 0.0;
    g22[i] = 0.0;
  }
  uint64_t g11_N = 0;
  uint64_t g12_N = 0;
  uint64_t g22_N = 0;
  double dist, distx, disty, distz;
  int pos;
  int t1, t2;
  for (auto it1 = part.begin(); it1 != part.end(); ++it1) {
    for (auto it2 = part.begin(); it2 != part.end(); ++it2) {
      if (it1 == it2) { continue; }
      t1 = static_cast<int>((*it1)[4]);
      t2 = static_cast<int>((*it2)[4]);
      if ((t1 < 3) || (t2 < 3)) { continue; }
      distx = (*it1)[1] - (*it2)[1];
      while (distx < -0.5*box_len_x) { distx = distx + box_len_x; }
      while (distx >= 0.5*box_len_x) { distx = distx - box_len_x; }
      disty = (*it1)[2] - (*it2)[2];
      while (disty < -0.5*box_len_y) { disty = disty + box_len_y; }
      while (disty >= 0.5*box_len_y) { disty = disty - box_len_y; }
      distz = (*it1)[3] - (*it2)[3];
      while (distz < -0.5*box_len_z) { distz = distz + box_len_z; }
      while (distz >= 0.5*box_len_z) { distz = distz - box_len_z; }
      dist = sqrt(distx*distx + disty*disty + distz*distz);
      if (dist >= g_len) { continue; }
      pos = static_cast<int>(floor(dist / g_d));
      if (t1 == 3) {
        if (t2 == 3) {
          g11long[pos] = g11long[pos] + 1;
          g11_N++;
        }
        if (t2 == 4) {
          g12long[pos] = g12long[pos] + 1;
          g12_N++;
        }
      }
      if (t1 == 4) {
        // if (t2 == 3) {
        //  g12[pos] = g12[pos] + 1;
        //  g12_N++;
        // }
        if (t2 == 4) {
          g22long[pos] = g22long[pos] + 1;
          g22_N++;
        }
      }
    }
  }
  // normalize
  double vol_tot = box_len_x * box_len_y * box_len_z;
  uint64_t longN1 = static_cast<uint64_t>(N1);
  uint64_t longN2 = static_cast<uint64_t>(N2);
  for (int i = 0; i < g_bins; i++) {
    // normalization:
    // Tracked particle: partX1
    // no(partX1) * densX2 * vol = no(partX1) * no(partX2) * vol / vol_tot
    g11[i] = static_cast<double>(g11long[i]) * vol_tot /
        (static_cast<double>(longN1*longN1) * gvol[i]);
    g12[i] = static_cast<double>(g12long[i]) * vol_tot /
        (static_cast<double>(longN1*longN2) * gvol[i]);
    g22[i] = static_cast<double>(g22long[i]) * vol_tot /
        (static_cast<double>(longN2*longN2) * gvol[i]);
  }
  // Output data
  // First, set precision to maximum for double values
  typedef std::numeric_limits< double > dbl;
  std::cout.precision(dbl::max_digits10);
  // Write header
  std::cout << "# g(r) calculated from: " << args[1] << std::endl;
  std::cout << "# Number of particles species 1 (previously 3): " << N1
      << std::endl;
  std::cout << "# Number of particles species 2 (previously 4): " << N2
      << std::endl;
  std::cout << "# Total volume of box: " << vol_tot << " = " << box_len_x
      << " x " << box_len_y << " x " << box_len_z << std::endl;
  std::cout << "# Bins for g(r): " << g_bins << std::endl;
  std::cout << "# Length g is sampled on: " << g_len << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Column 1: r" << std::endl;
  std::cout << "# Column 2: g11 (previously 3,3)" << std::endl;
  std::cout << "# Column 3: g12 (previously 3,4)" << std::endl;
  std::cout << "# Column 4: g22 (previously 4,4)" << std::endl;
  // Write data
  for (int i = 0; i < g_bins; i++) {
    std::cout << gr[i] << " ";
    std::cout << g11[i] << " ";
    std::cout << g12[i] << " ";
    std::cout << g22[i] << std::endl;
  }
  std::cout << std::endl;
  // destroy
  for (auto it = part.begin(); it != part.end(); ++it) {
    delete[] (*it);
  }
  delete[] gr;
  delete[] gvol;
  delete[] g11;
  delete[] g12;
  delete[] g22;
  delete[] g11long;
  delete[] g12long;
  delete[] g22long;
  return 0;
}
