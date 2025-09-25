// Test the function charge_distribution for calculating the
// electrode charge distribution via application to simple problems with
// known solutions (planar electrodes, linear or constant potential...)
#include "../../../../parameter_handler/src/parameter_handler.hpp"
#include "../../../src/data_frame.hpp"
#include "../../../src/properties.hpp"
#include "../../../src/constants.hpp"
#include "../../../src/boundary_surface.hpp"
#include "../../../src/boundary_surface_planar.hpp"
#include "../../../src/surface_charge_distribution.hpp"
#include "../../../src/system.hpp"
#include <vector>
#include <string> // for parameter labels
#include <cmath> // std::exp
#include <math.h> // M_PI
#include <iostream>
#include <fstream> // write results to file

int main(int argc, char** argv) {
  // Get parameters from file/commandline input
  ParameterHandler params(argc, argv);
  params.process_parameters();
  System<3> system(&params);
  // System properties:
  std::vector<size_t> grid_counts = system.grid_counts;
  std::vector<double> bin_sizes = system.bin_sizes;
  std::vector<size_t> xy_counts{grid_counts.at(0), grid_counts.at(1)};
  double resolution = bin_sizes.at(0);
  std::vector<double> system_lengths = system.system_lengths;
  // Set boundary positions
  BoundarySurfacePlanar<3> surface(system, 0);
  std::cout << "Calculate electrode charge...\n";
  // potential on the electrodes
  double potential = params.get_double("potential");
  double temperature;
  system.get_property<double>("temperature", &temperature);
  std::cout << "potential = " << potential << " V "; // in Volt
  potential *= ELECTRON_CHARGE / (BOLTZMANN * temperature) * 1e4;
  std::cout << "= " << potential << " kT/e\n"; // in reduced units
  surface.set_boundary_value(potential);
  // define potential profile
  std::cout << "define potential profile\n";
  DataFrame<3,double> total_ES(grid_counts);
  double potential_step = -2. * potential / static_cast<double>(grid_counts.at(2) - 1); // i.e. value at position gc-1 is -potential
  std::cout << "potential step " << potential_step << std::endl;
  double value = potential - potential_step; // first potential step is added in for loop, yielding potential as value at z=0
  for (size_t k = 0; k < grid_counts.at(2); ++k) {
    value += potential_step;
    for (size_t i = 0; i < grid_counts.at(0); ++i) {
      for (size_t j = 0; j < grid_counts.at(1); ++j) {
	total_ES.at(i,j,k) = value;
      }
    }
  }
  std::vector<std::pair<std::vector<double>, double>> charge_normal(0);
  double total_charge_normal;
  std::vector<std::pair<std::vector<double>, double>> charge_poisson(0);
  double total_charge_poisson;
  std::cout << "calculate charge distr\n";
  charge_distribution(surface, system, total_ES,
		  &charge_normal, &total_charge_normal, resolution);
//  charge_distribution_poisson(surface, system, total_ES,
//		  &charge_poisson, &total_charge_poisson, resolution);
  // expected result (total)
  double diff_expected = 2. * potential * system_lengths.at(0) * system_lengths.at(1) /
	(system.bjerrum * 4. * M_PI * system_lengths.at(2));
  // Write results to file
  std::cout << "print results\n";
  std::string filename = "electrode_charge.dat";
  std::fstream file;
  file.open(filename, std::ios::out);
  file << "# [x] [y] [z] [surface charge density [e] (normal / Poisson)]\n";
  file << "# total charge: ";
//  for (double val : total_charge_diff) { file << std::to_string(val) + " "; }
  file << "normal: " << total_charge_normal << " ";
  file << "Poisson: " << total_charge_poisson << " ";
  file << "- expected: +-" << diff_expected;
  file << "\n";
  file << "# grid_counts=";
  for (size_t gc : grid_counts) {
    file << gc << ",";
  }
  file << std::endl;
  for (size_t i = 0; i  < charge_normal.size(); i++) {
    for (double pos : charge_normal.at(i).first) {
      file << pos << " ";
    }
    file << charge_normal.at(i).second << std::endl; //" ";
//    file << charge_poisson.at(i).second << std::endl;
  }
//  for (std::pair<std::vector<double>, double>& point : charge_diff) {
//    for (double pos : point.first) {
//      file << pos << " ";
//    }
//    file << point.second << std::endl;
//  }
  file.close();
  // If something is wrong: check potential gradients
//  double grad = -2. * potential / system_lengths.at(2);
//  std::cout << "expected potential gradient: " << grad << std::endl;
//  std::cout << "expected potential at bin size " << bin_sizes.at(2) << ": ";
//  std::cout << potential + grad * bin_sizes.at(2) << std::endl;
//  file.open("potential.dat", std::ios::out);
//  file << "# [x] [y] [z] [electrostatic potential]\n";
//  system.print_data(total_ES, file);
//  file.close();
  // all done
  std::cout << "see results in " << filename << std::endl;
  std::cout << "relative deviation from expected results: ";
  std::cout << (total_charge_normal - diff_expected) / diff_expected << std::endl;
}
