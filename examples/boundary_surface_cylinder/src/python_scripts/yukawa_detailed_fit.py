# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Fit a Yukawa model to the DFT results for the electrostatic potential.
# Here, we compose it of multiple independent point charges, instead of assuming
# a priori that the charge in the middle of the cylinder is constant.

# ______________________________________________________________________________
# Imports
# ______________________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import sys
import potential_models as model
sys.path.append("/home/fd160/Documents/capdft/examples/plots/src")
import load_data as DAT

# ______________________________________________________________________________
# Extract DFT result from data file
# Here we consider the total electrostatic potential produced by the colloid
# ('external potential') and the surrounding ion layer ('excess potential').
# ______________________________________________________________________________
print("Extract DFT result from data file")
electrolyte_default = ["full", "vacuum"][0]
if electrolyte_default == "vacuum":
    x,y,z,V_DFT = np.loadtxt("extpot_es.dat", skiprows=2).transpose()
else:
    x,y,z,V_DFT = np.loadtxt("total_ES.dat", skiprows=2).transpose()
# ______________________________________________________________________________
# Define system properties
# We need the system dimensions, height and radius of the cylinder, and the
# potential boundary value. These parameters can be extracted from the parameter
# file (params.txt) and the data file (total_ES.dat). Here we assume generally
# that the cylinder is positioned central within the system volume and aligned
# along the z-direction.
# ______________________________________________________________________________
print("Define system properties")
# system dimensions
grid_counts = np.array([DAT.get_param(f"{i}_grid_count") for i in range(3)])
system_lengths = np.array([DAT.get_param(f"{i}_length") for i in range(3)])
bin_sizes = system_lengths / grid_counts
# cylinder geometry
half_height = DAT.get_param("0_height") / 2
R = DAT.get_param("0_radius")
center = system_lengths / 2
center_indices = [round(center[i] / bin_sizes[i]) for i in range(3)]
center_index = (center_indices[2] + 
                center_indices[1] * grid_counts[2] + 
                center_indices[0] * grid_counts[2] * grid_counts[1])
# boundary value
boundary_potential = V_DFT[center_index] # boundary value, in kT/e
# ______________________________________________________________________________
# Estimate Yukawa parameters
# We apply a model where the cylinder is represented by a row of eleven charged
# spheres, each of which produces a Yukawa-type electrostatic potential. The
# screening length is assumed to be the same for all these Yukawa potentials.
# The interaction strength (amplitude of the Yukawa potential) may differ.
# However, we take advantage of the system's symmetry under inversion of z, so
# that we need only six, not eleven different parameters for the interaction
# strength.
# ______________________________________________________________________________
print("Estimate Yukawa parameters")
# model:
z_values = np.arange(0, half_height+R, 2*R) # z-positions of the spheres
def fit_function(positions, screening, A0, A1, A2, A3, A4, A5):
    V = model.yukawa1(positions, A0, screening, R) # centermost sphere
    for z,A in zip(z_values[1:], [A1,A2,A3,A4,A5]): # other spheres
        V += model.yukawa2(positions, A, screening, z, R)
    return V
# initial guess:
l = 1 # screening length, estimated [nm]
A = 0.4 * boundary_potential * np.exp(R/l) # potential amplitude, estimated
params = (1/l, *(len(z_values) * [A]))
print(f"# initial parameters: A={params[1]:.3}kT/e, "
      f"l={1/params[0]:.3}nm")
# ______________________________________________________________________________
# Choose suitable coordinate system and select fit region
# The implementation of the Yukawa potential used in the fit function above
# applies to a coordinate system centered around the cylinder center. I.e. to
# use it correctly, we need to shift the given positions appropriately.
# To save computational work, we can make use of the symmetries of the system:
# Due to (at least approximate) cylinder symmetry, we can perform the fit only
# within the x-z-plane instead of the whole 3d space. (Actually, we could limit
# the fit region even more strictly, e.g. to only positive x- & z-values, but
# most crucial is the dimensionality reduction.) Furthermore, we are not really
# concerned by what happens within the cylinder volume: Later when we add steric
# repulsion, close to the cylinder these steric interactions dominate over the
# electrostatic contributions. Therefore, we also exclude the cylinder-internal
# volume from the fit.
# ______________________________________________________________________________
print("Choose suitable coordinate system and select fit region")
# coordinate system:
positions = np.array([x,y,z]).transpose() - center
r = np.sqrt(positions[:,0]**2 + positions[:,1]**2)
z_centered = positions[:,2]
# fit region:
selection = np.where(y == y[center_index], True, False) # slice along middle x-z-plane
cylinder_volume = np.where(r <= R, True, False) * np.where(np.abs(z_centered) <= half_height, True, False)
top_volume = np.where(r**2 + (z_centered - half_height)**2 <= R**2, True, False)
bottom_volume = np.where(r**2 + (z_centered + half_height)**2 <= R**2, True, False)
cylinder_volume = (cylinder_volume | top_volume | bottom_volume)
selection *= (cylinder_volume == False) # exclude cylinder volume
# ______________________________________________________________________________
# Perform Yukawa fit
# Using the scipy.optimize.curve_fit algorithm, we fit the above described model
# to the electrostatic potential profile obtained from the DFT calculations.
# To quantify the quality of the fit, we calculate the mean squared and the
# maximum deviation between DFT and fit result. As explained previously, the
# potential profile within the cylinder volume is practically irrelevant for our
# purposes. Therefore, we do not include this region in our deviation calculus.
# ______________________________________________________________________________
print("Perform Yukawa fit")
params, covariance = scipy.optimize.curve_fit(fit_function,
                                              positions[selection],
                                              V_DFT[selection],
                                              p0 = params)
print(f"# optimal parameters: A = "
      "["+" ".join([f"{A:.3}" for A in params[1:]])+"]"
      f" nm kT/e, l={1/params[0]:.3}nm")
print("Apply fit to whole system volume")
# calculate potential profile in the Yukawa model:
V_MD = fit_function(positions, *params)
# set potential to constant boundary value within cylinder volume:
V_MD_adjusted = np.where(cylinder_volume, boundary_potential, V_MD)
# calculate mean squared / maximum deviation from DFT results (excluding the
# internal cylinder volume):
deviations = (V_MD_adjusted - V_DFT)[cylinder_volume == False]
deviation_rms = np.sqrt(sum(deviations**2)) / np.prod(deviations.shape)
deviation_max = max(abs(deviations))
print(f"Deviation from DFT result: root mean square {deviation_rms:.1}, "
      f"maximum {deviation_max:.1}")
# ______________________________________________________________________________
# Print results into file
# We print the optimized fit parameters and their (co-)variances. Also, we print
# the resulting potential profile, both in the full version and in the adjusted
# version where we imposed the fixed potential boundary condition.
# ______________________________________________________________________________
print("Print results into file")
positions = np.array([x,y,z]).transpose()
with open(electrolyte_default+"/yukawa_detailed_fit.dat", "w") as file:
    # header
    file.write(f"# Yukawa model\n")
    file.write(f"# parameters:\n")
    for i in range(len(z_values)):
        file.write(f"# A{i}: {params[i+1]:.3} kT/e\n")
    file.write(f"# l: {1/params[0]:.3} nm\n")
    try:
        file.write(f"# covariances: ")
        for line in covariance:
            file.write("(")
            for element in line:
                file.write(f"{element:.1},")
            file.write("),")
        file.write("\n")
        file.write(f"# deviation from DFT result:\n"
                   f"# root mean square deviation {deviation_rms:.1}\n"
                   f"# maximum deviation {deviation_max:.1}\n")
    except:
        pass
    file.write("# [x] [y] [z] [potential] [with constant boundary value]\n")
    file.write("# grid_counts=")
    for gc in grid_counts:
        file.write(f"{gc},")
    file.write("\n")
    # data
    for position, V, V_adj in zip(positions, V_MD, V_MD_adjusted):
        file.write(f"{position[0]} {position[1]} {position[2]} {V} {V_adj}\n")

