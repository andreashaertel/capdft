# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Plot Yukawa fit parameters A of the detailed fit (with six different charges)
# versus the more simple fit (with only two independent charge values).

# ______________________________________________________________________________
# Imports
# ______________________________________________________________________________
import matplotlib.pyplot as plt
import numpy as np
import sys
import load_data_local as dat
sys.path.append("/home/fd160/Documents/capdft/examples/plots/src")
import load_data as DAT

# ______________________________________________________________________________
# Initialization
# ______________________________________________________________________________
# The folders run<..> contain the data:
folders = sys.argv[1:]
# Choose a suitable unit for the charge plots:
name, plot_unit, unit_identifier = [
        ("amplitude", "nm kT/e", "A"),                              # option 0
        ("charge", "e", "Q"),                                       # option 1
        ("interaction strength", r"$\sqrt{\text{kT nm}}$", "Y")     # option 2
        ][1]                                                        # choice
# Choose with / without ion data:
electrolyte_default = ["full", "vacuum"][1]
title = "Charges in Yukawa model"
#if electrolyte_default == "vacuum":
#    title = "Yukawa parameters for cylinder in vacuum"
#else:
#    title = "Yukawa parameters for cylinder in electrolyte"
# Initialize data fields:
A_mid = [] # amplitude A of Yukawa potential on middle spheres (simple model)
A_tip = [] # amplitude A + A_offset of Yukawa potential on tips (simple model)
A_detailed = [] # list of amplitudes [A0,...,A5] in the detailed model 
screening_simple = [] # screening length of Yukawa potential (simple model)
screening_detailed = [] # screening length of Yukawa potential (detailed model)
bjerrum = [] # needed for the conversion from parameter A to sphere charge
sphere_radius = [] # needed for the charge-vs-position plot
# Initialize filters:
spherocylinders = [] # only spherocylinders
cylinders = [] # only cylinders with flat cap
valid = [] # omit folders where required parameters can't be found
# ______________________________________________________________________________
# Extract data
# ______________________________________________________________________________
for i,direc in enumerate(folders):
  subdirec = direc + "/" + electrolyte_default
  print(direc)
  # Extract parameters
  try:
    A, A_offset = dat.get_yukawa_amplitudes(subdirec, version="simple")
    l = dat.get_yukawa_screening(subdirec, version="simple")
    A_d = dat.get_yukawa_amplitudes(subdirec, version="detailed")
    l_d = dat.get_yukawa_screening(subdirec, version="detailed")
    valid.append(True)
    # Append values only after ensuring validity, to avoid that A is appended,
    # but l isn't:
    A_mid.append(A)
    A_tip.append(A + A_offset)
    screening_simple.append(l)
    A_detailed.append(A_d)
    screening_detailed.append(l_d)
    print(f"Screening lengths: {l:.3}nm - {l_d:.3}nm")
  except: 
    valid.append(False)
    print(f"Warning: Can't extract parameters from {subdirec}")
    continue
  # Determine sphero-/cylinder type
  cap_type = dat.get_cap_type(direc)
  if cap_type == "spherical":
    spherocylinders.append(True)
    cylinders.append(False)
  elif cap_type == "flat":
    cylinders.append(True)
    spherocylinders.append(False)
  # Determine Bjerrum length and sphere radius
  bjerrum.append(DAT.get_param("bjerrum_length", direc))
  sphere_radius.append(DAT.get_param("0_radius", direc))
A_mid = np.array(A_mid)
A_tip = np.array(A_tip)
screening_simple = np.array(screening_simple)
A_detailed = np.array(A_detailed)
screening_detailed = np.array(screening_detailed)
bjerrum = np.array(bjerrum)
sphere_radius = np.array(sphere_radius)
folders = np.array(folders)[valid]

# Conversion into requested unit:
if plot_unit == "nm kT/e":
    pass # parameters A already have the right unit
elif plot_unit == "e":
    # Convert from amplitude A [nm*kT/e] to charge Q [e]. This is defined
    # such that for r->0, the Yukawa potential becomes equivalent to the
    # Coulomb potential Q / (4*Pi*epsilon*r) = Q*kT*bjerrum_length / (e²*r):
    # Q[e] = A[nm*kT/e] / bjerrum_length[nm].
    conversion_factor = 1 / bjerrum
    A_tip *= conversion_factor
    A_mid *= conversion_factor
else:
    # The third option is chosen such that the pair interaction potential
    # between two spheres can be calculated in a straightforward way:
    # U_Yukawa(r) = Y1 * Y2 * e^(-r/l) / r, analogous to
    # U_Coulomb(r) = Q1 * Q2 / (4*Pi*epsilon*r).
    # Then the parameters Y have the unit sqrt(nm kT). Similar to the above
    # calculation, we can express the conversion in terms of the Bjerrum length:
    # Y[sqrt(nm kT)] = A[nm kT/e] / sqrt(bjerrum_length[nm]).
    conversion_factor = 1 / np.sqrt(bjerrum)
    A_tip *= conversion_factor
    A_mid *= conversion_factor
# Calculate total charge on cylinder:
A_total = 2 * A_tip + 9 * A_mid
# Pre-define plot labels:
def symbol(index):
    return f"${unit_identifier}_"+r"\text{"+index+"}$"
unit = f"[{plot_unit}]"

# ______________________________________________________________________________
# plot charge vs position:
# ______________________________________________________________________________
for i,direc in enumerate(folders[spherocylinders]):
    number_of_spheres = len(A_detailed[i])
    z_values = np.arange(number_of_spheres) * 2 * sphere_radius[i]
    plt.plot(z_values, A_detailed[i], f'C{i}+', label=direc)
for i,direc in enumerate(folders[spherocylinders]):
    number_of_spheres = len(A_detailed[i])
    z_values = np.arange(number_of_spheres) * 2 * sphere_radius[i]
    plt.plot(z_values, (number_of_spheres - 1) * [A_mid[i]] + [A_tip[i]], 
             f'C{i}-', label=direc)
plt.grid()
plt.legend(ncols=2, title="detailed model:   simple model:")
plt.xlabel(f"sphere position $z_i$ [nm]")
plt.ylabel(f"{name} {symbol('i')} per sphere {unit}")
plt.title(title)
plt.savefig(f"fit_comparison_{unit_identifier}.png",
            bbox_inches="tight")
plt.show()

# ______________________________________________________________________________
# plot charge ratio vs position:
# ______________________________________________________________________________
for i,direc in enumerate(folders[spherocylinders]):
    number_of_spheres = len(A_detailed[i])
    z_values = np.arange(number_of_spheres) * 2 * sphere_radius[i]
    A_values = A_detailed[i]
    plt.plot(z_values, A_values / A_values[0], f'C{i}+', label=direc)
for i,direc in enumerate(folders[spherocylinders]):
    number_of_spheres = len(A_detailed[i])
    z_values = np.arange(number_of_spheres) * 2 * sphere_radius[i]
    A_values = (number_of_spheres - 1) * [A_mid[i]] + [A_tip[i]]
    plt.plot(z_values, A_values / A_values[0], f'C{i}-', label=direc)
plt.grid()
plt.legend(ncols=2, title="detailed model:   simple model:")
plt.xlabel(f"sphere position $z_i$ [nm]")
plt.ylabel(f"{name} {symbol('i')}/{symbol('0')} per sphere\n"
           f"relative to {name} {symbol('0')} of centermost sphere")
plt.title(title)
plt.savefig(f"fit_comparison_ratio_{unit_identifier}.png",
            bbox_inches="tight")
plt.show()
