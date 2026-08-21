# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Plot Yukawa fit parameters for different runs.

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
        ][2]                                                        # choice
# Choose with / without ion data:
electrolyte_default = "vacuum"
if electrolyte_default == "vacuum":
    title = "Yukawa parameters for cylinder in vacuum"
else:
    title = "Yukawa parameters for cylinder in electrolyte"
# Initialize data fields:
A_mid = [] # amplitude A of Yukawa potential on middle spheres
A_tip = [] # amplitude A + A_offset of Yukawa potential on tips
screening = [] # screening length of Yukawa potential
bjerrum = [] # needed for the conversion from parameter A to sphere charge
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
    A, A_offset = dat.get_yukawa_amplitudes(subdirec)
    l = dat.get_yukawa_screening(subdirec)
    valid.append(True)
    # Append values only after ensuring validity, to avoid that A is appended,
    # but l isn't:
    A_mid.append(A)
    A_tip.append(A + A_offset)
    screening.append(l)
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
  # Determine Bjerrum length
  bjerrum.append(DAT.get_param("bjerrum_length", direc))
A_mid = np.array(A_mid)
A_tip = np.array(A_tip)
screening = np.array(screening)
bjerrum = np.array(bjerrum)
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
# plot charge on tip vs mid:
# ______________________________________________________________________________
x = A_mid
y = A_tip
dx = max(x) - min(x)
dy = max(y) - min(y)
plt.scatter(x[spherocylinders], y[spherocylinders], c="C1", label="spherocylinder")
plt.scatter(x[cylinders], y[cylinders], c="C0", label="cylinder")
for direc,xval,yval in zip(folders,x,y):
  i = int(direc[-2:])
  shift1 = i % 2 * 0.02 * dx 
  shift2 = (i % 4 - i % 2) * 0.02 * dy
  plt.text(xval + shift1, yval - shift2, i, ha='right')

# plot diagonal axis for visual orientation
xrange = (min(x) - 0.1*dx, max(x) + 0.1*dx)
c = np.linspace(*xrange)
plt.plot(c,c, 'k',lw=0.3)
plt.xlim(*xrange)
plt.grid()
plt.legend()
plt.ylabel(f"{name} {symbol("tip")} per cylinder tip {unit}")
plt.xlabel(f"{name} {symbol("mid")} on the remaining spheres {unit}")
plt.title(title)
plt.savefig(f"fit_parameters_charges_{unit_identifier}.png",
            bbox_inches="tight")
plt.close()

# ______________________________________________________________________________
# plot charge ratio on tip / mid vs total charge:
# ______________________________________________________________________________
x = A_total
y = A_tip / A_mid
dx = max(x) - min(x)
dy = max(y) - min(y)
plt.scatter(x[spherocylinders], y[spherocylinders], c="C1", label="spherocylinder")
plt.scatter(x[cylinders], y[cylinders], c="C0", label="cylinder")
for direc,xval,yval in zip(folders,x,y):
  i = int(direc[-2:])
  shift1 = i % 2 * 0.02 * dx 
  shift2 = (i % 4 - i % 2) * 0.02 * dy
  plt.text(xval - shift1, yval - shift2, i)

xrange = (min(x) - 0.1*dx, max(x) + 0.1*dx)
plt.xlim(*xrange)
plt.grid()
plt.legend()
plt.xlabel(f"total cylinder charge {unit}")
plt.ylabel(f"ratio between {name}s \n"
           f"on cylinder tip and middle {symbol("tip")}"
           f"/{symbol("mid")}")
plt.title(title)
plt.savefig(f"fit_parameters_charge_ratio_{unit_identifier}.png",
            bbox_inches="tight")
plt.close()

# ______________________________________________________________________________
# plot charge difference on tip - mid vs total charge:
# ______________________________________________________________________________
x = A_total
y = A_tip - A_mid
dx = max(x) - min(x)
dy = max(y) - min(y)
plt.scatter(x[spherocylinders], y[spherocylinders], c="C1", label="spherocylinder")
plt.scatter(x[cylinders], y[cylinders], c="C0", label="cylinder")
for direc,xval,yval in zip(folders,x,y):
  i = int(direc[-2:])
  shift1 = i % 2 * 0.02 * dx 
  shift2 = (i % 4 - i % 2) * 0.02 * dy
  plt.text(xval - shift1, yval - shift2, i)

xrange = (min(x) - 0.1*dx, max(x) + 0.1*dx)
plt.xlim(*xrange)
plt.grid()
plt.legend()
plt.xlabel(f"total cylinder charge {unit}")
plt.ylabel(f"difference between {name}s\n"
           f"on cylinder tip and mid {symbol("tip")} "
           f"- {symbol("mid")} {unit}")
plt.title(title)
plt.savefig(f"fit_parameters_charge_offset_{unit_identifier}.png",
            bbox_inches="tight")
plt.close()

# ______________________________________________________________________________
# plot screening length vs total charge:
# ______________________________________________________________________________
x = A_total
y = screening
dx = max(x) - min(x)
dy = max(y) - min(y)
plt.scatter(x[spherocylinders], y[spherocylinders], c="C1", label="spherocylinder")
plt.scatter(x[cylinders], y[cylinders], c="C0", label="cylinder")
for direc,xval,yval in zip(folders,x,y):
  i = int(direc[-2:])
  shift1 = i % 2 * 0.02 * dx 
  shift2 = (i % 4 - i % 2) * 0.02 * dy
  plt.text(xval - shift1, yval - shift2, i)

xrange = (min(x) - 0.1*dx, max(x) + 0.1*dx)
plt.xlim(*xrange)
plt.grid()
plt.legend()
plt.xlabel(f"total cylinder charge {unit}")
plt.ylabel(r"screening length $\lambda$ [nm]")
plt.title(title)
plt.savefig(f"fit_parameters_screening_{unit_identifier}.png",
            bbox_inches="tight")
plt.close()
