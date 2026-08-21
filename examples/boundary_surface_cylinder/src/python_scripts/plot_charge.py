# SPDX-FileCopyrightText: 2026 Fabienne Dressler <fab.dressler@web.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Plot average colloid charge (grouped by surface region: tip / mantle) for
# different runs.

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
folders = sys.argv[1:]
electrolyte_default = "full"
if electrolyte_default == "vacuum":
    title = "Charge on cylinder surface in vacuum"
else:
    title = "Charge on cylinder surface in electrolyte"
charge_total = []
charge_cap = []
charge_mantle = []
spherocylinders = []
cylinders = []
valid = []
# ______________________________________________________________________________
# Extract data
# ______________________________________________________________________________
for i,direc in enumerate(folders):
  print(direc)
  charge_total.append(dat.get_charges(direc, region="total"))
  charge_cap.append(dat.get_charges(direc, region="cap", density=True))
  charge_mantle.append(dat.get_charges(direc, region="mantle", density=True))
  if abs(charge_total[-1]) == np.inf or abs(charge_cap[-1]) == np.inf or abs(charge_mantle[-1]) == np.inf:
      valid.append(False)
      print(f"Warning: Infinite charge value in {direc}")
      continue
  else:
      valid.append(True)
      cap_type = dat.get_cap_type(direc)
      if cap_type == "spherical":
          spherocylinders.append(True)
          cylinders.append(False)
      elif cap_type == "flat":
          cylinders.append(True)
          spherocylinders.append(False)
charge_cap = np.array(charge_cap)[valid]
charge_mantle = np.array(charge_mantle)[valid]
charge_total = np.array(charge_total)[valid]
folders = np.array(folders)[valid]

# ______________________________________________________________________________
# plot charge on cap vs mantle:
# ______________________________________________________________________________
x = charge_mantle
y = charge_cap
dx = max(x) - min(x)
dy = max(y) - min(y)
print(dx, dy)
plt.scatter(x[spherocylinders], y[spherocylinders], c="C1", label="spherocylinder")
plt.scatter(x[cylinders], y[cylinders], c="C0", label="cylinder")
for direc,xval,yval in zip(folders,x,y):
  i = int(direc[-2:])
  shift1 = i % 2 * 0.02 * dx 
  shift2 = (i % 4 - i % 2) * 0.02 * dy
  plt.text(xval + shift1, yval - shift2, i, ha='right')

xrange = (min(x) - 0.1*dx, max(x) + 0.1*dx)
c = np.linspace(*xrange)
plt.plot(c,c, 'k',lw=0.3)
plt.xlim(*xrange)
plt.grid()
plt.legend()
plt.xlabel(r"average charge density on cylinder mantle [e/nm²]")
plt.ylabel("average charge density on cylinder cap [e/nm²]")
plt.title(title)
plt.show()

# ______________________________________________________________________________
# plot charge ratio on cap / mantle vs total charge:
# ______________________________________________________________________________
x = charge_total
y = charge_cap / charge_mantle
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
plt.xlabel("total cylinder charge [e]")
plt.ylabel("ratio between average charge densities\n"+
           r"on cylinder cap and mantle $\langle \sigma \rangle_\text{cap}$"
           r"/$\langle \sigma \rangle_\text{mantle}$")
plt.title(title)
plt.show()

# ______________________________________________________________________________
# plot charge difference on cap - mantle vs total charge:
# ______________________________________________________________________________
x = charge_total
y = charge_cap - charge_mantle
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
plt.xlabel("total cylinder charge [e]")
plt.ylabel("difference between average charge densities\n"+
           r"on cylinder cap and mantle $\langle \sigma \rangle_\text{cap}$ "
           r"- $\langle \sigma \rangle_\text{mantle}$ [e/nm²]")
plt.title(title)
plt.show()

