# Show differences between two data sets
import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import potentials as POT
import numpy as np
import matplotlib.pyplot as plt
import sys

plots = "" # "current_plots/"
figsize = (8,5)
d_range = (0.7,1.4)

# process command line arguments
arg = sys.argv
if len(arg) == 1:
  print("Usage: give two data files via\n"
        f"     {arg[0]} <filename 1> <filename 2>\n")
  exit()
files = [arg[1], arg[2]]
X = 0
Y = 2
Z = 3
zrange = None
save = None
plane_pos = None
for a in arg[3:]:
  if a[:2] == "Z=": Z = eval(a[2:])
  elif a[:2] == "X=": X = eval(a[2:])
  elif a[:2] == "Y=": Y = eval(a[2:])
  elif a[:6] == "range=": zrange = eval(a[6:])
  elif a[:5] == "save=": save = a[5:]
  elif a[:10] == "plane_pos=": plane_pos = eval(a[10:])
  else: print("invalid argument ", a)
  
def label(index):
  if index == 0: return 'x'
  if index == 1: return 'y'
  if index == 2: return 'z'
  else: return "value (column "+str(index)+")"

x_mesh = []
z_mesh = []
total_ES_xz = []
for filename in files:
  x,z,phi = DAT.load_slice_plane(filename, X,Y,[Z], plane_pos)
  x_mesh.append(x)
  z_mesh.append(z)
  total_ES_xz.append(phi[0])
if 1: 
  # plot data seperately
  for d in range(len(files)):
    plt.figure(figsize=figsize)
    POT.heatmap(x_mesh[d],z_mesh[d],total_ES_xz[d],crange=zrange,labels=[label(X), label(Y), label(Z)])
    if save is not None:
      plt.savefig(f"{d}_{save}", bbox_inches='tight')
    plt.show()
# check compability
if np.any(x_mesh[0] != x_mesh[1]) or np.any(z_mesh[0] != z_mesh[1]):
  print(arg[0], ": data grids not compatible")
  exit(1)
# plot data difference
if 1:
  plt.figure(figsize=figsize)
  POT.heatmap(x_mesh[0],z_mesh[0],total_ES_xz[1] - total_ES_xz[1],crange=zrange,labels=[label(X), label(Y), label(Z)])
  if save is not None:
    plt.savefig(f"diff_{save}.pdf", bbox_inches='tight')
    plt.show()
