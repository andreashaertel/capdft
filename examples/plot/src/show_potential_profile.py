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
  print("Usage: give a data file via\n"
        f"     {arg[0]} <filename> Z=<columns>\n")
  exit()
files = [arg[1]]
X = 0
Y = 2
Z = 3
zrange = None
save = None
plane_pos = None
for a in arg[2:]:
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

if 1:
  x_mesh = []
  z_mesh = []
  total_ES_xz = []
  for filename in files:
    x,z,phi = DAT.load_slice_plane(filename, X,Y,[Z], plane_pos)
    x_mesh.append(x)
    z_mesh.append(z)
    total_ES_xz.append(phi[0])
  print("density profiles (structured, heatmaps)")
  for d in range(len(files)):
    plt.figure(figsize=figsize)
    POT.heatmap(x_mesh[d],z_mesh[d],total_ES_xz[d],crange=zrange,labels=[label(X), label(Y), label(Z)])
    if plane_pos is not None:
      projected_axis = label((X+Y)%3)
      plt.title(f"${projected_axis} = {plane_pos} \cdot L_{projected_axis}$")
    if save is not None:
      plt.savefig(f"{d}_{save}.pdf", bbox_inches='tight')
    plt.show()
