import matplotlib.pyplot as plt
import numpy as np
import load_data as dat
import data_structure as struct

def eff_range(z, epsilon=1e-4):
  z = np.array(z)
  epsilon *= (z.max() - z.min())
  if abs(z.min()) < epsilon:
    eff_min = z[z > 0].min()
    eff_max = z[z < z.max() - epsilon].max()
  else:
    eff_min = z[z > z.min() + epsilon].min()
    eff_max = z[z < z.max() - epsilon].max()
  print("plot z range: ", eff_min, eff_max) 
  return eff_min, eff_max

def var_range(z, nsigma=2):
  z = np.array(z)
  mean = np.mean(z)
  delta_z = nsigma * np.std(z)
  return max([mean - delta_z, z.min()]), min([mean + delta_z, z.max()])

def heatmap(a,b,c,crange,labels):
  if type(crange) == str and crange == "auto": crange = eff_range([c])
  if crange is None: plt.pcolormesh(a,b,c)
  else: plt.pcolormesh(a,b,c,vmin=crange[0],vmax=crange[1])
  plt.colorbar(label=labels[2])
  plt.xlabel(labels[0])
  plt.ylabel(labels[1])
  
def heatmap_xz(x,z,c,crange=None,clabel=""):
  heatmap(z,x,c,crange,["$z$ [nm]", "$x$ [nm]", clabel])

