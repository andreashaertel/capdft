# Produce all the single-potential-plots in ../plots/thesis
# (density profiles, potential profiles...)

#______________________________________________________________________________
# Imports, settings, definitions
import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import numpy as np
import matplotlib.pyplot as plt
import potentials as POT
import sys

# plot settings
plots = "thesis/"
figsize = (8,5)
def marker(i):
  if i==-1: return '+', 'k'
  return '+', f"C{i}"

# physical parameters
species = [0,1]
d_range = (0.6,1.4) # range densities

e = 1.602176634e-19 # elementary charge [C]

# process command line arguments
arg = sys.argv
if len(arg) == 1:
  print("Usage: give one or multiple potential values via\n"
        f"   {arg[0]} <potentials>\n"
        f"or {arg[0]} <potentials (undulated)> <potentials (planar)>")
  exit()
if len(arg) == 2:
  potentials_m = eval(arg[1])
  potentials_p = eval(arg[1])
elif len(arg) == 3:
  potentials_m = eval(arg[1])
  potentials_p = eval(arg[2])
else:
  print("too many command line arguments...")
  exit()
if type(potentials_m) == float:
  potentials_m = [potentials_m]
if type(potentials_p) == float:
  potentials_p = [potentials_p]

#______________________________________________________________________________
# data 
sample = STR.data_set(*STR.potential_data(potentials_m).transpose(), potentials_m)
planar_sample = STR.data_set(*STR.potential_data(potentials_p, base=STR.planar_data).transpose(), potentials_p)

# load data from files: densities, potentials, electrode charge densities (ECD)
for dat in [sample]:
    dat.get_densities_xz()
    dat.get_total_ES_xz()
for dat in [planar_sample]:
    dat.get_densities_z()
    dat.get_total_ES_z()

# check that all files have been found
if sample.len != len(potentials_m): print("warning sample")
if planar_sample.len != len(potentials_p): print("warning planar sample")

#______________________________________________________________________________
# Plots
print("Plot:")

# density profiles (structured, heatmaps)
if 1:
  print("density profiles (structured, heatmaps)")
  dat = sample
  for nu,drange in zip(species,[(0.6,1.2),(0.1,1.7)]):
    for d in range(dat.len):
      plt.figure(figsize=figsize)
      POT.heatmap_xz(dat.x_mesh[d],dat.z_mesh[d],dat.densities_xz[nu][d],crange=drange,clabel=r"$\rho$ [nm$^{-3}$]")
      plt.savefig(plots+f"densities/undulated/heatmap/nu{nu}_{dat.potentials[d]}V.pdf", bbox_inches='tight')

# density profiles (z direction, with both species in one plot)
if 1:
  print("density profiles (structured, z direction)")
  plt.figure(figsize=figsize)
  for pos in [0,0.5]:
    for dat in [sample]:
      for d in range(dat.len):
        i = 0
        for nu in species:
          z = dat.z(d)
          rho = dat.densities_z(nu,pos)[d]
          pot = dat.potentials[d]
          plt.plot(z, rho, marker(i)[0], color=marker(i)[1],
                  label=f"species {nu+1}")
          i+=1
        PLOT.layout_densities_z()
        plt.legend()
        plt.ylim(*d_range)
        plt.title(f"$\psi^0=${dat.potentials[d]}V - density profiles "
                f"at $x = {pos} L_x$")
        plt.savefig(plots+f"densities/undulated/slice_z/comp_speciies_xpos{pos}_{dat.potentials[d]}V.pdf")
        #plt.show()
if 1:
  print("density profiles (planar, z direction)")
  plt.figure(figsize=figsize)
  for dat in [planar_sample]:
    for d in range(dat.len):
      i = 0
      for nu in species:
        z = dat.z(d)
        rho = dat.densities_z(nu,0)[d]
        pot = dat.potentials[d]
        plt.plot(z, rho, marker(i)[0], color=marker(i)[1],
                label=f"species {nu+1}")
        i+=1
      PLOT.layout_densities_z()
      plt.legend()
      plt.ylim(*d_range)
      plt.title(f"planar electrodes - density profile at $\psi^0=${dat.potentials[d]}V")
      plt.savefig(plots+f"densities/planar/comp_species_{dat.potentials[d]}V.pdf")
      #plt.show()
