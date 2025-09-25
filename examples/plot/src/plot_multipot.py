# Source code to produce the plots with multiple potential values in the thesis
# (comparison of potential / density profiles for different applied potentials)

import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import numpy as np
import matplotlib.pyplot as plt
import sys

plots = "thesis/" # "current_plots/"
figsize = (8,5)
species = [0,1]
d_range = (0.7,1.4)

# process command line arguments
arg = sys.argv
if len(arg) == 1:
  print("Usage: give a sample of potential values via\n"
        f"     {arg[0]} -f <file with list of (multiple) samples>\n"
        f"  or {arg[0]} <label> <potentials>\n"
        f"  or {arg[0]} <label> <potentials (undulated)> <potentials (planar)>")
  exit()
samples = []
if arg[1] == "-f":
  file = open(arg[2], "r")
  for line in file:
    if line[0] == "#" or line == "": pass
    else:
      sample = line.split(' ')
      if len(sample) == 3:
        sample = [sample[0], eval(sample[1]), eval(sample[2])]
        samples.append(sample)
        print(sample)
      else:
        print(f"Warning: need three columns in {arg[2]}")
        exit()
elif len(arg) == 3:
  samples.append([arg[1], eval(arg[2]), eval(arg[2])])
elif len(arg) == 4:
  samples.append([arg[1], eval(arg[2]), eval(arg[3])])
else:
  print("too many command line arguments...")
  exit()

for sample_label, potentials_m, potentials_p in samples:
  suffix = "_" + sample_label
  #potentials_m = [0.0,0.1,0.2,0.33, 0.4, 0.55,0.6]
  #potentials_p = [0.0,0.1,0.2,0.3, 0.4, 0.5,0.6]
  colors = plt.cm.viridis(np.linspace(0, 1, len(potentials_m)))
  def marker(i):
    try: return '.', colors[i]
    except: return DAT.markers[i], f"C{i}"
  sample = STR.data_set(*STR.potential_data(potentials_m).transpose(), potentials_m)
  planar_sample = STR.data_set(*STR.potential_data(potentials_p, base=STR.planar_data).transpose(), potentials_p)
  
  for dat in [sample]:
    dat.get_densities_xz()
    dat.get_total_ES_xz()
  for dat in [planar_sample]:
    dat.get_densities_z()
    dat.get_total_ES_z()
  
  # density profiles
  if 1:
    plt.figure(figsize=figsize)
    for pos in [0,0.5]:
      for nu in species:
        i = 0
        for dat in [sample]:
          for d in range(dat.len):
            z = dat.z(d)
            rho = dat.densities_z(nu,pos)[d]
            pot = dat.potentials[d]
            plt.plot(z, rho, marker(i)[0], color=marker(i)[1], label=pot)
            i+=1
        PLOT.layout_densities_z()
        plt.legend(title='$\psi^0$ [V]')
        plt.ylim(*d_range)
        plt.title(f"density profile of species {nu+1} "
                f"at $x = {pos} L_x$")
        plt.savefig(plots+f"densities/undulated/slice_z/comp_pot_xpos{pos}_nu{nu}"+suffix+".pdf")
        #plt.show()
  if 1:
    plt.figure(figsize=figsize)
    for nu in species:
      i = 0
      for dat in [planar_sample]:
        for d in range(dat.len):
          z = dat.z(d)
          rho = dat.densities_z(nu,0)[d]
          pot = dat.potentials[d]
          plt.plot(z, rho, marker(i)[0], color=marker(i)[1], label=pot)
          i+=1
      PLOT.layout_densities_z()
      plt.legend(title='$\psi^0$ [V]')
      plt.ylim(*d_range)
      plt.title(f"planar electrodes - density profile of species {nu+1}")
      plt.savefig(plots+f"densities/planar/comp_pot_nu{nu}"+suffix+".pdf")
      #plt.show()
  
  
  # potential profiles
  if 1:
    plt.figure(figsize=figsize)
    for pos in [0,0.5]:
      i = 0
      for dat in [sample]:
        for d in range(dat.len):
          z = dat.z(d)
          vals = dat.total_ES_z(pos)[d]
          pot = dat.potentials[d]
          plt.plot(z, abs(vals), marker(i)[0], color=marker(i)[1], label=pot)
          if pot > 0 and dat.len <= 3: PLOT.layout_total_ES_z(True, (PLOT.f(pos),8-PLOT.f(pos)), PLOT.volt_reduce_unit(pot), num=i)
          elif pot > 0: PLOT.layout_total_ES_z(False)
          i+=1
      plt.legend(title='$\psi^0$ [V]')
      plt.title(f"electrostatic potential "
              f"at $x = {pos} L_x$")
      plt.savefig(plots+f"potential/undulated/comp_pot_xpos{pos}"+suffix+".pdf")
      #plt.show()
  if 1:
    plt.figure(figsize=figsize)
    i = 0
    for dat in [planar_sample]:
      for d in range(dat.len):
        z = dat.z(d)
        vals = np.array(dat.total_ES_z(0)[d])
        pot = dat.potentials[d]
        plt.plot(z, abs(vals), marker(i)[0], color=marker(i)[1], label=pot)
        if pot > 0 and dat.len <= 3: PLOT.layout_total_ES_z(True, (0,6.4), PLOT.volt_reduce_unit(pot), num=i)
        elif pot > 0: PLOT.layout_total_ES_z(False)
        i+=1
    plt.legend(title='$\psi^0$ [V]')
    plt.title(f"planar electrodes - electrostatic potential")
    plt.savefig(plots+"potential/planar/comp_pot"+suffix+".pdf")
    #plt.show()
  
  methods = ECD.main_ECD_methods
  for dat in [sample]:
    dat.get_local_ECD(methods)
  
  # local ECD
  if 1:
    plt.figure(figsize=figsize)
    for side in ["left","right"]:
      i = 0
      data = [sample]
      methods = ECD.main_ECD_methods
      for m in range(len(methods)):
        for dat in data:
          for d in range(dat.len):
            ecd = dat.local_ECD(side)[m][d]
            pot = dat.potentials[d]
            plt.plot(dat.x(d), ecd, marker(i)[0]+'-', color=marker(i)[1], label=pot)
            i += 1
      PLOT.layout_local_ECD()
      plt.legend(title="$\psi^0$ [V]")
      plt.title(f"charge distribution on the {side} electrode")
      plt.savefig(plots+f"electrode_charge/local_distribution/comp_pot_{side}"+suffix+".pdf")
      #plt.show()
  
  # local ECD - change with potential
  if 1:
    plt.figure(figsize=figsize)
    for side in ["left","right"]:
      data = [sample]
      methods = ECD.main_ECD_methods
      for m in range(len(methods)):
        for dat in data:
          for d in range(1,dat.len):
            delta_sigma = dat.local_ECD(side)[m][d]-dat.local_ECD(side)[m][d-1]
            delta_pot = dat.potentials[d]-dat.potentials[d-1]
            pot = dat.potentials[d]
            plt.plot(dat.x(d), delta_sigma/delta_pot, 'x-', label=f"{pot}V, {side} el.")
    PLOT.layout_local_ECD()
    plt.ylabel("$d\sigma/d\psi^0$ [e/nm²/V]")
    plt.legend()
  #  plt.title(f"charge distribution on the {side} electrode")
    plt.savefig(plots+f"electrode_charge/local_distribution/derivative"+suffix+".pdf")
    #plt.show()
