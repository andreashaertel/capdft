# Plots for comparison of different numerical resolutions (appendix A in the
# bachelor thesis)
import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import numpy as np
import matplotlib.pyplot as plt

plots = "thesis/" #"current_plots/"
figsize = (8,5)
species = [0,1]
main = STR.data_set(*STR.get_reduced_dict(STR.main_data, 
            ["potential"]).transpose())
main_gc = STR.data_set(*STR.get_reduced_dict(STR.gc_data, 
            ["potential"]).transpose())
example = STR.data_set([STR.example[0]], [STR.example[1]], [0.1])
planar = STR.data_set(*STR.get_reduced_dict(STR.planar_data, 
            ["potential"]).transpose())
#planar_gc = STR.data_set(*STR.get_reduced_dict(STR.planar_gc_data, 
#            ["potential"]).transpose())

for dat in [main_gc, example]:
  dat.get_densities_xz()
  dat.get_total_ES_xz()
  dat.get_gc()

# density profiles
if 1:
  for pos in [0,0.5]:
    for nu in species:
      plt.figure(figsize=figsize)
      i = 0
      for dat in [example, main_gc]:
        for d in range(dat.len):
          z = dat.z(d)
          rho = dat.densities_z(nu,pos)[d]
          gc = dat.gc[d]
          plt.plot(z, rho, f"C{i}"+DAT.markers[i], label=gc)
          i+=1
      PLOT.layout_densities_z()
      plt.legend(title='grid counts')
      plt.title(f"$\psi^0 = 0.1$ V - density profile of species {nu+1} "
              f"at $x = {pos} L_x$")
      plt.savefig(plots+f"densities/undulated/slice_z/comp_gc_xpos{pos}_nu{nu}.pdf")
      #plt.show()

# potential profiles
if 1:
  for pos in [0,0.5]:
    plt.figure(figsize=figsize)
    i = 0
    for dat in [example, main_gc]:
      for d in range(dat.len):
        z = dat.z(d)
        vals = dat.total_ES_z(pos)[d]
        gc = dat.gc[d]
        plt.plot(z, abs(vals), f"C{i}"+DAT.markers[i], label=gc)
        i+=1
    PLOT.layout_total_ES_z()
    plt.legend(title='grid counts')
    plt.title(f"$\psi^0 = 0.1$ V - electrostatic potential "
            f"at $x = {pos} L_x$")
    plt.savefig(plots+f"potential/undulated/comp_gc_xpos{pos}.pdf")
    #plt.show()

methods = ECD.main_ECD_methods
for dat in [example, main_gc]:
  dat.get_local_ECD(methods)

# local ECD
if 1:
  for side in ["left","right"]:
    plt.figure(figsize=figsize)
    i = 0
    for m in range(len(methods)):
      for dat in [example, main_gc]:
        for d in range(dat.len):
          ecd = dat.local_ECD(side)[m][d]
          gc = dat.gc[d]
          plt.plot(dat.x(d), ecd, f"C{i}"+DAT.markers[i]+'-', label=gc)
          i += 1
    PLOT.layout_local_ECD()
    plt.legend(title="grid counts")
    plt.title(f"$\psi^0 = 0.1$ V - charge distribution on the {side} electrode")
    plt.savefig(plots+f"electrode_charge/local_distribution/comp_gc_{side}.pdf")
    #plt.show()

methods = ECD.all_ECD_methods
for dat in [main]:
  dat.get_mean_ECD(methods)
methods = ECD.all_planar_ECD_methods
for dat in [planar]:
  dat.get_mean_ECD(methods)

# mean ECD
if 1:
  for side in ["left","right"]:
    plt.figure(figsize=figsize)
    i = 0
    methods = ECD.all_ECD_methods
    for m in range(len(methods)):
      for dat in [main]:
        ecd = dat.mean_ECD(side)[m]
        plt.plot(dat.potentials, ecd, f"C{i}"+DAT.markers[i], 
                  label=ECD.ECD_methods[methods[m]])
        i += 1
    PLOT.layout_mean_ECD()
    plt.legend(title="Method")
    plt.title(f"mean charge density on the {side} electrode")
    plt.savefig(plots+f"electrode_charge/potential_dependence/both_structures/comp_methods_macro_{side}.pdf")
    #plt.show()

# differential capacitance
if 1:
  for side in ["left","right"]:
    plt.figure(figsize=figsize)
    methods = ECD.all_ECD_methods
    dat = main
    ECD.plot_c_diff(dat.potentials, dat.mean_ECD(side), labels=[f"structured electrodes ({ECD.ECD_methods[m]})" for m in methods], 
            colors=[f"C{i}"+DAT.markers[i] for i in range(len(methods))])
    plt.grid()
    plt.title(f"differential capacitance at the {side} electrode")
    plt.savefig(plots+f"differential_capacitance/undulated/comp_methods_macro_{side}.pdf")
    methods = ECD.all_ECD_methods
    dat = planar
    ECD.plot_c_diff(dat.potentials, dat.mean_ECD(side), labels=["planar electrodes" for m in methods], 
        colors=[f"k"+DAT.markers[i] for i in range(len(methods))])
    plt.legend()
    plt.title(f"differential capacitance at the {side} electrode")
    plt.savefig(plots+f"differential_capacitance/both_structures/comp_methods_macro_{side}.pdf")
    #plt.show()

#methods = ECD.main_ECD_methods
#for dat in [main]:
#  dat.get_local_ECD(methods)

#if 1:
#  for side in ["left","right"]:
#    i = 0
#    for m in range(len(methods)):
#      for dat in [planar]:
#        ecd = dat.mean_ECD(side)[m]
#        plt.plot(dat.potentials, ecd, f"C{i}"+DAT.markers[i], 
#                  label=ECD.ECD_methods[methods[m]])
#        i += 1
#    PLOT.layout_mean_ECD()
#    plt.legend(title="Method")
#    plt.title(f"mean charge density on the {side} (planar) electrode")
#    plt.savefig(plots+f"planar_mECD_comp_methods_{side}.pdf")
#    plt.show()

