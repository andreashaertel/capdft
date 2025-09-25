# Source code to produce the potential-dependence plots in the thesis
# (electrode charge density / differential capacitance ...)

import data_structure as STR
import electrode_charge_density as ECD
import load_data as DAT
import plot_basics as PLOT
import numpy as np
import matplotlib.pyplot as plt

plots = "thesis/" # "current_plots/"
figsize = (8,5)

main = STR.data_set(*STR.get_reduced_dict(STR.main_data, 
            ["potential"]).transpose())
planar = STR.data_set(*STR.get_reduced_dict(STR.planar_data, 
            ["potential"]).transpose())

methods = ECD.main_ECD_methods
for dat in [main]:
  dat.get_mean_ECD(methods)
methods = ECD.planar_ECD_methods
for dat in [planar]:
  dat.get_mean_ECD(methods)

# mean ECD
if 1:
  for side in ["left","right"]:
    plt.figure(figsize=figsize)
    i = 0
    methods = ECD.main_ECD_methods
    for m in range(len(methods)):
      for dat in [main]:
        ecd = dat.mean_ECD(side)[m]
        plt.plot(dat.potentials, ecd, f"C{i}"+DAT.markers[i], 
                  label="structured electrodes")
        i += 1
    PLOT.layout_mean_ECD()
    plt.title(f"mean charge density on the {side} electrode")
    plt.savefig(plots+f"electrode_charge/potential_dependence/undulated/macro_{side}.pdf")
    methods = ECD.planar_ECD_methods
    for m in range(len(methods)):
      for dat in [planar]:
        ecd = dat.mean_ECD(side)[m]
        plt.plot(dat.potentials, ecd, f"C{i}"+DAT.markers[i], 
                  label="planar electrodes")
        i += 1
    plt.legend()
    plt.savefig(plots+f"electrode_charge/potential_dependence/both_structures/macro_{side}.pdf")
    #plt.show()

# differential capacitance
if 1:
  for side in ["left","right"]:
    plt.figure(figsize=figsize)
    i = 0
    ECD.plot_c_diff(main.potentials, main.mean_ECD(side), labels=["structured electrodes"], colors=[f"C{i}"+DAT.markers[i]])
    plt.savefig(plots+f"differential_capacitance/undulated/macro_{side}.pdf")
    i += 1
    ECD.plot_c_diff(planar.potentials, planar.mean_ECD(side), labels=["planar electrodes"], colors=[f"C{i}"+DAT.markers[i]])
    plt.legend()
    plt.savefig(plots+f"differential_capacitance/both_structures/macro_{side}.pdf")
    #plt.show()

methods = ECD.main_ECD_methods
for dat in [main]:
  dat.get_local_ECD(methods)
