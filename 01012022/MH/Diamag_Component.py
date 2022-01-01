import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import csv
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import codecs
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import math
from parameters_MH import *


# USING ARGUMENT PARSER (argparse package) to use arguments (in case of more than 2 arguments passed from terminal to script):
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--i', nargs = 1)
args = parser.parse_args()
i = int(args.i[0])

######################################################################################################################
#########################  	 DATA ACQUISITION AND PLOTTING: 	######################################################
######################################################################################################################

################## 	DATA ACQUISITION from files: 	#########################
num_skipped = 1  # no. of rows skipped from top

x1, y1 = np.loadtxt(filepath[i]+'5K.dat', dtype='double', delimiter='\t', skiprows=num_skipped, unpack=True, usecols=[0, 2])
x2, y2 = np.loadtxt(filepath[i]+'100K.dat', dtype='double', delimiter='\t', skiprows=num_skipped, unpack=True, usecols=[0, 2])
x3, y3 = np.loadtxt(filepath[i]+'300K.dat', dtype='double', delimiter='\t', skiprows=num_skipped, unpack=True, usecols=[0, 2])

#################### finding diamagnetic contribution ###################
x_max = np.max(x1)
x_min = np.min(x1)
y_max = np.max(y1)
y_min = np.min(y1)

for i in range(0,len(x3)):
    if x3[i]==x_min:
        min_ele=x3[i]
        pos_x_min = i
    if x3[i]==x_max:
        max_ele=x3[i]
        pos_x_max = i
print(pos_x_max, pos_x_min)
print(pos_x_max-150, pos_x_min-150)

slope_left = (y3[pos_x_min]-y3[pos_x_min-150])/(x3[pos_x_min]-x3[pos_x_min-150])
slope_right = (y3[pos_x_max]-y3[pos_x_max-150])/(x3[pos_x_max]-x3[pos_x_max-150])
print(slope_left, slope_right)

avg_slope = (slope_left + slope_right)*0.5

y_diam1 = []
y_diam2 = []
y_diam3 = []

for i in range(0,len(x1)):
	y_diam1.append(avg_slope*x1[i])
for i in range(0,len(x2)):
	y_diam2.append(avg_slope*x2[i])
for i in range(0,len(x3)):
	y_diam3.append(avg_slope*x3[i])

################# PLOTTING FROM DATA ACQUIRED from files: #######################
fig, ax = plt.subplots()
#plt.plot(x1, y1, '-', color='red', label='aBVDGC exp at 300K')
#plt.plot(x1, y_diam, '-', color='blue', label='diamagnetic part')
plt.plot(x1, y1-y_diam1, '-', color='black', label='5K data with diam. part removed')
plt.plot(x2, y2-y_diam2, '-', color='blue', label='100K data with diam. part removed')
plt.plot(x3, y3-y_diam3, '-', color='green', label='300K data with diam. part removed')
plt.xlim()
plt.ylim()
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel(r"$B$/T$")
plt.ylabel(r"$M$/Am2 kg$^{-1}}$.")
#plt.ylabel(r"$\Delta$""$M$/emu g$^{-1}}$.")
plt.show()
plt.close()


########################### WRITING DATA TO FILE #############################
f = open(filepath_subtracted_data[i], "w")
f.write("{},{}\n".format("H", "M"))
for x in zip(x1, y1-y_diam1):
	f.write("{},{}\n".format(x[0], x[1]))
f.close()
