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

############################# IMPORT code ################################
'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--J1', nargs = 1)
parser.add_argument('--J2', nargs = 1)
#parser.add_argument('--J3', nargs = 1)
args = parser.parse_args()
#print(args)
J1_2 = float(args.J1[0])
J2_2 = float(args.J2[0])
#J3 = float(args.J3[0])
print('[J1, J2] = ', [J1_2, J2_2])
#print('[J1, J2, J3] = ', [J1_2, J2_2, J3])
'''
##########################################################################
# USING ARGUMENT PARSER (argparse package) to use arguments (in case of more than 2 arguments passed from terminal to script):
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--i', nargs = 1)
args = parser.parse_args()
i = int(args.i[0])
#########################################################################

fc = np.subtract(1, np.subtract(fa, fb))		# fa and fb are imported from parameters_MH and then used to calc fc


######################################################################################################################
#########################  	 DATA ACQUISITION AND PLOTTING: 	######################################################
######################################################################################################################
mu_0 = 4*math.pi*10e-6

################## 	DATA ACQUISITION from files: 	#########################
num_skipped = 1  # no. of rows skipped from top

Ba, Ma = np.loadtxt('./results/data_'+file_name[0]+'_T'+'5'+'_Nitt'+str(Nitt)+'_warmsteps'+str(warm)+'_measure'+str(measure)+'_N_2'+str(N_superparamag[i])+'_J1_2'+str(J1_superparamag[i])+'_J2_2'+str(J2_superparamag[i])+'.csv', dtype='double', delimiter=',', skiprows=num_skipped, unpack=True, usecols=[0, 2])
Bb, Mb = np.loadtxt('./results/data_'+file_name[1]+'_T'+'5'+'_Nitt'+str(Nitt)+'_warmsteps'+str(warm)+'_measure'+str(measure)+'_N_2'+str(N_ferromag[i])+'_J1_2'+str(J1_ferromag[i])+'_J2_2'+str(J2_ferromag[i])+'.csv', dtype='double', delimiter=',', skiprows=num_skipped, unpack=True, usecols=[0, 2])
Bc, Mc = np.loadtxt('./results/data_'+file_name[2]+'_T'+'5'+'_Nitt'+str(Nitt)+'_warmsteps'+str(warm)+'_measure'+str(measure)+'_N_2'+str(N_ferrimag[i])+'_J1_2'+str(J1_ferrimag[i])+'_J2_2'+str(J2_ferrimag[i])+'.csv', dtype='double', delimiter=',', skiprows=num_skipped, unpack=True, usecols=[0, 2])

Ha = (Ba/mu_0) - Ma
Hb = (Bb/mu_0) - Mb
Hc = (Bc/mu_0) - Mc

##########################		 M-H data				#########################
#x, y = np.loadtxt(filepath_MH, dtype='double', delimiter=',', skiprows=1, unpack=True, usecols=[0, 2])

##########################  	(DELTA)M-H data			   ######################
B, delM = np.loadtxt(filepath_deltaMH[i], dtype='double', delimiter=',', skiprows=1, unpack=True, usecols=[0, 1])
## convert B to H:
H = (B/mu_0) - delM   # units of H are Am-1   1 Am-1 = 0.01256637061436 Oe

################# PLOTTING FROM DATA ACQUIRED from files: #######################
fig, ax = plt.subplots()

plt.plot(B, delM, '-', color='black', label=label_)
plt.plot(Ba, ((fa[i]*Ma*mu_B*N_Avo/(N_superparamag[i]*N_superparamag[i]))+(fb[i]*Mb*mu_B*N_Avo/(N_ferromag[i]*N_ferromag[i]))+(fc[i]*Mc*mu_B*N_Avo/(N_ferrimag[i]*N_ferrimag[i])))/divisor[i], '-', color='red', label='sum of signals')
plt.xlim()
plt.ylim()
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.legend()
plt.xlabel(r"$B$/T$")
plt.ylabel(r"$\Delta$""$M$/Am2 kg$^{-1}}$.")
fig.savefig("./results/"+ fig_label[i] +"_summed_plot"+"_Nitt"+str(Nitt)+"_"+str(warm)+"_"+str(measure)+"_spmfr_"+str(fa[i])+"fmfr_"+str(fb[i])+"fimfr_"+str(fc[i])+"_"+str(J1_superparamag[i])+"_"+str(J2_superparamag[i])+"_"+str(J1_ferromag[i])+"_"+str(J2_ferromag[i])+"_"+str(J1_ferrimag[i])+"_"+str(J2_ferrimag[i])+"_"+str(N_superparamag[i])+"_"+str(N_ferromag[i])+"_"+str(N_ferrimag[i])+".png", format="png",dpi=600)
plt.show()
plt.close()

