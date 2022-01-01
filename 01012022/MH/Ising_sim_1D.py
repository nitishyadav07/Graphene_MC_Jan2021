#!/usr/bin/env python
"""
Monte Carlo simulation of the 1D Ising model

http://www.physics.rutgers.edu/~haule/681/src_MC/python_codes/ising.py
"""
from matplotlib.pyplot import imshow
import numpy as np
from scipy import *
from pylab import *
import math as math
import sys # to write data to file
import csv # for handling csv files, data appending
import os
from os import path		# for checking if data file already exists or not
from numpy import array, savetxt # for saving csv file
from parameters_MH import *	# parameter definitions are stored in this file
from func_mag import *		# function definitions are stored in this file
import sys, getopt
'''
print(str(sys.argv))
print(sys.argv[0]) # prints python_script.py
print(sys.argv[2]) # prints var2
'''
#J1_2 = float(sys.argv[1])
#J2_2 = float(sys.argv[2])

#############################################################################################################
#																											#
#																											#
#																											#
#############################################################################################################
# USING ARGUMENT PARSER (argparse package) to use arguments (in case of more than 2 arguments passed from terminal to script):
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--J1', nargs = 1)
parser.add_argument('--J2', nargs = 1)
parser.add_argument('--N', nargs = 1)
#parser.add_argument('--J3', nargs = 1)
args = parser.parse_args()
#print(args)
J1_2 = float(args.J1[0])
J2_2 = float(args.J2[0])
#J3 = float(args.J3[0])
N_2 = int(args.N[0])
print('[J1, J2, N] = ', [J1_2, J2_2, N_2])
#print('[J1, J2, J3] = ', [J1_2, J2_2, J3])
if N_2 == 80:
	g = 16*2.002
if N_2 == 20:
	g = 4*2.002
#print("J1_2 = ", str(J1_2), "J2_2 = ", str(J2_2), "Nitt = "+str(Nitt), "N = "+str(N_2))
#print("Started Calculation")
#############################################################################################################
#																											#
#																											#
#																											#
#############################################################################################################
if g > 2.002:
	flag_spm = True      	# if 'true' then calculation is to be done for Superparamagnetic case
else:
	flag_spm = False

if J2_2 >= 0:
	FERRO = True
elif J2_2 < 0:
	FERRO = False
FERRI = 1-FERRO
#############################################################################################################
#																											#
#																											#
#																											#
#############################################################################################################
if flag_spm == True:
	system_name = 'superparamag'
	Mag_label = r'M$_{superparamag}$'
	plot_color = ['black', 'red']
	file_flag = 0
elif flag_spm == False and FERRO == True:
	system_name = 'ferromag'
	Mag_label = r'M$_{ferromag}$'
	plot_color = ['blue', 'green']
	file_flag = 1
elif flag_spm == False and FERRI == True:
	system_name = 'ferrimag'
	Mag_label = r'M$_{ferrimag}$'
	plot_color = ['brown', 'navy']
	file_flag = 2
#############################################################################################################
#																											#
#																											#
#																											#
#############################################################################################################
print("results/data_"+system_name+"_T"+str(T)+"_Nitt"+str(Nitt)+"_warmsteps"+str(warm)+"_measure"+str(measure)+"_N_2"+str(N_2)+"_J1_2"+str(J1_2)+"_J2_2"+str(J2_2)+".csv")
if path.exists("results/data_"+system_name+"_T"+str(T)+"_Nitt"+str(Nitt)+"_warmsteps"+str(warm)+"_measure"+str(measure)+"_N_2"+str(N_2)+"_J1_2"+str(J1_2)+"_J2_2"+str(J2_2)+".csv") == False:

	(edge1, edge2) = RandomL(N_2)
	#edge1 = ones((N), dtype=int)
	#edge2 = ones((N), dtype=int)
	pc_status = 0
	pc_multiplier = 1 #multiplies by 2 the total values in case of Ferrimagnetic case
	(wMag, wEne, wfield, wChi) = MH_loop(Nitt, edge1, edge2, J1_2, J2_2, T, J3_divisor, J3_high, J3_low, J3_int, pc_status, pc_multiplier)
	f = open("results/data_"+system_name+"_T"+str(T)+"_Nitt"+str(Nitt)+"_warmsteps"+str(warm)+"_measure"+str(measure)+"_N_2"+str(N_2)+"_J1_2"+str(J1_2)+"_J2_2"+str(J2_2)+".csv", "w")
	f.write("{},{},{},{}\n".format("Field", "Energy", "Magnetization", "Susceptibility"))
	for x in zip(wfield, wEne, wMag, wChi):
		f.write("{},{},{},{}\n".format(x[0], x[1], x[2], x[3]))
	f.close()
	#print([aE, aM])
#############################################################################################################
#																											#
#																											#
#																											#
#############################################################################################################
