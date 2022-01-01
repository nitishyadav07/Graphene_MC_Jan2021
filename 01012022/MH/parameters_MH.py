#########################################################################################################################
###########################################                           ###################################################
########################################### PARAMETERS FOR M-H plots: ###################################################
###########################################                           ###################################################
#########################################################################################################################

Nitt = 50000  				# total number of Monte Carlo steps
warm = 1000    				# Number of warmup steps
measure=100     			# How often to take a measurement
############################ physical constants: ##########################
k_B = 1.381e-23				# Boltzmann constant (Joules/Kelvin)
mu_0 = 1.25663706212e-6		# magnetic permeability of vacuum (H/m)
mu_B = 9.274e-24 			# Bohr magneton (Joules/Tesla)
#g = 16*2.002 					# spin g-factor for superparamag model
g = 2.002 					# spin g-factor for free electrons (for ferromag electrons)
h_bar = 1.054e-34			# reduced Planck's constant
m_e = 9.109e-31				# electron rest mass (kg)
N_Avo = 6.022e23			# Avogadro number
#################  TEMPERATURE and ARRAY SIZE: ##############################
T = 5
#N_2 = 200	      		# linear dimension of the lattice components, lattice-size= N_2 x N_2 (for both edge1 and edge2)
#N = 0


'''
						# N_2 = 50 for SPM case, 200 for ferro/ferri-mag case
'''
######### parameters for SUPERPARAMAGNETIC RUN: ######################
#J1 = 1.5
#J2 = 0
######### parameters for FERRO-/FERRIMAGNETIC RUN: ##########################

#J1_2 = 0.01
#J2_2 = 0.005			
'''
					# Good values are (for each case): 
									#J1_superparamag = 1.7
									#J2_superparamag = 0.0
									#J1_ferromag = 0.4
									#J2_ferromag = 0.003
									#J1_ferrimag = 0.4
									#J2_ferrimag = -0.001

									#N_superparamag = 50
									#N_ferromag = 200
									#N_ferrimag = 200
						# for sum of three types (for BVDGC fitting, SPM+ferro+ferri)
									#fa = 0.01
									#fb = 0.1
									#fc = 1-fa-fb		
									#divisor = 5
'''
######################################################################
J3_low = -31		#-31		# plus/minus3.13 Tesla, unit of J3 (magnetic field) used is Tesla,  B = 3.135 T on either side
J3_high = 31		#31
J3_int = 1
J3_divisor=10		# to make data more resolved, we divide the interval of the list/array J3 by this number, 10, 100, or so.
J3_multiplier = 1			# for superparamagnetic calculations, is 1 for ferrimagnetic/ferromagnetic calculationss
######################################################################
#f1 = 0.005			# fraction of signal from each measurement to be added to get final signal
#f2 = 1-f1
#divisor = 4			# factor to divide final graph values to obtain perfect fit with good y-axis match of experimental plot (because not all of the material in the real/synthesized carbon is graphitic/graphenic. Infact, this factor gives and idea about the fraction of material that is graphenic  (in the current case, 30May2021, BVDGC, divisor = 4)
######################################################################
_N_rows = 500   	# No. of max rows from top to be selected for plotting
num_skipped = 1  	# no. of rows skipped from top
################## range of axes: #######################
xlim_low =	-3.135
xlim_high = 3.135
ylim_low = -0.011
ylim_high = 0.011




#############################################################################################################
#																											#
#								filename, filepaths, optimized parameters:									#
#																											#
#############################################################################################################

#in the lists below, i = 0 for BVPGC, 1 for aBVPGC, 2 for BVDGC, 3 for aBVDGC case
#i is passed as argument for 'Diamag_Component.py', 'Combined_My_ising_1D_GNP_GNR.py' and 'MH_add_and_plot.py'

file_name = ['superparamag', 'ferromag', 'ferrimag']

## 1 ## Diamag_Component.py -->
filepath =['/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_BVPGC/BVPGC', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_aBVPGC/aBVPGC', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVDGC/BVDGC', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVDGC/aBVDGC']
filepath_subtracted_data = ['/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_BVPGC/Diam_subtr_BVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_aBVPGC/Diam_subtr_aBVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVDGC/Diam_subtr_BVDGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVDGC/Diam_subtr_aBVDGC.csv']


## 2 ##	Combined_My_ising_1D_GNP_GNR.py -->
filepath_deltaMH = ['/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_BVPGC/Diam_subtr_BVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_aBVPGC/Diam_subtr_aBVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVDGC/Diam_subtr_BVDGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVDGC/Diam_subtr_aBVDGC.csv']



## 3 ## MH_add_and_plot.py -->
#Parameters used at present:
J1_superparamag = [0.0, 0.0, 0.0, 0.0]
J2_superparamag = [0.0, 0.0, 0.0, 0.0]
J1_ferromag = [0.2, 0.2, 0.2, 0.2]
J2_ferromag = [0.001, 0.001, 0.001, 0.001]
J1_ferrimag = [0.2, 0.2, 0.2, 0.2]
J2_ferrimag = [-0.001, -0.001, -0.001, -0.001]
N_superparamag = [80, 80, 80, 20]
N_ferromag = [400, 400, 400, 400]
N_ferrimag = [200, 200, 200, 500]
fa = [0.02, 0.0, 0.02, 0.0]
fb = [0.04, 0.001, 0.04, 0.002]
#fc = 1-fa-fb		 # this calc is shifted to MH_add_and_plot.py
divisor = [0.015, 0.034, 0.0105, 0.00032]
filepath_deltaMH=['/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_BVPGC/Diam_subtr_BVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_petal/PPMS_aBVPGC/Diam_subtr_aBVPGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVDGC/Diam_subtr_BVDGC.csv', '/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVDGC/Diam_subtr_aBVDGC.csv']
label_ = ['BVPGC exper at 5K', 'aBVPGC exper at 5K', 'BVDGC exper at 5K', 'aBVDGC exper at 5K']
fig_label = ['BVPGC_5K', 'aBVPGC_5K', 'BVDGC_5K', 'aBVDGC_5K']


#############################################################################################################
#																											#
#												End of code													#
#																											#
#############################################################################################################














################################  EXTRA CODE (NOT USEFUL as on 31st December, 2021)  #########################################


# BVPGC is SPM+Ferro+Ferri
# aBVPGC is pure Ferro + Ferri



######## best paramaters for BVPGC (spm+ferro+ferri):
'''
J1_superparamag = 1.7
J2_superparamag = 0.0
J1_ferromag = 0.3
J2_ferromag = 0.003
J1_ferrimag = 0.4
J2_ferrimag = -0.001
N_superparamag = 50
N_ferromag = 200
N_ferrimag = 200
fa = 0.02
fb = 0.04
fc = 1-fa-fb		
divisor = 5
filepath_MH='/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVPGC/BVPGC5K.dat'
filepath_deltaMH='/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_BVPGC/MH_BVPGC_Diam_subt.csv'
'''


######## best paramaters for aBVPGC (ferro+ferri):
'''
J1_superparamag = 1.7
J2_superparamag = 0.0
J1_ferromag = 0.3
J2_ferromag = 0.0001
J1_ferrimag = 0.4
J2_ferrimag = -0.002
N_superparamag = 50
N_ferromag = 200
N_ferrimag = 200
fa = 0.0
fb = 0.001
fc = 1-fa-fb		
divisor = 1
filepath_MH='/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVPGC/aBVPGC5K.dat'
filepath_deltaMH='/home/nitish/Documents/Work/JNU/JNU_Manuscripts/MC_paper/PPMS/MH_detritus/PPMS_aBVPGC/MH_aBVPGC_Diam_subt.csv'
'''



######## best paramaters for BVDGC:
'''

'''


######## best paramaters for aBVDGC:
'''
J1_superparamag = 1.7
J2_superparamag = 0.0
J1_ferromag = 0.3
J2_ferromag = 0.0001
J1_ferrimag = 0.2
J2_ferrimag = -0.004
N_superparamag = 50
N_ferromag = 200
N_ferrimag = 100
fa = 0.0
fb = 0.001
fc = 1-fa-fb		
divisor = 1
'''




