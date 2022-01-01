FROM FILE: "./Detrit_MH/add_and_plot.py"   -->


######## best paramaters for BVPGC (spm+ferro+ferri):

J1_superparamag = 0.0
J2_superparamag = 0.0
N_superparamag = 80

J1_ferromag = 0.2
J2_ferromag = 0.001
N_ferromag1 = 400
N_ferromag2 = 1000

J1_ferrimag = J1_ferromag
J2_ferrimag = -0.001
N_ferrimag1 = 200
N_ferrimag2 = 500

fspm = 0.1
fferro1 = 0.4
fferro2 = 0.01 #0.03
fferri1 = 0.2
fferri2 = 1-(fspm + fferro1 + fferro2 + fferri1)



######## best paramaters for a-BVPGC (spm+ferro+ferri):

J1_superparamag = 0.0
J2_superparamag = 0.0
N_superparamag = 20

J1_ferromag = 0.2
J2_ferromag = 0.001
N_ferromag1 = 400
N_ferromag2 = 1000

J1_ferrimag = J1_ferromag
J2_ferrimag = -0.001
N_ferrimag1 = 200
N_ferrimag2 = 500

fspm = 0.9
fferro1 = 0.1
fferro2 = 0.0
fferri1 = 0.0
fferri2 = 0.0 #1-(fspm + fferro1 + fferro2 + fferri1)






######## best paramaters for BVDGC (spm+ferro+ferri):

J1_superparamag = 0.0
J2_superparamag = 0.0
N_superparamag = 80

J1_ferromag = 0.2
J2_ferromag = 0.001
N_ferromag1 = 400
N_ferromag2 = 1000

J1_ferrimag = J1_ferromag
J2_ferrimag = -0.001
N_ferrimag1 = 200
N_ferrimag2 = 500

fspm = 0.15
fferro1 = 0.1
fferro2 = 0.03
fferri1 = 0.6
fferri2 = 1-(fspm + fferro1 + fferro2 + fferri1)







######## best paramaters for aBVDGC (ferro+ferri):
J1_superparamag = 0.0
J2_superparamag = 0.0
N_superparamag = 80

J1_ferromag = 0.2
J2_ferromag = 0.001
N_ferromag1 = 400
N_ferromag2 = 1000

J1_ferrimag = J1_ferromag
J2_ferrimag = -0.001
N_ferrimag1 = 200
N_ferrimag2 = 500

fspm = 0.0
fferro1 = 0.01
fferro2 = 0.01
fferri1 = 0.8
fferri2 = 1-(fspm + fferro1 + fferro2 + fferri1)

















-------------------------------------------------------------------------------------------------

FROM FILE: "./Detrit_MH/MH_add_and_plot_SINGLE.py"   -->

BVDGC is SPM+Ferro+Ferri
aBVDGC is pure Ferro + Ferri
'''
######## best paramaters for aBVDGC (ferro+ferri):
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
filepath_MH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVDGC/aBVDGC5K.dat'
filepath_deltaMH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVDGC/MH_aBVDGC_Diam_subt.csv'
'''


Best paramaters for BVDGC (spm+ferro+ferri):
'''
######## best paramaters for aBVPGC (ferro+ferri):
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
filepath_MH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVPGC/aBVPGC5K.dat'
filepath_deltaMH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVPGC/MH_aBVPGC_Diam_subt.csv'
'''

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
filepath_MH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_BVDGC/BVDGC5K.dat'
filepath_deltaMH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_BVDGC/MH_BVDGC_Diam_subt.csv'
'''





---------------------------------------------------------------------------------------------------------------------
Best params for a-BVDGC as seen in file "MH_add_and_plot.py"-->

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

Best params for BVPGC as seen in file "MH_add_and_plot.py"-->

'''
######## best paramaters for BVPGC (spm+ferro+ferri):
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
filepath_MH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_BVPGC/BVPGC5K.dat'
filepath_deltaMH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_BVPGC/MH_BVPGC_Diam_subt.csv'
'''


Best params for a-BVPGC as seen in file "MH_add_and_plot.py"-->

'''
######## best paramaters for aBVPGC (ferro+ferri):
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
filepath_MH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVPGC/aBVPGC5K.dat'
filepath_deltaMH='/home/nitish/Desktop/JNU/Work_/My_manuscripts/Analysis/Magnetization/MH_detritus/PPMS_a-BVPGC/MH_aBVPGC_Diam_subt.csv'
'''
