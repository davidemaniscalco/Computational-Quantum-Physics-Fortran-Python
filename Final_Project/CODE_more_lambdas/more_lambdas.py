#############################################################################
# FINAL EXERCISE INFORMATION THEORY AND COMPUTATION COURSE
# Master degree in Physics of Data, University of Padova, January 2020
# Davide Maniscalco
#--------------------------------------------------------------------------
# PYTHON SCRPIT FOR AUTOMATION
# -------------------------------------------------------------------------
# This script is essentially made of one loop that at each iteration generates two
#  values of lambda, save them, launch the fortran code for the computation of epsilon,
#  and the gap; then calculates the estimator. All the results (one estimator for each
#  pair of lambdas) are stored in the "total_vec.dat" file
#####################################################################################


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import subprocess

aa = np.linspace(-1,1,21)
dd = zip(np.repeat(aa,21),np.tile(aa,21))
rr = np.array(list(dd))

nn = 3
nsteps = 1000

total_vector = np.zeros((len(aa))**2)
gap = np.zeros((nsteps+1,2**nn-1))

for qq in range((len(aa))**2):
    print(qq)
    np.savetxt('input_list.dat', rr[qq].reshape(1,2), fmt='%1.1f')
    subprocess.call('./Exfinale.out')
    spectrum = np.genfromtxt('true_spectrum.dat') #g.d.s. "spectrum"
    epsilon_gs = np.genfromtxt('epsilon_gs.dat')  #g.d.s. "true_evol"
    
    spectrum_notime = spectrum[:,1:]
    epsilon_notime = epsilon_gs[:,1:]
    
    for ii in range(0,2**nn-1): gap[:,ii] = spectrum_notime[:,ii+1] - spectrum_notime[:,0]
    gap[gap < 1e-6] = 0.0
    epsilon_notime_nozero = epsilon_notime[:,1:]
    
    estimator = epsilon_notime_nozero/gap**2
    x_index = int(np.nanargmax(estimator[estimator != np.inf]) // (2**nn-1))
    y_index = int(np.nanargmax(estimator[estimator != np.inf]) % (2**nn-1))
    
    total_vector[qq] = estimator[x_index, y_index]
    print(total_vector[qq])
print(total_vector)
np.savetxt('totalvec.dat',total_vector)
