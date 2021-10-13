#IMPORTING PACKAGES
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import itertools
from scipy.optimize import curve_fit
from functools import reduce
import math as ma
import subprocess

#USER DECISION
def write(inputs):
      np.savetxt('inputdim.dat',inputs,fmt='%u')

inputs = np.array([3001])
inputs = inputs.reshape(len(inputs),)

print('choose one of the following: [a/b/c]')
print('a) Only generate \'inputdim.dat\' file with inputs: ',inputs,'\n')
print('b) Generate all the input files and run the program \'Ex5-Maniscalco.out\' for histograms and fits \n')
print('c) Only make the plots (you must already have all the files in the FITS folder)\n')
   
decision=input()
if decision == 'a':
    write(inputs)
elif decision == 'b':
    print('Processing... \n')
    write(inputs)
    subprocess.call('./Ex5-Maniscalco.out')   #this automatically launches the program	
    subprocess.call('gnuplot Ex5-Maniscalco-GNUPLOT', cwd='./FITS',shell=True)
    
elif decision == 'c':
    subprocess.call('gnuplot Ex5-Maniscalco-GNUPLOT', cwd='./FITS',shell=True)
else: print('ERROR: got ',decision,' invalid statement.')
    
print('done.')
