import numpy as np
import subprocess
import sys
import os

def write(Nmin,Nmax,step):
    N=np.arange(Nmin,Nmax,step)
    np.savetxt('inputdim.dat',N,fmt='%u')

Nstart = 100
Nend = 501
step = 50

print('Choose one of the following: [a/b/c]:')
print('a) Only generate \'inputdim.dat\' file with (start,stop,step) = (',Nstart,Nend,step,')\n')
print('b) Generate all the input files, runs the program \'Ex4-Maniscalco.out\' for time misurations, makes fits with gnuplot \n')
print('c) Only make the fits (you must already have all the files in the folder)\n')
   
decision=input()
if decision == 'a':
    write(Nstart,Nend,step)
elif decision == 'b':
    print('Processing for (start,stop,step) = (',Nstart,Nend,step,') \n')
    write(Nstart,Nend,step)
    subprocess.call('./Ex4-Maniscalco.out')   #this automatically launches the program	
    os.system('gnuplot Ex4-Maniscalco-gnuplot')
elif decision == 'c':
    os.system('gnuplot Ex4-Maniscalco-gnuplot')
else: print('ERROR: got ',decision,' invalid statement.')



