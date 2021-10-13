import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import pandas as pd
import math as ma



xmax = input('Choose xmax: ')
NN = input('Choose the number of spacings: ')
xmax = int(xmax)
NN = int(NN)
xmin = -xmax

np.savetxt('input.dat',np.c_[xmax,NN],fmt='%u')
subprocess.call('./Ex6-Maniscalco.out')
eigenvalues = np.genfromtxt('eigenvalues.dat')

f0 = plt.figure(0)
f1 = plt.figure(1)
f2 = plt.figure(2)
f3 = plt.figure(3)

hbar = 1
omega = 1
m = 0.5


def eigen0(x):
    return(((2*m*omega)/(ma.pi*hbar))**(0.25)*ma.exp(-m*omega/hbar*x**2))

def eigen1(x):
    return((((2*m*omega)/(ma.pi*hbar))**(0.25)*ma.exp(-m*omega/hbar*x**2)*x*ma.sqrt(4*m*omega/hbar)))

def eigen2(x):
    return(-(((m*omega)/(2*ma.pi*hbar))**(0.25)*ma.exp(-m*omega/hbar*x**2)*(4*m*omega*x**2/hbar-1)))

vect0 = np.vectorize(eigen0)
vect1 = np.vectorize(eigen1)
vect2 = np.vectorize(eigen2)
dom = np.arange(xmin, xmax, 0.1)

plt.figure(0)
plt.plot(eigenvalues[:,0],eigenvalues[:,1],label='Computed')
plt.plot(eigenvalues[:,0],eigenvalues[:,0]*2-1,label='Theoretical')
plt.xlabel('n')
plt.ylabel('f(n)')
plt.title('Eigenvalues')
plt.legend(loc=2)
plt.savefig('figures/N='+str(NN)+'_L='+str(xmax))

print('Mean absolute error of eigenvalues: ',sum(abs((eigenvalues[:,0]*2-1)-eigenvalues[:,1]))/len(eigenvalues))

eigenvectors = np.genfromtxt('eigenvectors.dat')

plt.figure(1)
plt.plot(eigenvectors[:,0],eigenvectors[:,1],label='Computed')
plt.plot(dom,vect0(dom),label='Theoretical')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('First eigenfunction')
plt.legend()
plt.savefig('figures/1f_'+'N='+str(NN)+'_L='+str(xmax))

plt.figure(2)
plt.plot(eigenvectors[:,0],eigenvectors[:,2],label='Computed')
plt.plot(dom,vect1(dom),label='Theoretical')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Second eigenfunction')
plt.legend()
plt.savefig('figures/2f_'+'N='+str(NN)+'_L='+str(xmax))

plt.figure(3)
plt.plot(eigenvectors[:,0],eigenvectors[:,3],label='Computed')
plt.plot(dom,vect2(dom),label='Theoretical')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Third eigenfunction')
plt.legend()
plt.savefig('figures/3f_'+'N='+str(NN)+'_L='+str(xmax))


