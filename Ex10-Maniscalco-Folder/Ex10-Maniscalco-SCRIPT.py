import numpy as np
import matplotlib.pyplot as plt
import subprocess

# lambda generation
lambda_min = float(input('Choose the minimum value of lambda: '))
lambda_max = float(input('Choose the maximum value of lambda: '))
num_lambda = int(input('Choose the number of lambdas: '))
lambdas = np.linspace(lambda_min,lambda_max,num_lambda)
np.savetxt('lambdas.dat',lambdas,fmt='%.8f')
print('\n')

#automation
subprocess.call('./Ex10-Maniscalco.out')

# data collection 
data = np.genfromtxt('for_graphics_eigens.dat')
niter = int(np.genfromtxt('niter.dat'))
dim = 2**(niter+1)

#plot
fig, ax = plt.subplots()
plt.plot(data[:,0],data[:,1],'.',label='RG groundstate')
plt.xlabel('$\lambda$')
plt.ylabel('Energy')
plt.grid()
plt.legend()
ax.text(0.45,0.05,'N = 2^'+str(niter),transform=ax.transAxes,fontsize=14,fontweight='bold')
plt.savefig('niter='+str(niter))
