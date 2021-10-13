import numpy as np
import os
import subprocess
import sys
import matplotlib.pyplot as plt

# lambda generation
lambda_min = float(input('Choose the minimum value of lambda: '))
lambda_max = float(input('Choose the maximum value of lambda: '))
num_lambda = int(input('Choose the number of lambdas: '))
lambdas = np.linspace(lambda_min,lambda_max,num_lambda)
np.savetxt('lambdas.dat',lambdas,fmt='%.8f')

#automation
subprocess.call('./Ex9-Maniscalco.out')

#data collection and plot
data = np.genfromtxt('for_graphics_eigens.dat')
Nandk = np.genfromtxt('N.dat')
NN = int(Nandk[0])
kk = int(Nandk[1])

x = np.linspace(lambda_min,lambda_max,200)
vec_nn = np.linspace(NN,NN+1,200)

def mean_field(lam,NN):
	y = np.where(lam<2.0, -1-lam**2/4, -abs(lam))
	return(NN*y)

fig, ax = plt.subplots()
plt.xlabel('$\lambda$')
plt.ylabel('Energy')
#plt.title('N = '+ str(NN))
plt.grid()
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.text(0.05,0.05,'N = '+ str(NN),transform=ax.transAxes,fontsize=14,
        bbox = dict(boxstyle='square', facecolor='wheat', alpha=0.3))


for ii in range(1,kk+1):
    plt.plot(data[:,0],data[:,ii],'.',label='Eigenvalue'+str(ii),alpha=0.5,ms=5.0)
#plt.plot(x, mean_field(x,NN))
#plt.legend(loc='best')
plt.savefig('N='+str(NN))

#fits
ii = 1
trtr= 2.0
p , residuals, rank, singular_values, rcond = np.polyfit(data[data[:,0]>trtr,0],data[data[:,0]>trtr,ii],1,full=True)
p2 , residuals2, rank2, singular_values2, rcond2 = np.polyfit(data[data[:,0]<trtr,0],data[data[:,0]<trtr,ii],2,full=True)
ig, ax = plt.subplots()
plt.xlabel('$\lambda$')
plt.ylabel('Energy')
plt.grid()
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.text(0.05,0.05,'N = '+ str(NN),transform=ax.transAxes,fontsize=14,
        bbox = dict(boxstyle='square', facecolor='wheat', alpha=0.3))

# warning: leave 'int' if lambda min and lambda max were not integers
plt.plot(data[:,0],data[:,ii],'r.',label='Eigenvalue'+str(ii),alpha=0.9,ms=5.0)
plt.plot(x, p[0]*x+p[1],label='Linear fit in ['+str(int(trtr))+':'+str(int(lambda_max))+']')
plt.plot(x,p2[0]*x**2+p2[1]*x+p2[2],label='Quadratic fit in ['+str(int(lambda_min))+':'+str(int(trtr))+']')
plt.legend(loc='best')
plt.savefig('N='+str(NN)+'_fits')

print('Threshold=',trtr,' sum of residuals linear fit: ',residuals)
print('Threshold=',trtr,' sum of residuals quadratic fit: ',residuals2)
