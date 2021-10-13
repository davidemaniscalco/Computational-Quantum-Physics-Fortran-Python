import numpy as np
import os
import subprocess
import sys
import matplotlib.pyplot as plt
from scipy import optimize

separable = np.genfromtxt('separable_times.txt')
general = np.genfromtxt('general_times.txt')


def g(x,beta):
    return beta*2*x
x = np.linspace(0,np.max(separable[:,1]),1000)
param2 = optimize.curve_fit(g, separable[:,1],separable[:,3])
print('beta: ',param2[0])
                
plt.figure(0)
plt.plot(separable[:,1],separable[:,3],color='b',label='separable state time')
plt.plot(x,g(x,param2[0]),label='linear fitting',color='r')
plt.xlabel('N')
plt.ylabel('time(s)')
plt.title('Separable state time, D fixed')
plt.legend(loc='best')
plt.savefig('sep_dfixed.png')
    
def f(x,alfa):
             return alfa*2**x
x = np.linspace(0,np.max(general[:,1]),1000)
param = optimize.curve_fit(f, general[:,1],general[:,3])
print('alfa: ',param[0])
                
plt.figure(1)
plt.plot(general[:,1],general[:,3],color='b',label='general state time')
plt.plot(x,f(x,param[0]),label='exponential fitting',color='r')
plt.xlabel('n')
plt.ylabel('time(s)')
plt.title('General state time, D fixed')
plt.legend(loc='best')
plt.savefig('gen_dfixed.png')
