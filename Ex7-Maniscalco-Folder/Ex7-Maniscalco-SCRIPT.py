import numpy as np
import os
import subprocess
import sys
import matplotlib.pyplot as plt

choice = input('Do you want to re-execute the FORTRAN program? (y/n)  ')
if choice == 'y': 
	run=True
elif choice == 'n': 
	run = False
else: 
	print('Wrong choice, breaking.')
	sys.exit()

if run:
	os.system('gfortran -I/usr/local/include  Ex7-Maniscalco-CODE.f90 -llapack -lfftw3 -o Ex7-Maniscalco.out')
	subprocess.call('./Ex7-Maniscalco.out')

evolution = np.genfromtxt('eigenstates.txt')
number = 4
last = np.shape(evolution)[1]
divisor = int(np.floor((last-1)/number))

plt.figure(0)
for ii in range(0,number):
    plt.plot(evolution[:,0],evolution[:,divisor*ii+1],color='b',linewidth=(1/number)*ii+0.1,label='$\psi$'+str(divisor*ii))
    plt.xlabel('n')
    plt.ylabel('f(n)')
    plt.title('Time evolution')
plt.plot(evolution[:,0],evolution[:,last-1],color='b',linewidth=1.1,label='$\psi$'+str(last-2))
plt.legend(loc=2)
plt.savefig('time_evolution.png')

plt.figure(1)
plt.plot(evolution[:,0],evolution[:,1],color='b',linewidth=(2/number),label='$\psi$ 0')
plt.plot(evolution[:,0],evolution[:,last-1],color='b',linewidth=1.1,label='$\psi$'+str(last-2))
plt.xlabel('n')
plt.ylabel('f(n)')
plt.title('First and last function')
plt.legend(loc=2)
plt.axvline(x=1,color='r',linestyle=':')
plt.xticks(np.arange(min(evolution[:,0]), max(evolution[:,0])+1, 1.0))
plt.savefig('first_and_last.png')

