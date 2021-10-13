import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fitting_func(x,cc):
        return(2*x**2/cc)
xmax = np.array([3,5,9,10,18,20])
N = np.array([15,40,127,156,500,619])
d = 2*xmax/N
const = 2*xmax**2/N
popt, pcov = curve_fit(fitting_func, xmax,N)

plt.plot(xmax,N,label='Computed values')
plt.plot(xmax,(2/popt[0])*xmax**2,label='Fit')
plt.xlabel('xmax')
plt.ylabel('N')
plt.legend(loc=2)
plt.savefig('Trial_fitting')

print(const)
print(d)
print('c = ',popt[0], 'error: ',pcov[0])