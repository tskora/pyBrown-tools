import numpy as np
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann

def msd(t, D):
	return 6 * D * t

T = 293.15
viscosity = 0.01005

results = np.genfromtxt('test_msd.txt')
ts = [ results[i][0] for i in range(1, len(results)) ]
msds = [ results[i][1] for i in range(1, len(results)) ]

D = curve_fit(msd, ts, msds)[0]
a = 10**17 * Boltzmann / (6 * np.pi) * T / ( D * viscosity ) * 100

print('D = {} A**2/ps'.format(D))
print('D = {} A**2/ns'.format(D*1000))
print('D = {} A**2/ps'.format(a))
print('D = {} A**2/ps'.format(a/10))

