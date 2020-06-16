import numpy as np
from scipy.constants import Boltzmann
import time
from tqdm import tqdm
import math

from pyBrown.box import Box
from pyBrown.bead import Bead

a = 51.0
box_length = 750.0
label = "TST"
filename = 'spiral_1.xyz'
dt = 50
T = 293.15
viscosity = 0.01005
n_particles = 1
n_steps = 10000000 # 1000000
n_chol = 1000
D = 10**19 * Boltzmann / (6 * np.pi) * T / ( a * viscosity ) # angstrom per picosecond squared

initial_coords = [ box_length * np.random.rand(3) - np.array([box_length/2, box_length/2, box_length/2]) for i in range(n_particles) ]
# print(initial_coords)
bs = [ Bead(initial_coords[i], a, label) for i in range(n_particles) ]
box = Box(bs, box_length, T, viscosity)
# print( box.beads )
print('a = {} A'.format(box.beads[0].a))
print('a = {} nm'.format(box.beads[0].a/10))
print('D = {} A**2/ps'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity ) ))
print('D = {} A**2/ns'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity )*1000 ))

import matplotlib.pyplot as plt
# with plt.xkcd():
plt.xlim((-3*box_length, 4*box_length))
plt.ylim((-3*box_length, 4*box_length))
for j in range(-3, 4):
	for k in range(-3, 4):
		plt.arrow(-box_length/2+j*box_length, -box_length/2+k*box_length, box_length, 0)
		plt.arrow(-box_length/2+j*box_length, -box_length/2+k*box_length, 0, box_length)
times = np.linspace(0, n_steps*dt, n_steps)
angles = 0.000001*times
angles_prime = np.array( [0.0] + list( 0.000001*dt*np.ones(n_steps-1) + np.arccos( 12*D*np.sqrt(dt*0.000001/times[1:]) ) ) )
xs = np.sqrt( 6*D*times ) * np.cos(angles)
ys = np.sqrt( 6*D*times ) * np.sin(angles)
xs_prime = np.sqrt( 6*D*times ) * np.cos(angles_prime)
ys_prime = np.sqrt( 6*D*times ) * np.sin(angles_prime)
plt.plot( xs[::1000], ys[::1000], '--', color = 'blue' )
plt.plot( xs_prime[::100], ys_prime[::100], 'o', color = 'red' )
print(angles_prime)
plt.show()
plt.plot( times, xs**2+ys**2 )
print( np.polyfit( times, xs**2+ys**2, 1 )[0]/6 )
print( D )
plt.show()

# A = np.random.rand(n_dim,n_dim)
# A = A @ A.transpose()

with open(filename, 'w') as output_file:
	for i in tqdm( range(n_steps) ):
		for bead in box.beads:
			while xs[i] < -box_length / 2:
				xs[i] += box_length
			while xs[i] >= box_length / 2:
				xs[i] -= box_length
			while ys[i] < -box_length / 2:
				ys[i] += box_length
			while ys[i] >= box_length / 2:
				ys[i] -= box_length
		if i%n_chol==0:
			output_file.write('{}\n'.format(len(box.beads)))
			output_file.write('{} time [ps] {}\n'.format(filename, i*dt))
			output_file.write('{} {} {} 0.0\n'.format(bead.label, xs[i], ys[i]))
			# while xs_prime[i] < -box_length / 2:
			# 	xs_prime[i] += box_length
			# while xs_prime[i] >= box_length / 2:
			# 	xs_prime[i] -= box_length
			# while ys_prime[i] < -box_length / 2:
			# 	ys_prime[i] += box_length
			# while ys_prime[i] >= box_length / 2:
			# 	ys_prime[i] -= box_length
			# output_file.write('{} {} {} 0.0\n'.format(bead.label, xs_prime[i], ys_prime[i]))



	