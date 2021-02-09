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
filename = 'point_42_nopbc_2.xyz'
dt = 1.0

T = 298.15
viscosity = 0.01005
n_particles = 42
n_steps = 1000000#100000
n_write = 1
n_chol = n_steps

from pyBrown.sphere import overlap_pbc, Sphere

initial_coords = [ box_length * np.random.rand(3) - np.array([box_length/2, box_length/2, box_length/2]) ]
for i in range(n_particles):
	print(i+1)
	while(True):
		coord = box_length * np.random.rand(3) - np.array([box_length/2, box_length/2, box_length/2])
		# print(coord)
		overlaps = False
		for j in range(len(initial_coords)):
			if overlap_pbc( Sphere(coord, a), Sphere(initial_coords[j], a), 0.0, box_length ): overlaps = True
		if not overlaps: break
	initial_coords.append(coord)

# print(initial_coords)
# 1/0

# print(initial_coords)
bs = [ Bead(initial_coords[i], a, label) for i in range(n_particles) ]
box = Box(bs, box_length, T, viscosity)
# print( box.beads )
print('a = {} A'.format(box.beads[0].a))
print('a = {} nm'.format(box.beads[0].a/10))
print('D = {} A**2/ps'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity ) ))
print('D = {} A**2/ns'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity )*1000 ))

# A = np.random.rand(n_dim,n_dim)
# A = A @ A.transpose()

box.box_length *= 100

with open(filename, 'w') as output_file:
	start = time.time()
	# for i in tqdm( range(n_steps) ):
	for i in range(n_steps):
		if i % n_write == 0:
			output_file.write('{}\n'.format(len(box.beads)))
			output_file.write('{} time [ps] {}\n'.format(filename, i*dt))
			for bead in box.beads:
				output_file.write('{} {} {} {}\n'.format(bead.label, *bead.r))
		box.propagate(dt, i%n_chol == 0, i%n_chol == 0)

end = time.time()

print('{} seconds elapsed'.format(end-start))	