import numpy as np
from scipy.constants import Boltzmann

class Bead():

	def __init__(self, coords, radius, label):

		self.r = np.array(coords)
		self.a = radius
		self.label = label

	def set_D(self, D):

		self.D = D

	def propagate(self, dt):

		self.r += np.sqrt( 2 * self.D * dt ) * np.random.normal(0.0, 1.0, 3)

	def keep_in_box(self, box_length):

		for i in range(3):
			while self.r[i] < -box_length / 2:
				self.r[i] += box_length
			while self.r[i] >= box_length / 2:
				self.r[i] -= box_length

	def __str__(self):

		return "{}, radius = {}".format(self.r, self.a)

	def __repr__(self):

		return self.__str__()

class Box():

	def __init__(self, beads, box, T, viscosity):

		self.beads = beads
		self.box_length = box_length
		self.T = T
		self.viscosity = viscosity
		for bead in self.beads:
			bead.set_D(10**17 * Boltzmann / (6 * np.pi) * self.T / ( bead.a * self.viscosity ) * 100)

	def propagate(self, dt):

		for bead in self.beads:
			bead.propagate(dt)
			bead.keep_in_box(self.box_length)

	def __str__(self):

		return 'a'

a = 10.0
box_length = 2500000.0
label = "TST"
filename = 'test_1.xyz'
initial_coord = [ 0.0, 0.0, 0.0 ]
dt = 0.5
# D = 10**17 * Boltzmann / (6 * np.pi) * 293.15 / ( a * viscosity ) * 100 # angstrom per picosecond squared
# print('D = {} A**2/ps'.format(D))
# print('D = {} A**2/ns'.format(D*1000))
# print('a = {} A'.format(a))
# print('a = {} nm'.format(a/10))
T = 293.15
viscosity = 0.01005

bs = [ Bead(initial_coord, a, label) for i in range(10) ]
box = Box(bs, box_length, T, viscosity)
print( box.beads )
print('D = {} A**2/ps'.format(box.beads[0].D))
print('D = {} A**2/ns'.format(box.beads[0].D*1000))
print('a = {} A'.format(box.beads[0].a))
print('a = {} nm'.format(box.beads[0].a/10))

with open(filename, 'w') as output_file:

	for i in range(40000000):
		if i % 10000 == 0:
			output_file.write('{}\n'.format(len(box.beads)))
			output_file.write('{} time [ps] {}\n'.format(filename, i*dt))
			for bead in box.beads:
				output_file.write('{} {} {} {}\n'.format(bead.label, *bead.r))
		box.propagate(dt)




	