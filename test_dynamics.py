import numpy as np
from scipy.constants import Boltzmann
import time
from tqdm import tqdm

class Bead():

	def __init__(self, coords, radius, label):

		self.r = np.array(coords)
		self.a = radius
		self.label = label

	def translate(self, v):

		self.r += v

	def propagate(self, dt, deterministic = False):

		if not deterministic:
			self.r += np.sqrt( 2 * self.D * dt ) * np.random.normal(0.0, 1.0, 3)
		else:
			self.r += self.D * dt * np.array([0.0, 0.0, 1.0])

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

	def __eq__(self, p):

		if isinstance( p, Bead ):
			return ( np.all( self.r == p.r ) ) and ( self.a == p.a )
		return False

class Box():

	def __init__(self, beads, box, T, viscosity, deterministic = False):

		self.beads = beads
		self.box_length = box_length
		self.T = T
		self.viscosity = viscosity
		self.deterministic = deterministic

	def propagate(self, dt, cholesky = True):

		if cholesky: self.D = self.compute_Dmatrix()
		# print(D)
		# print(np.linalg.eig(D))
		if cholesky: self.B = np.linalg.cholesky(self.D)
		# print(self.B)
		BX = self.B @ np.random.normal(0.0, 1.0, 3 * len(self.beads)) * np.sqrt(2 * dt)
		# print(BX)

		for i, bead in enumerate( self.beads ):
			bead.translate( BX[3 * i: 3 * (i + 1)] )
			bead.keep_in_box(self.box_length)

	def compute_Dmatrix(self):

		return Boltzmann * self.T * 10**19 * self.M_rpy()

	def pointer(self, p1, p2, nearest = True):

		dist0 = np.linalg.norm(p2.r - p1.r)

		versors = [ self.box_length * np.array( [i, j, k] ) for i in range(-1, 2, 1)
							  								for j in range(-1, 2, 1)
							  								for k in range(-1, 2, 1) ]

		if nearest and dist0 > self.box_length / 2:

			for versor in versors:

				p1.translate(versor)

				dist1 = np.linalg.norm(p2.r - p1.r)

				if dist1 <= dist0:

					return p2.r - p1.r

				p1.translate(-versor)

		return p2.r - p1.r

	def Mii_rpy(self, p, viscosity):

		return np.identity(3) / ( 6 * np.pi * p.a * viscosity )

	def Mij_rpy(self, pi, pj, viscosity):

		Rh_larger = max( pi.a, pj.a )
		Rh_smaller = min( pi.a, pj.a )

		point = self.pointer(pi, pj)
		dist = np.linalg.norm( point )
		outer = np.outer(point/dist, point/dist)
	
		if dist > ( pi.a + pj.a ):
	
			coef_1 = 1.0 / ( 8 * np.pi * dist * viscosity )
			coef_2 = 1.0 + ( pi.a**2 + pj.a**2) / ( 3 * dist**2 )
			coef_3 = 1.0 - ( pi.a**2 + pj.a**2) / dist**2
	
			answer = coef_2 * np.identity(3)
			answer += coef_3 * outer
			answer *= coef_1
	
			return answer
        
		elif dist <= ( Rh_larger - Rh_smaller ):

			return np.identity(3) / ( 6 * np.pi * Rh_larger * viscosity )

		else:

			coef_1 = 1.0 / ( 6 * np.pi * pi.a * pj.a * viscosity )
			coef_2 = 16 * dist**3 * ( pi.a + pj.a )
			coef_3 = ( (pi.a - pj.a )**2 + 3 * dist**2 )**2
			coef_4 = ( coef_2 - coef_3 ) / ( 32 * dist**3 )
			coef_5 = 3 * ( (pi.a - pj.a)**2 - dist**2 )**2
			coef_6 = coef_5 / ( 32 * dist**3 )

			answer = coef_4 * np.identity(3)
			answer += coef_6 * outer
			answer *= coef_1
        
			return answer

	def M_rpy(self):

		M = [ [ None for j in range( len(self.beads) ) ] for i in range( len(self.beads) ) ]

		for i, pi in enumerate(self.beads):

			M[i][i] = self.Mii_rpy(pi, self.viscosity)

			for j in range(i):

				pj = self.beads[j]

				M[i][j] = self.Mij_rpy(pi, pj, self.viscosity)
				M[j][i] = np.transpose( M[i][j] )

		return np.block(M)

	def __str__(self):

		return 'a'

a = 51.0
box_length = 7500.0
label = "TST"
filename = 'test_15.xyz'
dt = 0.5
# D = 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( a * viscosity ) # angstrom per picosecond squared
# print('D = {} A**2/ps'.format(D))
# print('D = {} A**2/ns'.format(D*1000))
# print('a = {} A'.format(a))
# print('a = {} nm'.format(a/10))
T = 293.15
viscosity = 0.01005
n_particles = 1
n_steps = 1000000
# n_dim = 100
n_chol = 100

initial_coords = [ box_length * np.random.rand(3) - np.array([box_length/2, box_length/2, box_length/2]) for i in range(n_particles) ]
print(initial_coords)
bs = [ Bead(initial_coords[i], a, label) for i in range(n_particles) ]
box = Box(bs, box_length, T, viscosity, True)
print( box.beads )
print('a = {} A'.format(box.beads[0].a))
print('a = {} nm'.format(box.beads[0].a/10))
print('D = {} A**2/ps'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity ) ))
print('D = {} A**2/ns'.format( 10**19 * Boltzmann / (6 * np.pi) * 293.15 / ( box.beads[0].a * viscosity )*1000 ))


# A = np.random.rand(n_dim,n_dim)
# A = A @ A.transpose()

with open(filename, 'w') as output_file:
	start = time.time()
	for i in tqdm( range(n_steps) ):
		# np.linalg.cholesky(A)
		# np.linalg.inv(A)
		# np.linalg.inv(A)
		if i % 10000 == 0:
			output_file.write('{}\n'.format(len(box.beads)))
			output_file.write('{} time [ps] {}\n'.format(filename, i*dt))
			for bead in box.beads:
				output_file.write('{} {} {} {}\n'.format(bead.label, *bead.r))
		box.propagate(dt, i%n_chol == 0)

end = time.time()

print('{} seconds elapsed'.format(end-start))




	