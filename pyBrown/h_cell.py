# pyBrown is a bundle of tools useful for Brownian and Stokesian dynamics simulations
# Copyright (C) 2018  Tomasz Skora (tskora@ichf.edu.pl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses.

import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

class H_cell:

	def __init__(self, a, b, diff_coef, grid_points = 100, cs = None):

		self.a = a
		self.b = b
		self.diff_coef = diff_coef
		self.grid_points = grid_points

		self.dy = 2 * self.a / ( self.grid_points - 1 )

		if cs is None:
		
			self.cs = np.array( [1.0] + [1 for i in range( int(grid_points / 2) - 1 )] +
						    	[0 for i in range( grid_points - int(grid_points / 2) -1 )] +
						    	[0.0])

		else:

			self.cs = np.array( cs )

		self.ys = np.arange(-self.a, self.a + self.dy / 2, self.dy)

	#---------------------------------------------------------------------------

	def __str__(self):

		return '{}x{}; D={}; n={}'.format(self.a,
										  self.b,
										  self.diff_coef,
										  self.grid_points)

	#---------------------------------------------------------------------------

	def __repr__(self):

		return self.__str__()

	#---------------------------------------------------------------------------

	def flow_velocity(self, y):

		m = 1.7 + 0.5 * ( self.b / self.a )**( -1.4 )

		return ( m + 1 ) / m * ( 1 - ( np.abs(y) / self.a )**m )

	#---------------------------------------------------------------------------

	def D_eff(self, y):

		return self.diff_coef / self.flow_velocity(y)

	#---------------------------------------------------------------------------

	def c(self, i):

		if i in range(self.grid_points):

			return self.cs[i]

		else:

			return 0.0

	#---------------------------------------------------------------------------

	def alpha(self, i, dx):

		if i in range(1, self.grid_points - 1):

			return self.D_eff( self.ys[i] ) * dx / 2 / (self.dy)**2

		else:

			return 0.0

	#---------------------------------------------------------------------------

	def propagate(self, dx):

		b = np.array( [ ( ( 1 - 2 * self.alpha(i, dx) ) * self.c(i) +
						  self.alpha(i, dx) * self.c(i - 1) +
						  self.alpha(i, dx) * self.c(i + 1) )
						for i in range(self.grid_points) ] )

		A = np.zeros( (self.grid_points, self.grid_points), float )

		for i in range(self.grid_points):

			A[i, i] = 1 + 2 * self.alpha(i, dx)

			if i > 0: A[i, i - 1] = -self.alpha(i, dx)

			if i < self.grid_points - 1: A[i, i + 1] = -self.alpha(i, dx)

		self.cs = np.linalg.solve(A, b)

	#---------------------------------------------------------------------------

	def save_cell_to_file(self, filename):

		with open(filename + '.txt', 'w') as output:

			output.write( '{} {} {} {}\n'.format(self.a, self.b,
												 self.diff_coef,
												 self.grid_points) )

			for i in range(self.grid_points):

				output.write( '{} {}'.format( self.ys[i], self.cs[i] ) )

	#---------------------------------------------------------------------------

	def read_cell_from_file(filename):

		with open(filename, 'r') as input:

			a = float( input.readline().split()[0] )
			b = float( input.readline().split()[1] )
			diff_coef = float( input.readline().split()[2] )
			grid_points = int( input.readline().split()[0] )

			for line in input:

				print('')

		return self.__init__(a, b, diff_coef, grid_points)






# ---

h = H_cell(1, 1, 1)

plt.plot(h.ys, h.cs, '--')

for j in range(100):

	for i in range(1000):

		h.propagate(0.0001)

	plt.plot(h.ys, h.cs, '-', label = str(j))

plt.legend()

plt.savefig('hcell.jpg', dpi = 300)

# ---

# vels = [ h.flow_velocity( i ) for i in np.linspace(-1, 1, 21) ]

# print( vels )

# print( np.mean(vels) )
