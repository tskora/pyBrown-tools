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
import matplotlib.patches as mpatches
from scipy.linalg import solve_banded

from tqdm import tqdm

import plot_config

#-------------------------------------------------------------------------------

class H_cell:

	def __init__(self, a, b, diff_coef, grid_points, cs = None):

		self.a = a
		self.b = b
		self.diff_coef = diff_coef
		self.grid_points = grid_points

		self.dy = 2 * self.a / ( self.grid_points + 1 )

		if cs is None:
		
			self.cs = np.array( [1.0 for i in range( int(grid_points / 2) )] +
						    	[0.0 for i in range( grid_points - int(grid_points / 2) )] )

		else:

			self.cs = np.array( cs )

		self.ys = np.linspace(-self.a + self.dy, self.a - self.dy, grid_points)

	#---------------------------------------------------------------------------

	def __str__(self):

		return '{}x{}; D={}; n={};\n{}'.format(self.a,
										  self.b,
										  self.diff_coef,
										  self.grid_points,
										  self.cs)

	#---------------------------------------------------------------------------

	def __repr__(self):

		return self.__str__()

	#---------------------------------------------------------------------------

	def flow_velocity(self, y):

		m = 1.7 + 0.5 * ( self.b / self.a )**( -1.4 )

		return 1.0#( m + 1 ) / m * ( 1 - ( np.abs(y) / self.a )**m )

	#---------------------------------------------------------------------------

	def D_eff(self, y):

		return self.diff_coef / self.flow_velocity(y)

	#---------------------------------------------------------------------------

	def c(self, i):

		if i in range(self.grid_points):

			return self.cs[i]

		elif i < 0:

			return self.cs[0]

		else:

			return self.cs[-1]

	#---------------------------------------------------------------------------

	def alpha(self, i, dx):

		if i in range(self.grid_points):

			return self.D_eff( self.ys[i] ) * dx / 2 / (self.dy)**2

		else:

			return 0.0

	#---------------------------------------------------------------------------

	def propagate_euler(self, dx):

		for i in range(self.grid_points):

			incr =  2 * self.alpha(i, dx) * \
					( self.c(i + 1) - 2 * self.c(i) + self.c(i - 1) )

			self.cs[i] += incr

	#---------------------------------------------------------------------------

	def propagate_cn(self, dx):

		# --- CRANK-NICOLSON VERSION ---

		# b = np.array( [ ( ( 1 - 2 * self.alpha(i, dx) ) * self.c(i) +
		# 				  self.alpha(i, dx) * self.c(i - 1) +
		# 				  self.alpha(i, dx) * self.c(i + 1) )
		# 				for i in range(self.grid_points) ] )

		# A = np.zeros( (self.grid_points, self.grid_points), float )

		# for i in range(self.grid_points):

		# 	A[i, i] = 1 + 2 * self.alpha(i, dx)

		# 	if i > 0: A[i, i - 1] = -self.alpha(i, dx)

		# 	if i < self.grid_points - 1: A[i, i + 1] = -self.alpha(i, dx)

		# ab0 = np.array( [0.0] + [ A[i, i + 1] for i in range(self.grid_points - 1) ] ) 
		# ab1 = np.array( [ A[i, i] for i in range(self.grid_points) ] )
		# ab2 = np.array( [ A[i, i - 1] for i in range(1, self.grid_points) ] + [0.0] )
		# ab = np.array( [ab0, ab1, ab2] )

		# self.cs = solve_banded((1, 1), ab, b)

		# --- CRANK-NICOLSON WITH BOUNDARIES

		# assert len(self.cs) == self.grid_points, 'Concentration array has incorrect size'

		# b = np.array( [1.0] + [ ( ( 1 - 2 * self.alpha(i, dx) ) * self.c(i) +
		# 				  self.alpha(i, dx) * self.c(i - 1) +
		# 				  self.alpha(i, dx) * self.c(i + 1) )
		# 				for i in range(self.grid_points) ] + [0.0] )

		# A = np.zeros( (self.grid_points + 2, self.grid_points + 2), float )

		# A[0, 0] = 1.0
		# A[self.grid_points + 1, self.grid_points + 1] = 1.0

		# for i in range(self.grid_points):

		# 	A[i + 1, i + 1] = 1 + 2 * self.alpha(i, dx)

		# 	A[i + 1, i] = -self.alpha(i, dx)

		# 	A[i + 1, i + 2] = -self.alpha(i, dx)

		# ab0 = np.array( [0.0] + [ A[i, i + 1] for i in range(self.grid_points + 1) ] ) 
		# ab1 = np.array( [ A[i, i] for i in range(self.grid_points + 2) ] )
		# ab2 = np.array( [ A[i, i - 1] for i in range(1, self.grid_points + 2) ] + [0.0] )
		# ab = np.array( [ab0, ab1, ab2] )

		# self.cs = solve_banded((1, 1), ab, b)[1:-1]

		return 0

	#---------------------------------------------------------------------------

	def save_cell_to_file(self, filename):

		with open(filename, 'w') as output:

			output.write( '{} {} {} {}\n'.format(self.a, self.b,
												 self.diff_coef,
												 self.grid_points) )

			for i in range(self.grid_points):

				output.write( '{} {}\n'.format( self.ys[i], self.cs[i] ) )

	#---------------------------------------------------------------------------

	@classmethod
	def read_cell_from_file(cls, filename):

		with open(filename, 'r') as input:

			first_line = input.readline()

			a = float( first_line.split()[0] )
			b = float( first_line.split()[1] )
			diff_coef = float( first_line.split()[2] )
			grid_points = int( first_line.split()[3] )

			cs = [ ]

			for line in input:

				cs.append( float( line.split()[1] ) )

		return cls(a, b, diff_coef, grid_points, cs)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	h = H_cell.read_cell_from_file('h_cell.txt')
	print(h)

	# a = 1.0
	# b = 1.0
	# D = 1.0
	# step = 0.001
	# n_snapshots = 5
	# x_max = 1
	# grid_points = 40

	# h = H_cell(a, b, D, grid_points)
	# h1 = H_cell(a, b, D, grid_points)

	# # a = 1.0
	# # b = 1.0
	# # D = 1.0
	# # D_enh = 1.3
	# # step = 0.000001
	# # n_snapshots = 5
	# # x_max = 0.5
	# # grid_points = 50

	# # h = H_cell(a, b, D, grid_points)
	# # h1 = H_cell(a, b, D, grid_points)
	# # h_up = H_cell(a, b, D_enh, grid_points)

	# handles, labels = plt.gca().get_legend_handles_labels()
	# h_patch = mpatches.Patch(color='red', label=r'cn')
	# h_up_patch = mpatches.Patch(color='blue', label=r'euler')
	# handles.append(h_patch)
	# handles.append(h_up_patch)
	# labels.append(r'$D = 1$')
	# labels.append(r'$D = 1.3$')

	# plt.legend(handles, labels)
	
	# plt.plot(h.ys, h.cs, '--', color='black')

	# for j in tqdm( range(n_snapshots) ):
	
	# 	for i in tqdm( range( int( x_max / n_snapshots / step ) ) ):
	
	# 		h.propagate(step)

	# 		h1.propagate_euler(step)

	# 		# h_up.propagate(step)
	
	# 	plt.plot( h.ys, h.cs, '--', color = 'red' )

	# 	plt.plot( h1.ys, h1.cs, '-', color = 'blue'  )

	# 	# plt.plot( h_up.ys, h_up.cs, '-', label = str( ( j * 10000 + i ) * step ), color = 'blue' )
	
	# plt.savefig('hcell.jpg', dpi = 300)

	# plt.close()

	# h.save_cell_to_file('h_cell.txt')
