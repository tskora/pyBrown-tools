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
# from scipy.linalg import solve_banded

from tqdm import tqdm

from pyBrown.plot_config import plot_config

#-------------------------------------------------------------------------------

class H_cell:

	def __init__(self, a, b, diff_coef, grid_points, cs = None, v = 1):

		self.a = a
		self.b = b
		self.diff_coef = diff_coef
		self.grid_points = grid_points
		self.v = v

		self.dy = 2 * self.a / ( self.grid_points + 1 )

		if cs is None:
		
			self.cs = np.array( [0.0 for i in range( int(self.grid_points / 2) )] + 
				[1.0 for i in range( self.grid_points - int(self.grid_points / 2) )],
				float )

		else:

			self.cs = np.array( cs )
			assert len( self.cs ) == self.grid_points, 'Wrong size of concentration array'

		self.ys = np.linspace(-self.a + self.dy, self.a - self.dy, self.grid_points)

	#---------------------------------------------------------------------------

	def __str__(self):

		return '{}x{}; D={}; n={}; v={};\n{}'.format(
												self.a,
										  		self.b,
										  		self.diff_coef,
										  		self.grid_points,
										  		self.v,
										  		self.cs
										  	  )

	#---------------------------------------------------------------------------

	def __repr__(self):

		return self.__str__()

	#---------------------------------------------------------------------------

	def flow_velocity(self, y):

		m = 1.7 + 0.5 * ( self.b / self.a )**( -1.4 )

		return self.v * ( m + 1 ) / m * ( 1 - ( np.abs(y) / self.a )**m )

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

	def _propagate_euler_step(self, dx):

		new_cs = []

		for i in range(self.grid_points):

			incr =  2 * self.alpha(i, dx) * \
					( self.c(i + 1) - 2 * self.c(i) + self.c(i - 1) )

			new_c = self.cs[i] + incr

			assert new_c >= 0.0, 'Concentration becomes negative during propagation'

			new_cs.append( new_c )

		self.cs = np.array( new_cs, float )

	#---------------------------------------------------------------------------

	def _propagate_cn_step(self, dx):

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

		return None

	#---------------------------------------------------------------------------

	def propagate_euler(self, dx, x_max):

		for j in tqdm( range( int( x_max / dx ) ) ):

			self._propagate_euler_step(dx)

	#---------------------------------------------------------------------------

	def save_cell_to_file(self, filename):

		with open(filename, 'w') as output:

			output.write( '{} {} {} {} {}\n'.format(self.a, self.b,
												 self.diff_coef,
												 self.grid_points,
												 self.v) )

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
			v = float( first_line.split()[4] )

			cs = [ ]

			for line in input:

				cs.append( float( line.split()[1] ) )

		return cls(a, b, diff_coef, grid_points, cs, v)

#-------------------------------------------------------------------------------

def init_plot_concentration_profile(a, b, grid_points, dx, v):

	colors, symbols = plot_config()

	plt.grid()

	plt.xlabel(r'$y$')

	plt.ylabel(r'$c$/$c_{input}$')

	plt.title( r'$a = $ {}; $b = $ {}; $n_g = $ {}; $dx = $ {}; $v = $ {}'.format(
				a, b, grid_points, dx, v) )

	return colors, symbols

#-------------------------------------------------------------------------------

def finish_plot_concentration_profile(output):

	plt.legend()

	plt.ylim((0, 1))

	plt.savefig(output + '.jpg', dpi = 300)

	plt.close()

#-------------------------------------------------------------------------------

def h_cell_simulation(h_cells, dx, x_max, output, snapshots = []):

	if not isinstance( h_cells, list ): h_cells = [ h_cells ]

	colors, symbols = init_plot_concentration_profile(h_cells[0].a, h_cells[0].b,
						h_cells[0].grid_points, dx, h_cells[0].v)

	for snapshot in snapshots: assert x_max > snapshot, 'x_max has to be greater than x for any snapshot'
	snapshots.append(x_max)
	snapshots.sort()

	for i, h in enumerate( h_cells ):

		x0 = 0.0

		for j, x_snapshot in enumerate( snapshots ):

			x = x_snapshot - x0
			x0 = x_snapshot

			h.propagate_euler(dx, x)
			h.save_cell_to_file(output + '_' + str(i) + '_snapshot_' + str(x_snapshot) + '.txt')

			plt.plot(h.ys, h.cs,
					 symbols[j % len(symbols)],
					 color = colors[i % len(colors)],
					 label = r'$D = $ {}; $x = $ {}'.format(h.diff_coef, x_snapshot))

	finish_plot_concentration_profile(output)

#-------------------------------------------------------------------------------

def read_cells_from_files(input_template, snapshots, number_of_substances):

	hs = []

	for s in snapshots:

		hs.append(
			[ H_cell.read_cell_from_file(
				'{}_{}_snapshot_{}.txt'.format(input_template, i, s) )
			for i in range(number_of_substances) ]
			)

	Ds = [ h.diff_coef for h in hs[0] ]

	return hs, Ds

#-------------------------------------------------------------------------------

def _integrate_concentration_profile(cs, outlets = 2):

	assert len(cs) % outlets == 0, 'number of grid points should be a multiple of number of outputs'

	segment_length = len(cs) // outlets

	c_int = np.array( [ np.mean( cs[ i * segment_length: (i + 1) * segment_length ] )
		for i in range(outlets) ], float )

	return c_int

#-------------------------------------------------------------------------------

def perform_differential_analysis_for_2_substances(h_cells, diff_coefs, snapshots, number_of_outlets, input_template):

	dcps = integrate_differential_concentration_profiles(h_cells, len(snapshots), number_of_outlets)

	a = h_cells[0][0].a
	b = h_cells[0][0].b
	v = h_cells[0][0].v

	plot_differential_concentration_profiles(snapshots, dcps, diff_coefs, a, b, v, input_template)

#-------------------------------------------------------------------------------

def perform_analysis_for_many_substances(h_cells, diff_coefs, snapshots, number_of_outlets, outlet_index, input_template):

	a = h_cells[0][0].a
	b = h_cells[0][0].b
	v = h_cells[0][0].v

	cs = np.array( [ [ h_cells[i][j].cs for j in range( len(diff_coefs) ) ]
					for i in range( len(snapshots) ) ], float )

	colors, symbols = plot_config()

	plt.grid()

	plt.xlabel(r'$D [\frac{m^2}{s}]$')

	plt.ylabel(r'$c$/$c_{input}$')

	plt.ylim((0.0,0.5))

	plt.title( 'Concentration (in outlet {}/{}) dependence on diffusion coefficient\n'.format(outlet_index, number_of_outlets) +
				r'$a = $ {}; $b = $ {}; $v = $ {}'.format(
				a, b, v) )

	for i in range( len(snapshots) ):

		cs_int_D = []

		for j in range( len(diff_coefs) ):

			cs_int_D.append( _integrate_concentration_profile(cs[i][j], number_of_outlets)[outlet_index - 1] )

		plt.plot( diff_coefs, cs_int_D, '-o', label = r'snapshot $x =$ {}'.format(snapshots[i]), color = colors[i] )

	plt.legend()

	plt.savefig(input_template + '_{}_of_{}_outs.jpg'.format(outlet_index, number_of_outlets), dpi = 300)

	plt.close()

#-------------------------------------------------------------------------------

def integrate_differential_concentration_profiles(hs, num_of_snapshots, num_of_outlets):

	num=of_subs = len(hs)

	dcps = np.array(
		[ _integrate_concentration_profile(hs[i][0].cs - hs[i][1].cs, num_of_outlets)
		for i in range( num_of_snapshots ) ], float )

	return dcps

#-------------------------------------------------------------------------------

def plot_differential_concentration_profiles(ss, dcps, Ds, a, b, v, input_template):

	colors, symbols = plot_config()

	plt.grid()

	plt.xlabel(r'$x [m]$')

	plt.ylabel(r'$\Delta c$/$c_{input}$')

	plt.title( 'Differential concentration profiles\n' +
				r'$a = $ {}; $b = $ {}; $v = $ {}; $D = $ {}'.format(
				a, b, v, Ds) )

	for i in range(len(dcps[0])):
		plt.plot( ss, dcps[:,i], '-o', label = 'outlet nr {}'.format(i + 1),
									   color = colors[i % len(colors)] )
	
	plt.legend()

	plt.savefig(input_template + '_dcp.jpg', dpi = 300)

	plt.close()
