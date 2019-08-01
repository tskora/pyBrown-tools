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

import click
import numpy as np

from pyBrown.h_cell import H_cell, h_cell_simulation

#-------------------------------------------------------------------------------

@click.command()
@click.option(
	'--a-len', '-a',
	default = 1.0,
	show_default = True,
	help = 'half of the cell width in y direction'
	)
@click.option(
	'--b-len', '-b',
	default = 1.0,
	show_default = True,
	help = 'half of the cell width in z direction'
	)
@click.option(
	'--diff-coefs', '-d',
	default = "1.0",
	show_default = True,
	help = 'diffusion coefficients'
	)
@click.option(
	'--flow-velocity', '-v',
	default = 1.0,
	show_default = True,
	help = 'flow velocity')
@click.option(
	'--grid-points', '-n',
	help = 'number of grid points along y direction'
	)
@click.option(
	'--step', '-x',
	help = 'flow propagation step along x direction'
	)
@click.option(
	'--x-max', '-m',
	help = 'cell length in the x direction'
	)
@click.option(
	'--snapshots', '-s',
	default = "",
	show_default = True,
	help = 'x values for concentration profiles to be plotted')
@click.option(
	'--output-filename', '-o',
	help = 'output filename'
	)
@click.option(
	'--profile', '-p',
	default = 'step',
	show_default = True,
	help = 'initial concentration profile (step and well are implemented)')
def main(a_len, b_len, diff_coefs, grid_points, step, x_max,
		 output_filename, snapshots, flow_velocity, profile):

	Ds = parse_diff_coefs(diff_coefs)
	ss = parse_snapshots(snapshots)

	a = float(a_len)
	b = float(b_len)
	n_grid = int(grid_points)
	dx = float(step)
	x = float(x_max)
	v = float(flow_velocity)

	hs = []

	for D in Ds:

		if profile == 'step':
			hs.append( H_cell(a, b, D, n_grid, v = v) )
		elif profile == 'well':
			if n_grid % 3 == 0:
				cs = [1.0 for i in range(n_grid // 3)] +\
				[0.0 for i in range(n_grid // 3)] +\
				[1.0 for i in range(n_grid // 3)]
			elif n_grid %3 == 1:
				cs = [1.0 for i in range(n_grid // 3)] +\
				[0.0 for i in range(n_grid // 3 + 1)] +\
				[1.0 for i in range(n_grid // 3)]
			else:
				cs = [1.0 for i in range(n_grid // 3 + 1)] +\
				[0.0 for i in range(n_grid // 3)] +\
				[1.0 for i in range(n_grid // 3 + 1)]
		elif profile == 'quarter':
			if n_grid % 4 == 0:
				cs = [1.0 for i in range(n_grid//4)] +\
				[0.0 for i in range(3*n_grid//4)]
			elif n_grid % 4 == 1:
				# TODO
				cs = [0.0 for i in range(n_grid)]
			elif n_grid % 4 == 2:
				# TODO
				cs = [0.0 for i in range(n_grid)]
			else:
				# TODO
				cs = [0.0 for i in range(n_grid)]

			hs.append( H_cell(a, b, D, n_grid, cs, v) )

	h_cell_simulation( hs, dx, x, output_filename, ss )

#-------------------------------------------------------------------------------

def parse_diff_coefs(diff_coefs):

	if not ( ":" in diff_coefs ):

		return np.array( [ float( element )
			for element in diff_coefs.split() ] )

	else:

		return np.linspace( *[ float( diff_coefs.split(':')[i] )
			for i in range(3) ] )

#-------------------------------------------------------------------------------

def parse_snapshots(snapshots):

	if not ( ":" in snapshots ):

		return [ float( element ) for element in snapshots.split() ]

	else:

		return list( np.linspace( *[ float( snapshots.split(':')[i] )
			for i in range(3) ] ) )

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
