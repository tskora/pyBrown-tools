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
	'--output-filename', '-o',
	help = 'output filename'
	)
def main(a_len, b_len, diff_coefs, grid_points, step, x_max, output_filename):

	Ds = np.array( [ float( element ) for element in diff_coefs.split() ] )
	a = float(a_len)
	b = float(b_len)
	n_grid = int(grid_points)
	dx = float(step)
	x = float(x_max)

	hs = []

	for D in Ds:

		hs.append( H_cell(a, b, D, n_grid) )

	h_cell_simulation( hs, dx, x, output_filename )

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
