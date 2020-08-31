# pyBrown is a bundle of tools useful for Brownian and Stokesian dynamics
# simulations. Copyright (C) 2018  Tomasz Skora (tskora@ichf.edu.pl)
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
# along with this program. If not, see https://www.gnu.org/licenses.

import click

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pyBrown.messaging import timestamp

#-------------------------------------------------------------------------------

def read_digitized_grid_from_file(input_voxels_filename):

	with open(input_voxels_filename) as voxels_file:
		first_line = voxels_file.readline().split()
		grid_density = int( first_line[0] )
		box_size = float( first_line[1] )
		digitized_grid = np.zeros((grid_density, grid_density, grid_density))
		for line in voxels_file:
			i, j, k, val = line.split()
			digitized_grid[int(i)][int(j)][int(k)] = int(val)

	return digitized_grid

#-------------------------------------------------------------------------------

def plot_digitized_grid(digitized_grid):

	voxels = np.array( digitized_grid )

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	ax.voxels(voxels, facecolors='blue', edgecolor='k')

	plt.show()
	plt.close()

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
@click.option('-s', '--slice-index',
			  type = int, required = False)
def main(input_filename, slice_index):

	digitized_grid = read_digitized_grid_from_file( input_filename )

	plot_digitized_grid(digitized_grid)

	if slice_index != None:

		plane = digitized_grid[slice_index]

		plt.imshow(plane)
		
		plt.show()
	
#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
