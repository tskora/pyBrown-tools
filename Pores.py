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

from pyBrown.input_Pores import InputDataPores
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads
from pyBrown.sphere import Sphere, overlap

#-------------------------------------------------------------------------------

def what_radius_based_on_label(label, labels, radii):

	i = 0

	for l, r in zip( labels, radii ):

		if label == l: return radii[i][0]

		i += 1

#-------------------------------------------------------------------------------

def populate_box(input_xyz_filename, labels, radii, snapshot_time):

	populated_box = []
	start_snapshot = False
	with open(input_xyz_filename, 'r') as xyz_file:
		for line in xyz_file:
			if start_snapshot:
				if 'time' in line.split():
					start_snapshot = False
				else:
					if len( line.split() ) == 4:
						label = line.split()[0]
						radius = what_radius_based_on_label(label, labels, radii)
						s = Sphere( [ line.split()[i] for i in range(1,4) ], radius )
						populated_box.append( s )
			if 'time' in line.split():
				if float(line.split()[-1]) == snapshot_time:
					start_snapshot = True
	return populated_box

#-------------------------------------------------------------------------------

def create_grid(grid_density, box_size):

	dx = box_size / grid_density

	grid = [ [x, y, z] for x in np.linspace( -box_size/2, box_size/2, grid_density )
					   for y in np.linspace( -box_size/2, box_size/2, grid_density )
					   for z in np.linspace( -box_size/2, box_size/2, grid_density ) ]

	return grid

#-------------------------------------------------------------------------------

def digitize_grid(input_data):

	import multiprocessing
	from multiprocessing import Pool
	from functools import partial

	input_xyz_filename = input_data["input_xyz_filename"]
	input_labels = input_data["labels"]
	radii = input_data["hydrodynamic_radii"]
	snapshot_time = input_data["snapshot_time"]
	grid_density = input_data["grid_density"]
	box_size = input_data["box_size"]
	output_mode = input_data["output_mode"]

	populated_box = populate_box(input_xyz_filename, input_labels, radii, snapshot_time)
	grid = create_grid(grid_density, box_size)

	digitized_grid_ones = []
	digitized_grid_zeros = []

	nproc = multiprocessing.cpu_count()
	print('You have {0:1d} CPUs'.format(nproc))

	points_per_proc = len(grid) // nproc
	grid_for_proc = [ grid[m*points_per_proc:(m+1)*points_per_proc] for m in range(nproc - 1) ] \
					+ [ grid[(nproc-1)*points_per_proc:] ]

	pool = Pool(processes=nproc)

	_digitize_grid_partial = partial( _digitize_grid, populated_box = populated_box, dx = box_size / grid_density, box_size = box_size, output_mode = output_mode )
	digitized_grid_results = pool.map( _digitize_grid_partial, grid_for_proc )	

	for element in digitized_grid_results:
		digitized_grid_ones += element[0]
		digitized_grid_zeros += element[1]

	return digitized_grid_ones, digitized_grid_zeros

#-------------------------------------------------------------------------------

def _digitize_grid(grid, populated_box, dx, box_size, output_mode):

	if output_mode == 'both':
		ones = True
		zeros = True
	elif output_mode == 'ones':
		ones = True
		zeros = False
	elif output_mode == 'zeros':
		ones = False
		zeros = True
	else:
		print('Error')
		return None

	ones_grid = []
	zeros_grid = []

	all_points = 0
	overlapping_points = 0

	for point in grid:
		all_points += 1
		radius = 0.5 * dx
		with Sphere(point, radius) as probe:
			if overlap( probe, populated_box, 0.0 ):
				overlapping_points += 1
				print('{} {} {} 1'.format(*point))
				ones_grid.append(point)
			else:
				versors = [ np.array([nx * box_size,
									  ny * box_size,
									  nz * box_size])
							for nx in np.arange(-1, 2, 1)
							for ny in np.arange(-1, 2, 1)
							for nz in np.arange(-1, 2, 1) ]
				if_overlap = False
				for versor in versors:
					probe.translate( versor )
					if overlap(probe, populated_box, 0.0):
						if_overlap = True
					probe.translate( -versor )
				if if_overlap:
					overlapping_points += 1
					print('{} {} {} 1'.format(*point))
					ones_grid.append(point)
				else:
					print('{} {} {} 0'.format(*point))
					zeros_grid.append(point)

	return ones_grid, zeros_grid

#-------------------------------------------------------------------------------

def write_digitized_grid_to_file(input_data, digitized_grid_ones, digitized_grid_zeros):

	input_xyz_filename = input_data["input_xyz_filename"]

	with open(input_xyz_filename[:-4]+'_pores.txt', 'w') as pores_file:
		for point in digitized_grid_ones:
			pores_file.write( '{} {} {} 1\n'.format(*point) )
		for point in digitized_grid_zeros:
			pores_file.write( '{} {} {} 0\n'.format(*point) )

	if input_data["verbose"]: print( len(digitized_grid_ones) / ( len(digitized_grid_ones) + len(digitized_grid_zeros) ) )

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "box_size", "input_xyz_filename", "snapshot_time",
						 "grid_density", "hydrodynamic_radii"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "float_type": 32,
				"output_mode": "both" }

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataPores(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Digitizing grid' )
	digitized_grid_ones, digitized_grid_zeros = digitize_grid(i.input_data)

	timestamp( 'Writing digitized grid to file' )
	write_digitized_grid_to_file( i.input_data, digitized_grid_ones, digitized_grid_zeros )
	
#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
