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

def digitize_grid(grid, populated_box, dx, box_size):

	all_points = 0
	overlapping_points = 0

	for point in grid:

		all_points += 1
		radius = 0.5 * dx

		with Sphere(point, radius) as probe:

			if overlap( probe, populated_box, 0.0 ):
				overlapping_points += 1
				print('{} {} {} 1'.format(*point))
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
				else:
					print('{} {} {} 0'.format(*point))

	return overlapping_points, all_points, overlapping_points/all_points

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
	defaults = {"debug": False, "verbose": False, "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataPores(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	xyz_file = open(i.input_data["input_xyz_filename"], 'r')

	populated_box = []

	start_snapshot = False

	for line in xyz_file:
		if start_snapshot:
			if 'time' in line.split():
				start_snapshot = False
			else:
				if len( line.split() ) == 4:
					label = line.split()[0]
					radius = what_radius_based_on_label(label, i.input_data["labels"], i.input_data["hydrodynamic_radii"])
					s = Sphere( [ line.split()[i] for i in range(1,4) ], radius )
					populated_box.append( s )
		if 'time' in line.split():
			if float(line.split()[-1]) == i.input_data["snapshot_time"]:
				start_snapshot = True

	N = i.input_data["grid_density"]
	dx = i.input_data["box_size"] / N

	grid = [ [x, y, z] for x in np.linspace( -i.input_data["box_size"]/2, i.input_data["box_size"]/2, N )
					   for y in np.linspace( -i.input_data["box_size"]/2, i.input_data["box_size"]/2, N )
					   for z in np.linspace( -i.input_data["box_size"]/2, i.input_data["box_size"]/2, N ) ]

	import multiprocessing
	from multiprocessing import Pool
	from functools import partial

	nproc = multiprocessing.cpu_count()
	print('You have {0:1d} CPUs'.format(nproc))

	points_per_proc = len(grid) // nproc
	grid_for_proc = [ grid[m*points_per_proc:(m+1)*points_per_proc] for m in range(nproc - 1) ] \
					+ [ grid[(nproc-1)*points_per_proc:] ]

	pool = Pool(processes=nproc)

	digitize_grid_partial = partial( digitize_grid, populated_box = populated_box, dx = dx, box_size = i.input_data["box_size"] )

	count = pool.map( digitize_grid_partial, grid_for_proc )
	
	print(count)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
