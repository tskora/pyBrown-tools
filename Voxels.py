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

from pyBrown_tools.input_Voxels import InputDataVoxels
from pyBrown_tools.messaging import timestamp
from pyBrown_tools.trajectories import read_trajectories, add_auxiliary_data_multibeads
from pyBrown_tools.sphere import Sphere, overlap

#-------------------------------------------------------------------------------

def read_radii_from_str_file(input_str_filename, mode):

	with open(input_str_filename) as str_file:
		labels = []
		radii = []
		for line in str_file:
			if line.split()[0] == 'sub':
				label = line.split()[1]
				hydrodynamic_radius = float(line.split()[6])
				lennard_jones_radius = float(line.split()[8]) / 2

				if label not in labels:
					labels.append(label)
					if mode == 'hydrodynamic': radii.append( hydrodynamic_radius )
					elif mode == 'lennard-jones': radii.append( lennard_jones_radius )
					else: print('radii mode unknown')

	return labels, radii

#-------------------------------------------------------------------------------

def what_radius_based_on_label(label, labels, radii):

	i = 0

	for l, r in zip( labels, radii ):

		if label == l: return radii[i]

		i += 1

#-------------------------------------------------------------------------------

def populate_box(input_xyz_filename, labels, radii, snapshot_time):

    populated_box = []
    start_snapshot = False
    with open(input_xyz_filename, 'r') as xyz_file:
    	for line in xyz_file:
            if start_snapshot:
                # 'time' with start_snapshot means the next snapshot, so break
                if 'time' in line.split():
                    break;
                elif len( line.split() ) == 4:
                    label = line.split()[0]
                    radius = what_radius_based_on_label(label, labels, radii)
                    s = Sphere( [ line.split()[i] for i in range(1,4) ], radius )
                    populated_box.append( s )

            if 'time' in line.split():
                time = float(line.split()[-1])  
                if time >= snapshot_time:
                    print('Saving snapshot for time %f (requested: %f)' % (time, snapshot_time))
                    start_snapshot = True

    return populated_box

#-------------------------------------------------------------------------------

def create_grid(grid_density, box_size):

	dx = box_size / grid_density

	grid = [ [x, y, z] for x in np.linspace( -(box_size - dx)/2, (box_size - dx)/2, grid_density )
					   for y in np.linspace( -(box_size - dx)/2, (box_size - dx)/2, grid_density )
					   for z in np.linspace( -(box_size - dx)/2, (box_size - dx)/2, grid_density ) ]

	i = 0
	for x in np.linspace( -(box_size - dx)/2, (box_size - dx)/2, grid_density ):
		idx_float = ( x + ( box_size - dx ) / 2 ) / dx
		idx = int( round( idx_float, 0 ) )
		assert i == idx, 'float addition error'
		i += 1

	return grid

#-------------------------------------------------------------------------------

def digitize_grid(input_data):

    import multiprocessing
    from multiprocessing import Pool
    from functools import partial
    
    input_xyz_filename = input_data["input_xyz_filename"]
    input_labels = input_data["labels"]
    radii = input_data["radii"]
    snapshot_time = input_data["snapshot_time"]
    grid_density = input_data["grid_density"]
    box_size = input_data["box_size"]
    nproc = input_data["omp_cores"]

    populated_box = populate_box(input_xyz_filename, input_labels, radii, snapshot_time)
    grid = create_grid(grid_density, box_size)

    digitized_grid = np.zeros( (grid_density, grid_density, grid_density) )
    
    if nproc == 0: 
        nproc = multiprocessing.cpu_count()
    print('You have {0:1d} CPUs'.format(nproc))

    points_per_proc = len(grid) // nproc
    grid_for_proc = [ grid[m*points_per_proc:(m+1)*points_per_proc] for m in range(nproc - 1) ] \
					+ [ grid[(nproc-1)*points_per_proc:] ]

    pool = Pool(processes=nproc)

    _digitize_grid_partial = partial( _digitize_grid, populated_box = populated_box, dx = box_size / grid_density, box_size = box_size )
    overlap_indices_partial = pool.map( _digitize_grid_partial, grid_for_proc )	

    for element in overlap_indices_partial:
        for indices in element:
            digitized_grid[indices[0]][indices[1]][indices[2]] = 1.0

    print( np.sum(digitized_grid)/grid_density**3 )

    return digitized_grid

#-------------------------------------------------------------------------------

def _digitize_grid(grid, populated_box, dx, box_size):

	overlap_indices = []

	all_points = 0
	overlapping_points = 0

	for point in grid:
		all_points += 1
		radius = 0.5 * dx
		with Sphere(point, radius) as probe:
			if overlap( probe, populated_box, 0.0 ):
				overlapping_points += 1
				i = int( round( ( point[0] + ( box_size - dx ) / 2 ) / dx, 0 ) )
				j = int( round( ( point[1] + ( box_size - dx ) / 2 ) / dx, 0 ) )
				k = int( round( ( point[2] + ( box_size - dx ) / 2 ) / dx, 0 ) )
				# print( '{} {} {}'.format( i, j, k ) )
				overlap_indices.append( (i, j, k) )
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
					if if_overlap: break
				if if_overlap:
					overlapping_points += 1
					# print('{} {} {} 1'.format(*point))
					i = int( round( ( point[0] + ( box_size - dx ) / 2 ) / dx, 0 ) )
					j = int( round( ( point[1] + ( box_size - dx ) / 2 ) / dx, 0 ) )
					k = int( round( ( point[2] + ( box_size - dx ) / 2 ) / dx, 0 ) )
					# print( '{} {} {}'.format( i, j, k ) )
					overlap_indices.append( (i, j, k) )

	return overlap_indices

#-------------------------------------------------------------------------------

def write_digitized_grid_to_file(input_data, digitized_grid):
    
    #input_xyz_filename = input_data["input_xyz_filename"]
    output_filename = input_data["output_filename"]

    #with open(input_xyz_filename[:-4]+'_voxels.txt', 'w') as voxels_file:
    with open(output_filename, 'w') as voxels_file:
        voxels_file.write('{} {}\n'.format(input_data["grid_density"], input_data["box_size"]))
        for i in range( input_data["grid_density"] ):
            for j in range( input_data["grid_density"] ):
                for k in range( input_data["grid_density"] ):
                    voxels_file.write('{} {} {} {}\n'.format(i, j, k, int(digitized_grid[i][j][k])))

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

	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D

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
def main(input_filename):

    # here the list of keywords that are required for program to work is provided
    required_keywords = ["box_size", "input_xyz_filename", "input_str_filename",
						 "snapshot_time", "grid_density", "radii_mode"]
    # here the dict of keywords:default values is provided
    # if given keyword is absent in JSON, it is added with respective default value
    defaults = {"debug": False, "verbose": False, "float_type": 32, "omp_cores": 0 }

    timestamp( 'Reading input from {} file', input_filename )
    i = InputDataVoxels(input_filename, required_keywords, defaults)
    i.input_data["labels"], i.input_data["radii"] = read_radii_from_str_file(i.input_data["input_str_filename"], i.input_data["radii_mode"])
	
    timestamp( 'Input data:\n{}', i )

    timestamp( 'Digitizing grid' )
    digitized_grid = digitize_grid(i.input_data)

    timestamp( 'Writing digitized grid to file' )
    write_digitized_grid_to_file( i.input_data, digitized_grid )
	
#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
