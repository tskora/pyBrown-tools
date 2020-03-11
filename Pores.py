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

from pyBrown.input_Pores import InputDataPores
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads
from pyBrown.sphere import Sphere, overlap

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "box_size", "input_xyz_template", "input_xyz_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataPores(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )

	import numpy as np

	N = 100
	dx = i.input_data["box_size"] / N

	grid = [ [x, y, z] for x in np.linspace( -i.input_data["box_size"], i.input_data["box_size"], N )
					   for y in np.linspace( -i.input_data["box_size"], i.input_data["box_size"], N )
					   for z in np.linspace( -i.input_data["box_size"], i.input_data["box_size"], N ) ]


	for point in grid:
		radius = 0.5 * dx
		with Sphere(point, radius) as probe:

			if overlap( probe, populated_box, 0.0 ):
				print(1)
				continue

			else:
                versors = [ np.array([nx * i.input_data["box_size"][0],
                                      ny * i.input_data["box_size"][1],
                                      nz * i.input_data["box_size"][2]])
                            for nx in np.arange(-1, 2, 1)
                            for ny in np.arange(-1, 2, 1)
                            for nz in np.arange(-1, 2, 1) ]

                if_overlap = False

                for versor in versors:

                    probe.translate( versor )

                    if overlap(probe, populated_box, 0.0):
                        if_overlap = True

                    probe.translate( -versor )

                if if_overlap: print(1)
            	else: print(0)

	# timestamp( 'Computing radial distribution function' )
	# bins, rdfs, distinct_rdf = compute_rdfs( i.input_data, labels, auxiliary_data )

	# timestamp( 'Saving radial distribution function to a file' )
	# save_rdfs_to_file(i.input_data, bins, rdfs, distinct_rdf)

	del times
	# del msds

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
