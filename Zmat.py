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

from pyBrown_tools.input_Lengths import InputDataLengths
from pyBrown_tools.messaging import timestamp
from pyBrown_tools.trajectories import read_trajectories, add_auxiliary_data_multibeads, \
								 compute_zmat

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "input_xyz_template", "input_xyz_range",
						 "internal_coordinates"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "probing_frequency": 1,
				"min_time": 0.0, "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataLengths(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )
	timestamp( 'Computing zmat' )
	zmat_labels = compute_zmat( i.input_data, labels, auxiliary_data )
	# timestamp( 'Computing angle distribution' )
	# bins, angle_distribution = compute_angle_distribution( i.input_data, angle_labels, auxiliary_data )
	# save_angles_to_file( i.input_data, bins, angle_distribution )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
