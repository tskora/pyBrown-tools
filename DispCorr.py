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

from pyBrown.input_DispCorr import InputDataDispCorr
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads, \
								 separate_center_of_mass, compute_mean_displacement_autocorrelation, \
								 save_mean_displacement_autocorrelation_to_file

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "input_xyz_template", "input_xyz_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "probing_frequency": 1,
				"min_time": 0.0, "mode": "window", "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataDispCorr(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )
	timestamp( 'Separating the center of mass movement' )
	cm_labels = separate_center_of_mass( i.input_data, labels, auxiliary_data )
	del labels

	timestamp( 'Computing mean displacement autocorrelations' )
	mdas = compute_mean_displacement_autocorrelation( i.input_data, cm_labels, auxiliary_data )
	timestamp( 'Saving mean displacement autocorrelations to a file' )
	save_mean_displacement_autocorrelation_to_file(i.input_data, times, mdas)

	del times
	del mdas

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
