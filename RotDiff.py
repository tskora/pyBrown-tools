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

from pyBrown.input_RotDiff import InputDataRotDiff
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads, \
								 compute_orientations, compute_mean_orientation_autocorrelation, \
								 save_mean_orientation_autocorrelation_to_file, \
								 compute_mean_squared_angular_displacements, \
								 save_mean_squared_angular_displacements_to_file, \
								 plot_msads, plot_moas

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "input_xyz_template", "input_xyz_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "probing_frequency": 1, "fit_MOA": False,
				"min_time": 0.0, "mode": "window", "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataRotDiff(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )
	timestamp( 'Separating the rotational movement' )
	orientation_labels = compute_orientations( i.input_data, labels, auxiliary_data )
	del labels

	timestamp( 'Computing mean orientation autocorrelations' )
	moas = compute_mean_orientation_autocorrelation( i.input_data, orientation_labels, auxiliary_data )
	timestamp( 'Saving mean orientation autocorrelations to a file' )
	save_mean_orientation_autocorrelation_to_file(i.input_data, times, moas)

	timestamp( 'Computing mean squared angular displacements' )
	msads = compute_mean_squared_angular_displacements( i.input_data, orientation_labels, auxiliary_data )
	timestamp( 'Saving mean squared angular displacements to a file' )
	save_mean_squared_angular_displacements_to_file(i.input_data, times, msads)
	del orientation_labels

	timestamp( 'Plotting mean squared angular displacements' )
	plot_msads(i.input_data, times, msads)

	timestamp( 'Plotting mean orientation autocorrelations' )
	plot_moas(i.input_data, times, moas)

	del times
	del moas
	del msads

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
