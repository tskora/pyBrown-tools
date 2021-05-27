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
import os

from pyBrown.input_MSDD import InputDataMSDD
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads, \
								 separate_center_of_mass, \
								 compute_msdds, \
								 save_msdds_to_file
# from pyBrown.plotting import plot_msds

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "temperature", "viscosity",
						 "input_xyz_template", "input_xyz_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "fit_MSD": False,
				"probing_frequency": 1, "min_time": 0.0, "mode": "window",
				"float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataMSDD(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )
	timestamp( 'Separating the center of mass movement' )
	cm_labels = separate_center_of_mass( i.input_data, labels, auxiliary_data )
	del labels

	timestamp( 'Computing mean squared distance displacements' )
	msdds = compute_msdds( i.input_data, cm_labels, auxiliary_data )
	del cm_labels

	timestamp( 'Saving mean squared distance displacements to a file' )
	save_msdds_to_file(i.input_data, times, msdds)

	# timestamp( 'Plotting mean squared distance displacements' )
	# plot_msdds(i.input_data, times, msdds)

	del times
	del msdds

	timestamp( 'Deleting binary files' )
	delete_binary_files(auxiliary_data)

#-------------------------------------------------------------------------------

def delete_binary_files(aux):
		for key in aux.keys():
			if "temp_filename" in key: os.remove(aux[key])

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
