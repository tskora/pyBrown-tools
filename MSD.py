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

from pyBrown.input_MSD import InputDataMSD
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, read_energies, \
								 separate_center_of_mass, \
								 compute_msds, compute_menergies, \
								 save_msds_to_file, \
								 plot_msds, plot_menergies

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
				"probing_frequency": 1, "min_time": 0.0}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataMSD(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	if "input_enr_template" in i.input_data.keys() and \
	   "input_enr_range" in i.input_data.keys():

		timestamp( 'Reading energies')
		energies, times = read_energies(i.input_data)
		timestamp( 'Computing mean energies' )
		menergies = compute_menergies(energies)
		del energies
		timestamp( 'Plotting mean energies' )
		plot_menergies(i.input_data, times, menergies)
		del times
		del menergies

	timestamp( 'Reading trajectories' )
	temporary_filename, times, labels, min_time_index = read_trajectories(i.input_data)
	timestamp( 'Separating the center of mass movement' )
	temporary_filename_2, cm_labels = separate_center_of_mass(i.input_data, temporary_filename, labels)
	del labels
	timestamp( 'Computing mean square displacements' )
	msds = compute_msds(i.input_data, temporary_filename_2, cm_labels, min_time_index)
	del cm_labels
	timestamp( 'Saving mean square displacements to a file' )
	save_msds_to_file(i.input_data, times, msds)
	timestamp( 'Plotting mean square displacements' )
	plot_msds(i.input_data, times, msds)
	del times
	del msds

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
