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

from pyBrown.input_Energy import InputDataEnergy
from pyBrown.messaging import timestamp
from pyBrown.energies import read_energies, compute_menergies, plot_energies

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "input_enr_template", "input_enr_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "fit_MSD": False,
				"probing_frequency": 1, "mode": "window", "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataEnergy(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading energies' )
	times, energies = read_energies(i.input_data)

	timestamp( 'Averaging energies' )
	menergies = compute_mean_energies(i.input_data, energies)
	del energies

	# timestamp( 'Saving energies to a file' )
	# save_energies_to_file(i.input_data, times, energies)

	timestamp( 'Plotting energies' )
	plot_energies(i.input_data, times, energies)

	del times
	del menergies

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
