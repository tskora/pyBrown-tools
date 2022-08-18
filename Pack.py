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

from pyBrown_tools.input_Pack import InputDataPack
from pyBrown_tools.grid import pack_molecules
from pyBrown_tools.write import write_structure
from pyBrown_tools.parse import parse_input_filename
from pyBrown_tools.messaging import timestamp

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["packing_mode", "box_size", "hydrodynamic_radii",
						 "lennard-jones_radii", "lennard-jones_energies",
						 "charges", "masses", "bond_force_constants",
						 "angle_force_constants", "numbers_of_molecules",
						 "labels_of_molecules", "output_structure_filename"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"minimal_distance_between_surfaces":0.0, "max_bond_lengths":2.5e+07,
				"bond_lengths":'hydrodynamic_radii', "number_of_structures":1,
				"float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataPack(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	for file_count in range(1, i.input_data["number_of_structures"] + 1):

		coords = pack_molecules(i.input_data)

		output_structure_filename = i.input_data["output_structure_filename"].split('.')[0] +\
									'_{}.'.format(file_count) +\
									i.input_data["output_structure_filename"].split('.')[1]

		write_structure(i.input_data, coords, output_structure_filename)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
