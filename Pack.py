# pyBrown is a bound of tools useful for Brownian and Stokesian dynamics simulations
# Copyright (C) 2018  Tomasz Skora (tskora@ichf.edu.pl)
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
# along with this program.  If not, see https://www.gnu.org/licenses.

import click

from pyBrown.input import InputData
from pyBrown.grid import pack_molecules
from pyBrown.write import write_structure
from pyBrown.parse import parse_input_filename
from pyBrown.messaging import timestamp

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
				"bond_lengths":'hydrodynamic_radii'}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputData(input_filename, required_keywords, defaults).input_data

	### TO BE REFACTORED
	if not isinstance( i["max_bond_lengths"], list ):
		max_bond_lengths = [ ]
		for bfc in i["bond_force_constants"]:
			max_bond_lengths.append( [ i["max_bond_lengths"]
									for j in range( len( bfc ) ) ] )
		i["max_bond_lengths"] = max_bond_lengths

	if i["bond_lengths"] == 'hydrodynamic_radii':
		bond_lengths = [ ]
		for j, bfc in enumerate( i["bond_force_constants"] ):
			bond_lengths.append( [ ] )
			for k in range( len( bfc ) ):
				bond_lengths[j].append( i["hydrodynamic_radii"][j][k] +
										i["hydrodynamic_radii"][j][k+1] )
		i["bond_lengths"] = bond_lengths
	###

	coords = pack_molecules(i)

	write_structure(i, coords)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
