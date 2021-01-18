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
from pyBrown.messaging import timestamp

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided]
	required_keywords = ["input_xyz_filename", "input_str_filename", "output_str_filename",
						 "snapshot_time"]

	timestamp( 'Reading input from {} file', input_filename )
	i = InputData(input_filename, required_keywords, {}).input_data

	print(i)

	xyz_file = open(i["input_xyz_filename"], 'r')
	template_file = open(i["input_str_filename"], 'r')
	str_file = open(i["output_str_filename"], 'w')

	labels = []
	coords = []

	start_snapshot = False

	for line in xyz_file:
		if start_snapshot:
			if 'time' in line.split():
				start_snapshot = False
			else:
				if len( line.split() ) == 4:
					labels.append( line.split()[0] )
					coords.append( [ line.split()[i] for i in range(1,4) ] )
		if 'time' in line.split():
			if float(line.split()[-1]) == i["snapshot_time"]:
				start_snapshot = True

	for j, line in enumerate( template_file ):
		if line.split()[0] == 'sub':
			line_elements = line.split()
			line_elements[3] = coords[j][0]
			line_elements[4] = coords[j][1]
			line_elements[5] = coords[j][2]
			new_line = ''
			for i in range( len(line_elements) ): new_line += line_elements[i] + ' '
			new_line += '\n'
			str_file.write(new_line)
		else:
			str_file.write(line)


	xyz_file.close()
	template_file.close()
	str_file.close()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	main()
