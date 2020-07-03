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

import numpy as np
import matplotlib.pyplot as plt

from pyBrown.messaging import timestamp

#-------------------------------------------------------------------------------

def read_energies(input_data):

	input_enr_filenames = [ input_data["input_enr_template"] + str(i) + '.enr'
							for i in range( *input_data["input_enr_range"] ) ]

	number_of_enr_files = len( input_enr_filenames )

	numbers_of_timeframes = [ _file_length(input_enr_filenames[i]) - 1
							  for i in range( number_of_enr_files ) ]

	number_of_timeframes = min(numbers_of_timeframes)

	energies = []
	times = []

	for i, input_enr_filename in enumerate(input_enr_filenames):

		if input_data["verbose"]: timestamp('enr file {}/{}'.format( i + 1, len(input_enr_filenames) ))

		with open(input_enr_filename, 'r') as input_enr_file:

			_read_energies_from_enr_file(input_enr_file, energies, times)

	if input_data["verbose"]: timestamp('enrs read')

	for i in range(number_of_enr_files):

		energies[i] = energies[i][:number_of_timeframes]

	times = times[:number_of_timeframes]

	return np.array(times, input_data["float_type"]), np.array(energies, input_data["float_type"])

#-------------------------------------------------------------------------------

def compute_mean_energies(input_data, energies):

	menergies = np.zeros( len( energies[0] ), dtype = input_data["float_type"] )

	for energies_one_file in energies:

		menergies += energies_one_file

	return menergies / len( energies )

#-------------------------------------------------------------------------------

def _read_energies_from_enr_file(enr_file, energies, times):

	energies_from_current_file = []

	if len(times) == 0:
		add_times = True
	else:
		add_times = False
		len_i = len(times)

	enr_file.readline()

	for i, line in enumerate(enr_file):

		if add_times:
			times.append( float( line.split()[0] ) )

		else:
			if i >= len_i:
				break
			else:
				assert( float( line.split()[0] ) == times[i] )

		energies_from_current_file.append( float( line.split()[1] ) )

	energies.append( np.array( energies_from_current_file, float ) )

#-------------------------------------------------------------------------------

def save_menergies_to_file(input_data, times, energies):

	output_filename = input_data["input_enr_template"] + 'enr.txt'

	with open(output_filename, 'w') as output_file:

		first_line = 'time/ps energy/kcalmol-1\n'

		output_file.write(first_line)

		for i in range( len(times) ):

			output_file.write( '{} {}\n'.format(times[i], energies[i]) )

#-------------------------------------------------------------------------------

def _file_length(filename):

	with open(filename, 'r') as file:

		counter = 0

		for line in file:
			counter += 1

		return counter
