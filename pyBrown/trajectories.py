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

import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import Boltzmann
from scipy.linalg import eigh

import freud.box
import freud.msd

from pyBrown.messaging import timestamp
from pyBrown.plot_config import plot_config

#-------------------------------------------------------------------------------

def read_trajectories(input_data):

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]
	number_of_timeframes = _number_of_timeframes( input_xyz_filenames, input_data["probing_frequency"] )
	min_time_index = _find_min_time_index( input_xyz_filenames, input_data["min_time"], input_data["probing_frequency"] )
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_beads = _number_of_beads( input_xyz_filenames )
	number_of_beads_per_file = number_of_beads // number_of_xyz_files
	labels = [ '___' for _ in range( number_of_beads ) ]
	number_of_timeframes -= min_time_index
	times = np.zeros( number_of_timeframes, dtype = np.float32 )

	# temporary binary file which will contain the trajectories
	traj_temp_filename = input_data["input_xyz_template"] + 'trj_tmp.dat'

	if input_data["debug"]:
		print( 'Input xyz filenames: {}'.format(input_xyz_filenames) )
	if input_data["debug"] or input_data["verbose"]:
		print( 'Number of timeframes: {}'.format(number_of_timeframes) )
		print( 'Number of beads: {}'.format(number_of_beads) )
		print( 'Number of files: {}'.format(number_of_xyz_files) )

	auxiliary_data = {"input_xyz_filenames": input_xyz_filenames,
					  "number_of_timeframes": number_of_timeframes,
					  "min_time_index": min_time_index,
					  "number_of_beads": number_of_beads,
					  "number_of_beads_per_file": number_of_beads_per_file,
					  "number_of_timeframes": number_of_timeframes,
					  "traj_temp_filename": traj_temp_filename}

	trajectories = np.memmap( traj_temp_filename, dtype = np.float32,
							  mode = 'w+',
							  shape = ( number_of_beads,
									    number_of_timeframes, 3 ) )

	for i in range( number_of_beads ):
		trajectories[i] = np.zeros( (number_of_timeframes, 3), dtype = np.float32 )

	del trajectories

	for i, input_xyz_filename in enumerate(input_xyz_filenames):

		if input_data["verbose"] or input_data["debug"]:
			timestamp( 'XYZ file {} / {}', i + 1, number_of_xyz_files )

		with open(input_xyz_filename, 'r') as input_xyz_file:

			temp = np.memmap( traj_temp_filename, dtype = np.float32,
						   	  shape = ( number_of_beads,
						      			number_of_timeframes, 3 ) )

			trajectories = temp[i * number_of_beads_per_file : (i + 1) * number_of_beads_per_file, :, :]

			del temp

			_read_trajectories_from_xyz_file(input_xyz_file, trajectories, times,
											 				  labels, input_data["probing_frequency"],
											 				  input_data["min_time"])

	times = times[:] - times[0]

	return times, labels, auxiliary_data

#-------------------------------------------------------------------------------

def _second_Legendre_polynomial(orientation, director):

	return 1.5 * np.dot( orientation, director )**2 - 0.5


def _number_of_timeframes( input_xyz_filenames, probing_frequency ):

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  	probing_frequency )
							  	for i in range( number_of_xyz_files ) ]

	return min( numbers_of_timeframes )

def _min_time_index(filename, min_time, probing_frequency):

	with open(filename, 'r') as file:

		min_time_index = 0
		counter = 0

		for line in file:

			if 'xyz' in line.split()[0].split('.'):

				if ( counter % probing_frequency ) == 0:
					time = float( line.split()[3] )

					if time >= min_time:
						return min_time_index
					else:
						min_time_index += 1

				counter += 1

def _find_min_time_index(input_xyz_filenames, min_time, probing_frequency):

	number_of_xyz_files = len( input_xyz_filenames )

	min_time_index = _min_time_index( input_xyz_filenames[0], min_time, probing_frequency )

	for i in range( 1, number_of_xyz_files ):

		assert _min_time_index( input_xyz_filenames[i], min_time, probing_frequency ) == min_time_index

	print(min_time_index)

	return min_time_index

def _count_beads(filename):

	with open(filename, 'r') as file:
		return( int( file.readline() ) )

def _number_of_beads( input_xyz_filenames ):

	number_of_beads = _count_beads ( input_xyz_filenames[0] )

	for i in range( 1, len(input_xyz_filenames) ):

		assert _count_beads( input_xyz_filenames[i] ) == number_of_beads

	return number_of_beads * len( input_xyz_filenames )

def _molecule_sizes( input_labels, input_sizes ):

	molecule_sizes = {}

	for label, size  in zip( input_labels, input_sizes ):
		molecule_sizes[label] = size

	return molecule_sizes

def _molecule_numbers( input_labels, input_sizes, labels, molecule_sizes ):

	molecule_numbers = {}

	for label, size  in zip( input_labels, input_sizes ):
		molecule_numbers[label] = 0

	for label in labels:
		molecule_numbers[label] += 1

	for label in input_labels:
		molecule_numbers[label] = molecule_numbers[label] // molecule_sizes[label]

	return molecule_numbers

def _number_of_molecules( molecule_numbers, input_labels ):

	number_of_molecules = 0

	for label in input_labels:
		number_of_molecules += molecule_numbers[label]

	return number_of_molecules

def add_auxiliary_data_multibeads(input_data, labels, auxiliary_data):

	auxiliary_data["molecule_sizes"] = _molecule_sizes( input_data["labels"], input_data["sizes"] )
	auxiliary_data["molecule_numbers"] = _molecule_numbers( input_data["labels"], input_data["sizes"], labels, auxiliary_data["molecule_sizes"] )
	auxiliary_data["number_of_molecules"] = _number_of_molecules( auxiliary_data["molecule_numbers"], input_data["labels"] )

def _keep_bound_beads_in_the_same_box(r, r_ref, box_size):

	for i in range(0, 3):
		if ( r[i] - r_ref[i] >= box_size / 2 ):
			r[i] -= box_size
		elif ( r_ref[i] - r[i] >= box_size / 2 ):
			r[i] += box_size

def compute_orientations(input_data, labels, auxiliary_data):

	which_trajectory = 0
	which_orientation_trajectory = 0
	number_of_orientation_trajctories = 0

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_beads = auxiliary_data["number_of_beads"]
	number_of_beads_per_file = auxiliary_data["number_of_beads_per_file"]
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	traj_temp_filename = auxiliary_data["traj_temp_filename"]

	print( 'number of orientation trajectories: {}'.format(number_of_orientation_trajectories) )

	# temporary binary file which will contain the orientation trajectories
	orient_temp_filename = input_data["input_xyz_template"] + 'ort_tmp.dat'
	auxiliary_data["orient_temp_filename"] = orient_temp_filename

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes, 3 ) )

	for i in range(number_of_orientation_trajectories):
		orientation_trajectories[i] = np.zeros((number_of_timeframes, 3), dtype = np.float32)

	del orientation_trajectories

	orientation_labels = [ '___' for i in range(number_of_orientation_trajectories) ]

	trajectories = np.memmap( traj_temp_filename, dtype = np.float32,
						shape = ( number_of_beads, number_of_timeframes, 3 ) )

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						shape = ( number_of_orientation_trajectories, number_of_timeframes, 3 ) )

	while( which_trajectory < number_of_beads ):

		multiplicity = molecule_sizes[ labels[ which_trajectory ] ]

		if multiplicity == 1:

			orientation_labels[ which_orientation_trajectory ] = labels[ which_trajectory ]

			which_trajectory += 1
			which_orientation_trajectory += 1

		else:

			for i in range( number_of_timeframes ):

				r_ref = trajectories[which_trajectory, i, :]

				for j in range( 1, multiplicity ):

					r = trajectories[which_trajectory + j, i, :]

					_keep_bound_beads_in_the_same_box( r, r_ref, input_data["box_size"] )

					r_ref = r

					trajectories[which_trajectory + j, i, :] = np.array(r, dtype = np.float32)

			orientations = trajectories[which_trajectory + multiplicity - 1, :, :] - trajectories[which_trajectory, :, :]

			orientation_norms = np.linalg.norm(orientations, axis = 1)

			for v, vn in zip(orientations, orientation_norms):
				v /= vn
			
			orientation_trajectories[which_orientation_trajectory, :] = orientations

			orientation_labels[ which_orientation_trajectory ] = labels[ which_trajectory ]

			which_trajectory += multiplicity
			which_orientation_trajectory += 1

	del trajectories
	del orientation_trajectories

	if input_data["verbose"]: print('orientation computation performed')

	return orientation_labels

def compute_mean_orientation_autocorrelation(input_data, orientation_labels, auxiliary_data):

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	traj_temp_filename = auxiliary_data["traj_temp_filename"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	# temporary binary file which will contain the orientation autocorrelations
	oa_temp_filename = input_data["input_xyz_template"] + 'oa_tmp.dat'
	auxiliary_data["oa_temp_filename"] = oa_temp_filename

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	oas = np.memmap( oa_temp_filename, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes ) )

	for i in range(number_of_orientation_trajectories):
		oas[i] = np.zeros(number_of_timeframes, dtype = np.float32)

	moa = [ np.zeros(number_of_timeframes) for i in range(len(input_data["sizes"])) ]
	amounts = np.zeros(len(input_data["sizes"]))

	del orientation_trajectories

	del oas

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing autocorrelation: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   shape = ( number_of_orientation_trajectories,
						   	number_of_timeframes, 3 ) )

		orientation_trajectory = orientation_trajectories[i]

		del orientation_trajectories

		autocorrelation = _compute_autocorrelation(orientation_trajectory, mode = input_data["mode"])

		oas = np.memmap( oa_temp_filename, dtype = np.float32,
					 	   shape = ( len(orientation_labels), number_of_timeframes ) )

		oas[i] = autocorrelation

		del oas

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('averaging: {} / {}', i + 1, number_of_orientation_trajectories)

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		oas = np.memmap( oa_temp_filename, dtype = np.float32,
					 	   shape = ( number_of_orientation_trajectories, number_of_timeframes ) )

		autocorrelation = oas[i]

		moa[counter] += autocorrelation

		del oas

	for _moa, label in zip( moa, input_data["labels"] ):

		_moa /= molecule_numbers[ label ]

	return moa

def compute_mean_director( input_data, orientation_labels, auxiliary_data ):

	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	preQ = np.zeros( (number_of_orientation_trajectories, 3, 3), dtype = np.float32 )

	directors = np.zeros( (number_of_orientation_trajectories, 3), dtype = np.float32 )

	mean_director = [ np.zeros(3, dtype = np.float32) for i in range(len(input_data["sizes"])) ]

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing directors: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectory = orientation_trajectories[i]

		for j in range( number_of_timeframes ):

			for k in range(3):

				for l in range(3):

					preQ[i][k][l] += 1.5 * orientation_trajectory[j][k] * orientation_trajectory[j][l] - 0.5 * np.identity(3)[k][l]

		preQ[i] /= number_of_timeframes

	for i in range( number_of_orientation_trajectories ):

		_, v = eigh( preQ[i], eigvals = (2, 2) )
		directors[i] = np.transpose(v)[0]

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('averaging directors: {} / {}', i + 1, number_of_orientation_trajectories)

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		mean_director[counter] += directors[i]

	for i, label in enumerate( input_data["labels"] ):

		mean_director[i] /= molecule_numbers[ label ]

		mean_director[i] /=np.linalg.norm( mean_director[i] )

	return mean_director

def compute_nematic_order( input_data, orientation_labels, auxiliary_data, director ):

	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	mean_second_Legendre_polynomial = np.zeros( number_of_orientation_trajectories, dtype = np.float32 )
	nematic_order_parameter = np.zeros( len(input_data["sizes"]) )

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	amounts = np.zeros(len(input_data["sizes"]))

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('averaging second order Legendre polynomial over trajectory: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectory = orientation_trajectories[i]

		for j in range( number_of_timeframes ):

			mean_second_Legendre_polynomial[i] += _second_Legendre_polynomial(orientation_trajectory[j], director)

	mean_second_Legendre_polynomial /= number_of_timeframes
		
	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing nematic order parameter: {} / {}', i + 1, number_of_orientation_trajectories)

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		nematic_order_parameter[counter] += mean_second_Legendre_polynomial[i]

	for i, label in enumerate( input_data["labels"] ):

		nematic_order_parameter[i] /= molecule_numbers[ label ]

	return nematic_order_parameter

def save_mean_orientation_autocorrelation_to_file(input_data, times, moas):

	output_filename = input_data["input_xyz_template"] + 'moa.txt'

	with open(output_filename, 'w') as output_file:

		first_line = 'time/ps '
		line = '{} '

		for label in input_data["labels"]:

			first_line += ( label + ' ' )

			line += '{} '

		output_file.write(first_line + '\n')

		for i in range( len(times) ):

			line_values = [ times[i] ]

			for j in range( len(input_data["labels"]) ):

				line_values.append( moas[j][i] )

			output_file.write( line.format(*line_values) + '\n' )

def _compute_autocorrelation(orientation, mode = 'direct'):

	import numpy as np

	if mode == 'direct':

		return np.array( [ np.dot( orientation[0], orientation[j] ) for j in range(len(orientation)) ], np.float32 )

	elif mode == 'window':

		return np.array( [1.0] + [ np.mean( [ np.dot(orientation[k], orientation[k+j]) for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ], np.float32 )

def save_mean_squared_angular_displacements_to_file(input_data, times, msads):

	output_filename = input_data["input_xyz_template"] + 'msad.txt'

	with open(output_filename, 'w') as output_file:

		first_line = 'time/ps '
		line = '{} '

		for label in input_data["labels"]:

			first_line += ( label + ' ' )

			line += '{} '

		output_file.write(first_line + '\n')

		for i in range( len(times) ):

			line_values = [ times[i] ]

			for j in range( len(input_data["labels"]) ):

				line_values.append( msads[j][i] )

			output_file.write( line.format(*line_values) + '\n' )

#-------------------------------------------------------------------------------

def compute_mean_squared_angular_displacements(input_data, orientation_labels, auxiliary_data):

	import numpy as np

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   shape = ( number_of_orientation_trajectories, number_of_timeframes, 3 ) )

	# temporary binary file which will contain the squared angular displacements
	sad_temp_filename = input_data["input_xyz_template"] + 'sad_tmp.dat'
	auxiliary_data["sad_temp_filename"] = sad_temp_filename

	sads = np.memmap( sad_temp_filename, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes ) )

	for i in range(number_of_orientation_trajectories):
		sads[i] = np.zeros(number_of_timeframes, dtype = np.float32)

	msad = [ np.zeros(number_of_timeframes) for i in range(len(input_data["sizes"])) ]

	del orientation_trajectories

	del sads

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing sad: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectories = np.memmap( orient_temp_filename, dtype = np.float32,
						   shape = ( number_of_orientation_trajectories,
						   			 number_of_timeframes, 3 ) )

		orientation_trajectory = orientation_trajectories[i]

		del orientation_trajectories

		sad = _compute_sad(orientation_trajectory, mode = input_data["mode"])


		sads = np.memmap( sad_temp_filename, dtype = np.float32,
					 	   shape = ( number_of_orientation_trajectories, number_of_timeframes ) )

		sads[i] = sad

		del sads

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('averaging: {} / {}', i + 1, number_of_orientation_trajectories)

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		sads = np.memmap( sad_temp_filename, dtype = np.float32,
					 	   shape = ( number_of_orientation_trajectories, number_of_timeframes ) )

		sad = sads[i]

		msad[counter] += sad

		del sads

	print(msad)

	print(molecule_numbers)

	for _msad, label in zip( msad, input_data["labels"] ):

		_msad /= molecule_numbers[ label ]

	print(msad)

	print(molecule_numbers)

	return msad

def _compute_sad(orientation, mode = 'direct'):

	import numpy as np

	if mode == 'direct':

		return np.array( [ np.arccos( np.dot( orientation[0], orientation[j] ) )**2 for j in range(len(orientation)) ], np.float32 )

	elif mode == 'window':

		return np.array( [0.0] + [ np.mean( [ np.arccos( np.dot(orientation[k], orientation[k+j]) )**2 for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ], np.float32 )

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

	return np.array(energies), np.array(times, float)

#-------------------------------------------------------------------------------

def separate_center_of_mass(input_data, temporary_filename, labels):

	bead_dict = {}
	mol_dict = {}

	which_trajectory = 0
	which_cm_trajectory = 0

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  input_data["probing_frequency"] )
							  for i in range( number_of_xyz_files ) ]

	number_of_timeframes = min( numbers_of_timeframes )

	number_of_cm_trajectories = 0

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

		mol_dict[label] = 0

	for label in labels:

		mol_dict[label] += 1

	for label in input_data["labels"]:

		mol_dict[label] = mol_dict[label] // bead_dict[label]

		number_of_cm_trajectories += mol_dict[label]

	print( 'number of cm trajectories: {}'.format(number_of_cm_trajectories) )

	temporary_filename_2 = 'temp2.dat'

	cm_trajectories = np.memmap( temporary_filename_2, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_cm_trajectories,
									  number_of_timeframes, 3 ) )

	for i in range(number_of_cm_trajectories):
		for j in range(number_of_timeframes):
			cm_trajectories[i, j, :] = np.zeros(3, dtype = np.float32)

	del cm_trajectories

	cm_labels = [ '___' for i in range(number_of_cm_trajectories) ]

	while( which_trajectory < len( labels ) ):

		multiplicity = bead_dict[ labels[ which_trajectory ] ]

		if multiplicity == 1:

			temp = np.memmap( temporary_filename, dtype = np.float32,
							shape = ( len(labels), number_of_timeframes, 3 ) )

			temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
							shape = ( len(cm_labels), number_of_timeframes, 3 ) )

			temp2[which_cm_trajectory, :, :] = temp[which_trajectory, :, :]
			cm_labels[ which_cm_trajectory ] = labels[ which_trajectory ]

			which_trajectory += 1
			which_cm_trajectory += 1

			del temp
			del temp2

		else:

			temp = np.memmap( temporary_filename, dtype = np.float32,
							shape = ( len(labels), number_of_timeframes, 3 ) )

			temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
							shape = ( len(cm_labels), number_of_timeframes, 3 ) )

			# temp2[which_cm_trajectory, :, :] = temp[which_trajectory, :, :]

			# for i in range( 1, multiplicity ):

			# 	temp2[which_cm_trajectory, :, :] += temp[which_trajectory + i, :, :]

			# temp2[which_cm_trajectory, :, :] /= multiplicity

			###{

			temp2[which_cm_trajectory, :, :] = temp[which_trajectory, :, :]

			# print('bead')

			for i in range( number_of_timeframes ):

				# print('i: {}'.format(i))

				x, y, z = temp[which_trajectory, i, :]

				for j in range( 1, multiplicity ):

					# print('j: {}'.format(j))

					xi, yi, zi = temp[which_trajectory + j, i, :]
					shift = False

					if ( xi - x >= input_data["box_size"] / 2 ):
						xi -= input_data["box_size"]
						shift = True
					elif ( x - xi >= input_data["box_size"] / 2 ):
						xi += input_data["box_size"]
						shift = True

					if ( yi - y >= input_data["box_size"] / 2 ):
						yi -= input_data["box_size"]
						shift = True
					elif ( y - yi >= input_data["box_size"] / 2 ):
						yi += input_data["box_size"]
						shift = True

					if ( zi - z >= input_data["box_size"] / 2 ):
						zi -= input_data["box_size"]
						shift = True
					elif ( z - zi >= input_data["box_size"] / 2 ):
						zi += input_data["box_size"]
						shift = True

					# if shift: print('xi {} yi {} zi {}'.format(xi, yi, zi))

					x, y, z = xi, yi, zi

					temp[which_trajectory + j, i, :] = np.array([xi, yi, zi], dtype = np.float32)

			temp2[which_cm_trajectory, :, :] = temp[which_trajectory, :, :]

			for tframe in range(number_of_timeframes):

				print( np.sqrt( np.sum(( temp[which_trajectory + multiplicity - 1,tframe,:] - temp[which_trajectory,tframe,:] )**2) ) )

			print('')

			for i in range( 1, multiplicity ):

				temp2[which_cm_trajectory, :, :] += temp[which_trajectory + i, :, :]

			temp2[which_cm_trajectory, :, :] /= multiplicity

			#
			print( temp2[which_cm_trajectory, :, :] )
			#

			# for i in range( 1, multiplicity ):

			# 	for j in range(number_of_timeframes):

			# 		temp[which_trajectory + i, j, :]

			# for i in range(multiplicity):

			# 	for j in range(number_of_timeframes):
			# 		dist = np.sqrt( np.sum( ( temp2[which_cm_trajectory, j, :] - temp[which_trajectory + i, j, :] )**2 ) )
			# 		if dist > 160.0:
			# 			print('distance between cm and bead')
			# 			print( dist )
			# 			print( temp2[which_cm_trajectory, j, :] )
			# 			print( temp[which_trajectory + i, j, :] )

			###}

			cm_labels[ which_cm_trajectory ] = labels[ which_trajectory ]

			which_trajectory += multiplicity
			which_cm_trajectory += 1

			del temp
			del temp2

	#
	1/0
	#
	if input_data["verbose"]: print('cm separation performed')

	os.remove(temporary_filename)

	return temporary_filename_2, cm_labels

#-------------------------------------------------------------------------------

def compute_msds(input_data, temporary_filename_2, cm_labels, min_time_index):

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  input_data["probing_frequency"] )
							  for i in range( number_of_xyz_files ) ]

	number_of_timeframes = min( numbers_of_timeframes )

	bead_dict = {}

	temporary_filename_3 = 'temp3.dat'

	sds = np.memmap( temporary_filename_3, dtype = np.float32,
					 mode = 'w+',
					 shape = ( len(cm_labels), number_of_timeframes - min_time_index ) )

	if input_data["verbose"]: print("filling sds")

	for i in range( len(cm_labels) ):
		sds[i, :] = np.zeros( number_of_timeframes - min_time_index, dtype = np.float32 ) 

	del sds

	if input_data["verbose"]: print("end filling sds")

	msds = np.zeros( ( len( input_data["labels"] ), number_of_timeframes - min_time_index ), float )
	amounts = np.zeros( len( input_data["labels"] ), float )

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	### FREUD VERSION START

	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(cm_labels), number_of_timeframes, 3 ) )

	unify_coordinates(temp2, input_data["box_size"])

	trajs_all = []

	for k in range(len(input_data["labels"])):

		trajs_list = [ [ temp2[i][j + min_time_index] for i in range( len(cm_labels) ) if input_data["labels"][k] == cm_labels[i]  ] for j in range( number_of_timeframes - min_time_index ) ]

		trajs = np.array( trajs_list )

		trajs_all.append(trajs)

	msds = []

	for k in range(len(input_data["labels"])):

		box = freud.box.Box.cube( input_data["box_size"] )

		if (input_data["mode"] == "direct"):
			msd = freud.msd.MSD(box, 'direct')
		elif (input_data["mode"] == "window"):
			msd = freud.msd.MSD(box, 'window')

		msd.compute( positions = trajs_all[k] )
		msds.append(msd.msd)

	### FREUD VERSION END

	### OLD VERSION START

	# del trajs

	# del temp2

	###

	# for i in range( len(cm_labels) ):

	# 	if input_data["verbose"]: print('computing sd: {} / {}'.format(i, len(cm_labels)))

	# 	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
	# 					   shape = ( len(cm_labels), number_of_timeframes, 3 ) )

	# 	cm_trajectory = temp2[i]

	# 	del temp2

	# 	sd = _compute_sd(cm_trajectory, input_data["box_size"])

	# 	print(sd)

	# 	temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
	# 				 	   shape = ( len(cm_labels), number_of_timeframes ) )

	# 	temp3[i] = sd

	# 	del temp3

	# for i in range( len(cm_labels) ):

	# 	if input_data["verbose"]: print('averaging: {} / {}'.format(i, len(cm_labels)))

	# 	counter = 0

	# 	for input_label in input_data["labels"] :

	# 		if cm_labels[i] == input_label:

	# 			break

	# 		counter += 1

	# 	temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
	# 				 	   shape = ( len(cm_labels), number_of_timeframes ) )

	# 	sds = temp3[i]

	# 	msds[counter] += sds

	# 	del temp3

	# 	amounts[counter] += 1

	# for msd, amount in zip( msds, amounts ):

	# 	msd /= amount

	### OLD VERSION END

	os.remove(temporary_filename_2)
	os.remove(temporary_filename_3)

	return msds

#-------------------------------------------------------------------------------

def compute_menergies(energies):

	menergies = np.zeros( len( energies[0] ), float )

	for energies_one_file in energies:

		menergies += energies_one_file

	return menergies / len( energies )

#-------------------------------------------------------------------------------

def save_msds_to_file(input_data, times, msds):

	output_filename = input_data["input_xyz_template"] + 'msd.txt'

	with open(output_filename, 'w') as output_file:

		first_line = 'time/ps '
		line = '{} '

		for label in input_data["labels"]:

			first_line += ( label + ' ' )

			line += '{} '

		output_file.write(first_line + '\n')

		for i in range( len(times) ):

			line_values = [ times[i] ]

			for j in range( len(input_data["labels"]) ):

				line_values.append( msds[j][i] )

			output_file.write( line.format(*line_values) + '\n' )

#-------------------------------------------------------------------------------

def plot_msds(input_data, times, msds):

	colors, _ = plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'MSD [$\AA ^2$]')

	# plt.yticks([])

	for i, msd in enumerate(msds):

		plt.plot( times / 1000000, msd, '-', label = input_data["labels"][i], color = colors[i] )

		if input_data["fit_MSD"]:
			a, b = np.polyfit(times, msd, 1)
			# print('a = {}; b = {}'.format(a, b))
			D = a / 6 # in A**2 / ps
			# print('D = {}'.format(D))
			D_string = '{:7.5f}'.format(D * 1000) # in AA**2 / ns
			Rh = 10**17 * Boltzmann / (6 * np.pi) * input_data["temperature"] / ( D * 0.1 * input_data["viscosity"] ) # in nm
			Rh_string = '{:4.2f}'.format(Rh) # in nm

			plt.text( 0.0 * times[-1], (0.70 - i * 0.15) * np.ndarray.max( np.array( msds ) ), input_data["labels"][i] )
			plt.text( 0.0 * times[-1], (0.65 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$D =$' + D_string + r' $\frac{\AA ^2}{ns}$')
			plt.text( 0.0 * times[-1], (0.60 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$R_H =$' + Rh_string + ' nm')

			plt.plot(times / 1000000, a * np.array(times, float) + b * np.ones(len(times)), '--', label = 'linear fit for ' + input_data["labels"][i], color = colors[i]  )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'msd.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def plot_menergies(input_data, times, menergies):

	plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'E [$\frac{kcal}{mol}$]')

	plt.plot(times / 1000000, menergies, '-')

	plt.savefig(input_data["input_enr_template"] + 'enr.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def _read_trajectories_from_xyz_file(xyz_file, trajectories, times, labels, probing_frequency, min_time):

	number_of_beads_per_file = int( xyz_file.readline() )

	counter = 0
	min_time_index = 0

	label_index = 0
	for label in labels:
		if label != '___':
			label_index += 1
		else:
			break

	for i, line in enumerate(xyz_file):

		timeframe_index = counter // number_of_beads_per_file // probing_frequency

		bead_index = counter % number_of_beads_per_file

		is_included_as_timeframe_given_frequency = ( ( counter // number_of_beads_per_file ) % probing_frequency ) == 0

		is_line_unimportant = len( line.split() ) == 1

		if timeframe_index - min_time_index >= len(times): break

		if 'xyz' in line.split()[0].split('.'):

			if is_included_as_timeframe_given_frequency:

				time = float( line.split()[3] )

				if time >= min_time:

					times[ timeframe_index - min_time_index ] = time

				else:

					min_time_index += 1

		elif is_line_unimportant: continue

		else:

			if counter < number_of_beads_per_file:
				labels[label_index + counter] = line.split()[0]
			if is_included_as_timeframe_given_frequency:
				if time >= min_time:
					coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ], float )
					trajectories[bead_index, timeframe_index - min_time_index] = coords
			counter += 1

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

def _compute_sd(trajectory, box_size):

    sd = np.zeros( len(trajectory) )
    
    for i in range( 1, len(trajectory) ):

        jump = np.sqrt( np.sum( ( trajectory[i] - trajectory[i - 1] )**2 ) )
        
        if ( jump >= box_size * np.sqrt(3) / 2 ):
            # print('i: {}'.format(i))
            # print('present trajectory: {}'.format(trajectory[i]))
            # print('previous trajectory: {}'.format(trajectory[i-1]))
            # print('jump: {}'.format(jump))

            versors = box_size * np.array( [ [i, j, k]
                for i in [-1,0,1] for j in [-1,0,1] for k in [-1,0,1]
                if not ( i == 0 and j == 0 and k == 0 ) ], dtype = np.float32 )

            for versor in versors:
                # print('versor: {}'.format(versor))
                replica = versor + trajectory[i]
                # print('replica: {}'.format(replica))
                replica_jump = np.sqrt( np.sum( ( replica - trajectory[i - 1] )**2 ) )
                # print('replica jump: {}'.format(replica_jump))
                if replica_jump <= jump:
                    versor_chosen = versor
                    jump = replica_jump
            # print('revised jump: {}'.format(jump))
            # print('')

            assert jump <= box_size * np.sqrt(3) / 2, 'jump {} is impossibly large'.format(i)  

            for j in range(i, len(trajectory)):
                trajectory[j] += versor_chosen
        
        sd[i] = np.sum( ( trajectory[i] - trajectory[0] )**2 )
                          
    return sd

#-------------------------------------------------------------------------------

def _file_length(filename):

	with open(filename, 'r') as file:

		counter = 0

		for line in file:
			counter += 1

		return counter

#-------------------------------------------------------------------------------

def _count_beads(filename):

	with open(filename, 'r') as file:
		return( int( file.readline() ) )

#-------------------------------------------------------------------------------

def _count_timeframes(filename, frequency):

	number_of_beads = _count_beads(filename)
	file_length = _file_length(filename)

	if frequency == 1:
		return file_length // (number_of_beads + 2)
	else:
		return file_length // (number_of_beads + 2) // frequency

#-------------------------------------------------------------------------------

def unify_coordinates(trajectories, box_length):

	for particle in trajectories:
		for i in range(1, len(particle)):
			_coord_unify(particle[i-1], particle[i], box_length)

#-------------------------------------------------------------------------------

def _coord_unify(past, present, box_length):

	for i in range(3):
		while( present[i] - past[i] >= box_length / 2 ):
			present[i] -= box_length
		while( past[i] - present[i] >= box_length / 2 ):
			present[i] += box_length