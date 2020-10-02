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

from scipy.constants import Boltzmann
from scipy.linalg import eigh

import freud.box
import freud.msd

from pyBrown.messaging import timestamp

FREUD = True
CM = True

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
	times = np.zeros( number_of_timeframes, dtype = input_data["float_type"] )

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

	trajectories = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
							  mode = 'w+',
							  shape = ( number_of_beads,
									    number_of_timeframes, 3 ) )

	for i in range( number_of_beads ):
		trajectories[i] = np.zeros( (number_of_timeframes, 3), dtype = input_data["float_type"] )

	del trajectories

	for i, input_xyz_filename in enumerate(input_xyz_filenames):

		if input_data["verbose"] or input_data["debug"]:
			timestamp( 'XYZ file {} / {}', i + 1, number_of_xyz_files )

		with open(input_xyz_filename, 'r') as input_xyz_file:

			temp = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
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

def add_auxiliary_data_multibeads(input_data, labels, auxiliary_data):

	auxiliary_data["molecule_sizes"] = _molecule_sizes( input_data["labels"], input_data["sizes"] )
	auxiliary_data["molecule_numbers"] = _molecule_numbers( input_data["labels"], input_data["sizes"], labels, auxiliary_data["molecule_sizes"] )
	auxiliary_data["number_of_molecules"] = _number_of_molecules( auxiliary_data["molecule_numbers"], input_data["labels"] )

#-------------------------------------------------------------------------------

def separate_center_of_mass(input_data, labels, auxiliary_data):

	which_trajectory = 0
	which_cm_trajectory = 0
	number_of_cm_trajectories = 0

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_beads = auxiliary_data["number_of_beads"]
	number_of_beads_per_file = auxiliary_data["number_of_beads_per_file"]
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	number_of_cm_trajectories = auxiliary_data["number_of_molecules"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	traj_temp_filename = auxiliary_data["traj_temp_filename"]

	print( 'number of cm trajectories: {}'.format(number_of_cm_trajectories) )

	# temporary binary file which will contain the cm trajectories
	cm_temp_filename = input_data["input_xyz_template"] + 'cm_tmp.dat'
	auxiliary_data["cm_temp_filename"] = cm_temp_filename

	cm_trajectories = np.memmap( cm_temp_filename, dtype = input_data["float_type"],
							mode = 'w+',
							shape = ( number_of_cm_trajectories,
									  number_of_timeframes, 3 ) )

	for i in range(number_of_cm_trajectories):
		cm_trajectories[i] = np.zeros((number_of_timeframes, 3), dtype = input_data["float_type"])

	del cm_trajectories

	cm_labels = [ '___' for i in range(number_of_cm_trajectories) ]

	trajectories = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
						shape = ( number_of_beads, number_of_timeframes, 3 ) )

	cm_trajectories = np.memmap( cm_temp_filename, dtype = input_data["float_type"],
						shape = ( number_of_cm_trajectories, number_of_timeframes, 3 ) )

	while( which_trajectory < len( labels ) ):

		multiplicity = molecule_sizes[ labels[ which_trajectory ] ]

		if multiplicity == 1:

			cm_labels[ which_cm_trajectory ] = labels[ which_trajectory ]

			cm_trajectories[which_cm_trajectory] = trajectories[which_trajectory]

			which_trajectory += 1
			which_cm_trajectory += 1

		else:

			for i in range( number_of_timeframes ):

				r_ref = trajectories[which_trajectory, i, :]

				for j in range( 1, multiplicity ):

					r = trajectories[which_trajectory + j, i, :]

					_keep_bound_beads_in_the_same_box( r, r_ref, input_data["box_size"] )

					r_ref = r

					trajectories[which_trajectory + j, i, :] = np.array(r, dtype = input_data["float_type"])

			cm_trajectories[which_cm_trajectory] = trajectories[which_trajectory, :, :]

			for i in range( 1, multiplicity ):

				cm_trajectories[which_cm_trajectory, :, :] += trajectories[which_trajectory + i, :, :]

			cm_trajectories[which_cm_trajectory, :, :] /= multiplicity

			cm_labels[ which_cm_trajectory ] = labels[ which_trajectory ]

			which_trajectory += multiplicity
			which_cm_trajectory += 1

	if input_data["verbose"]: timestamp('cm separation performed')

	del trajectories

	del cm_trajectories

	return cm_labels

#-------------------------------------------------------------------------------

def compute_msds(input_data, cm_labels, auxiliary_data):

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_cm_trajectories = auxiliary_data["number_of_molecules"]
	cm_temp_filename = auxiliary_data["cm_temp_filename"]

	cm_trajectories = np.memmap( cm_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_cm_trajectories, number_of_timeframes, 3 ) )

	unify_coordinates(cm_trajectories, input_data["box_size"])

	del cm_trajectories

	if input_data["verbose"]: timestamp("filling sds")

	### FREUD VERSION START ###

	if FREUD:

		cm_trajectories = np.memmap( cm_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_cm_trajectories, number_of_timeframes, 3 ) )

		trajs_all = []

		for k in range(len(input_data["labels"])):

			trajs_list = [ [ cm_trajectories[i][j] for i in range( len(cm_labels) ) if input_data["labels"][k] == cm_labels[i]  ] for j in range( number_of_timeframes ) ]

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

	### FREUD VERSION END ###

	else:

		# temporary binary file which will contain the squared angular displacements
		sd_temp_filename = input_data["input_xyz_template"] + 'sd_tmp.dat'
		auxiliary_data["sd_temp_filename"] = sd_temp_filename

		sds = np.memmap( sd_temp_filename, dtype = input_data["float_type"],
								mode = 'w+',
								shape = ( number_of_cm_trajectories,
										  number_of_timeframes ) )

		for i in range(number_of_cm_trajectories):
			sds[i] = np.zeros(number_of_timeframes, dtype = input_data["float_type"])

		msds = [ np.zeros(number_of_timeframes) for i in range(len(input_data["sizes"])) ]

		del sds

		for i in range( number_of_cm_trajectories ):

			if input_data["verbose"]: timestamp('computing sd: {} / {}', i + 1, number_of_cm_trajectories)

			cm_trajectories = np.memmap( cm_temp_filename, dtype = input_data["float_type"],
							 shape = ( number_of_cm_trajectories,
							 			number_of_timeframes, 3 ) )

			cm_trajectory = cm_trajectories[i]

			del cm_trajectories

			# sd = _compute_sd( cm_trajectory, input_data["box_size"], mode = input_data["mode"] )
			sd = _compute_sd( cm_trajectory, input_data["box_size"] )

			sds = np.memmap( sd_temp_filename, dtype = input_data["float_type"],
							shape = ( number_of_cm_trajectories, number_of_timeframes ) )

			sds[i] = sd

			del sds

		for i in range( number_of_cm_trajectories ):

			if input_data["verbose"]: timestamp('averaging: {} / {}', i + 1, number_of_cm_trajectories)

			counter = 0

			for input_label in input_data["labels"] :

				if cm_labels[i] == input_label:

					break

				counter += 1

			sds = np.memmap( sd_temp_filename, dtype = input_data["float_type"],
							shape = ( number_of_cm_trajectories, number_of_timeframes ) )

			sd = sds[i]

			msds[counter] += sd

			del sds

		for _msd, label in zip( msds, input_data["labels"] ):

			_msd /= molecule_numbers[ label ]

	return msds

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

def compute_lengths(input_data, labels, auxiliary_data):

	which_trajectory = 0
	which_length_trajectory = 0
	number_of_length_trajctories = 0

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_beads = auxiliary_data["number_of_beads"]
	number_of_beads_per_file = auxiliary_data["number_of_beads_per_file"]
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	number_of_length_trajectories = auxiliary_data["number_of_molecules"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	traj_temp_filename = auxiliary_data["traj_temp_filename"]

	print( 'number of length trajectories: {}'.format(number_of_length_trajectories) )

	# temporary binary file which will contain the length trajectories
	length_temp_filename = input_data["input_xyz_template"] + 'len_tmp.dat'
	auxiliary_data["length_temp_filename"] = length_temp_filename

	length_trajectories = np.memmap( length_temp_filename, dtype = input_data["float_type"],
									 mode = 'w+',
									 shape = ( number_of_length_trajectories,
									  		   number_of_timeframes ) )

	for i in range(number_of_length_trajectories):
		length_trajectories[i] = np.zeros(number_of_timeframes, dtype = input_data["float_type"])

	del length_trajectories

	length_labels = [ '___' for i in range(number_of_length_trajectories) ]

	trajectories = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
						shape = ( number_of_beads, number_of_timeframes, 3 ) )

	length_trajectories = np.memmap( length_temp_filename, dtype = input_data["float_type"],
									 shape = ( number_of_length_trajectories, number_of_timeframes ) )

	while( which_trajectory < number_of_beads ):

		multiplicity = molecule_sizes[ labels[ which_trajectory ] ]

		if multiplicity == 1:

			length_labels[ which_length_trajectory ] = labels[ which_trajectory ]

			which_trajectory += 1
			which_length_trajectory += 1

		else:

			for i in range( number_of_timeframes ):

				r_ref = trajectories[which_trajectory, i, :]

				for j in range( 1, multiplicity ):

					r = trajectories[which_trajectory + j, i, :]

					_keep_bound_beads_in_the_same_box( r, r_ref, input_data["box_size"] )

					r_ref = r

					trajectories[which_trajectory + j, i, :] = np.array(r, dtype = input_data["float_type"])

			lengths = np.linalg.norm( trajectories[which_trajectory + multiplicity - 1, :, :] - trajectories[which_trajectory, :, :], axis = 1 )
			
			length_trajectories[which_length_trajectory] = lengths

			length_labels[ which_length_trajectory ] = labels[ which_trajectory ]

			which_trajectory += multiplicity
			which_length_trajectory += 1

	del trajectories
	del length_trajectories

	if input_data["verbose"]: print('length computation performed')

	return length_labels

#-------------------------------------------------------------------------------

def compute_length_distribution(input_data, length_labels, auxiliary_data):

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_length_trajectories = auxiliary_data["number_of_molecules"]
	length_temp_filename = auxiliary_data["length_temp_filename"]
	number_of_bins = input_data["number_of_bins"]

	length_trajectories = np.memmap( length_temp_filename, dtype = input_data["float_type"],
						   			 shape = ( number_of_length_trajectories, number_of_timeframes ) )

	lengths_separated = [ [ ] for i in range(len(input_data["sizes"])) ]

	bins = []
	length_distribution = []

	for i in range( number_of_length_trajectories ):

		# if input_data["verbose"]: timestamp('computing sad: {} / {}', i + 1, number_of_orientation_trajectories)

		length_trajectories = np.memmap( length_temp_filename, dtype = input_data["float_type"],
						   				 shape = ( number_of_length_trajectories, number_of_timeframes) )

		length_trajectory = length_trajectories[i]

		del length_trajectories

		counter = 0

		for input_label in input_data["labels"] :

			if length_labels[i] == input_label:

				break

			counter += 1

		lengths_separated[counter].append( length_trajectory )

	for i in range(len(input_data["sizes"])):

		l, b = np.histogram( lengths_separated[i], bins = input_data["number_of_bins"], range = input_data["bin_range"] )

		bins.append( b )

		length_distribution.append( l )

	return bins, length_distribution

#-------------------------------------------------------------------------------

def save_lengths_to_file(input_data, bins, lengths):

	output_filename = input_data["input_xyz_template"] + 'len.txt'

	with open(output_filename, 'w') as output_file:

		first_line = ''
		line = ''

		for i, label in enumerate( input_data["labels"] ):

			first_line += ( 'length/A ' + label + ' ' )

			line += '{} {} '

		output_file.write(first_line + '\n')

		for i in range( input_data["number_of_bins"] ):

			line_values = [ ]

			for j in range( len(input_data["labels"]) ):

				line_values.append( bins[j][i] )
				line_values.append( lengths[j][i] )

			output_file.write( line.format(*line_values) + '\n' )

#-------------------------------------------------------------------------------

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

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes, 3 ) )

	for i in range(number_of_orientation_trajectories):
		orientation_trajectories[i] = np.zeros((number_of_timeframes, 3), dtype = input_data["float_type"])

	del orientation_trajectories

	orientation_labels = [ '___' for i in range(number_of_orientation_trajectories) ]

	trajectories = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
						shape = ( number_of_beads, number_of_timeframes, 3 ) )

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
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

					trajectories[which_trajectory + j, i, :] = np.array(r, dtype = input_data["float_type"])

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

#-------------------------------------------------------------------------------

def compute_mean_squared_angular_displacements(input_data, orientation_labels, auxiliary_data):

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_orientation_trajectories, number_of_timeframes, 3 ) )

	# temporary binary file which will contain the squared angular displacements
	sad_temp_filename = input_data["input_xyz_template"] + 'sad_tmp.dat'
	auxiliary_data["sad_temp_filename"] = sad_temp_filename

	sads = np.memmap( sad_temp_filename, dtype = input_data["float_type"],
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes ) )

	for i in range(number_of_orientation_trajectories):
		sads[i] = np.zeros(number_of_timeframes, dtype = input_data["float_type"])

	msad = [ np.zeros(number_of_timeframes) for i in range(len(input_data["sizes"])) ]

	del orientation_trajectories

	del sads

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing sad: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_orientation_trajectories,
						   			 number_of_timeframes, 3 ) )

		orientation_trajectory = orientation_trajectories[i]

		del orientation_trajectories

		sad = _compute_sad(orientation_trajectory, mode = input_data["mode"])


		sads = np.memmap( sad_temp_filename, dtype = input_data["float_type"],
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

		sads = np.memmap( sad_temp_filename, dtype = input_data["float_type"],
					 	   shape = ( number_of_orientation_trajectories, number_of_timeframes ) )

		sad = sads[i]

		msad[counter] += sad

		del sads

	for _msad, label in zip( msad, input_data["labels"] ):

		_msad /= molecule_numbers[ label ]

	return msad

#-------------------------------------------------------------------------------

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

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	oas = np.memmap( oa_temp_filename, dtype = input_data["float_type"],
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes ) )

	for i in range(number_of_orientation_trajectories):
		oas[i] = np.zeros(number_of_timeframes, dtype = input_data["float_type"])

	moa = [ np.zeros(number_of_timeframes) for i in range(len(input_data["sizes"])) ]
	amounts = np.zeros(len(input_data["sizes"]))

	del orientation_trajectories

	del oas

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('computing autocorrelation: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_orientation_trajectories,
						   	number_of_timeframes, 3 ) )

		orientation_trajectory = orientation_trajectories[i]

		del orientation_trajectories

		autocorrelation = _compute_autocorrelation(orientation_trajectory, mode = input_data["mode"])

		oas = np.memmap( oa_temp_filename, dtype = input_data["float_type"],
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

		oas = np.memmap( oa_temp_filename, dtype = input_data["float_type"],
					 	   shape = ( number_of_orientation_trajectories, number_of_timeframes ) )

		autocorrelation = oas[i]

		moa[counter] += autocorrelation

		del oas

	for _moa, label in zip( moa, input_data["labels"] ):

		_moa /= molecule_numbers[ label ]

	return moa

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

def compute_mean_director( input_data, orientation_labels, auxiliary_data ):

	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	preQ = np.zeros( (number_of_orientation_trajectories, 3, 3), dtype = input_data["float_type"] )

	directors = np.zeros( (number_of_orientation_trajectories, 3), dtype = input_data["float_type"] )

	mean_director = [ np.zeros(3, dtype = input_data["float_type"]) for i in range(len(input_data["sizes"])) ]

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

		mean_director[i] /= np.linalg.norm( mean_director[i] )

	return mean_director

#-------------------------------------------------------------------------------

def compute_nematic_order( input_data, orientation_labels, auxiliary_data, director ):

	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	number_of_orientation_trajectories = auxiliary_data["number_of_molecules"]
	orient_temp_filename = auxiliary_data["orient_temp_filename"]

	mean_second_Legendre_polynomial = np.zeros( number_of_orientation_trajectories, dtype = input_data["float_type"] )
	nematic_order_parameter = np.zeros( len(input_data["sizes"]) )

	orientation_trajectories = np.memmap( orient_temp_filename, dtype = input_data["float_type"],
						   					shape = ( number_of_orientation_trajectories,
						   							  number_of_timeframes, 3 ) )

	amounts = np.zeros(len(input_data["sizes"]))

	for i in range( number_of_orientation_trajectories ):

		if input_data["verbose"]: timestamp('averaging second order Legendre polynomial over trajectory: {} / {}', i + 1, number_of_orientation_trajectories)

		orientation_trajectory = orientation_trajectories[i]

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		for j in range( number_of_timeframes ):

			mean_second_Legendre_polynomial[i] += _second_Legendre_polynomial(orientation_trajectory[j], director[counter])

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

#-------------------------------------------------------------------------------

def compute_rdfs( input_data, labels, auxiliary_data ):

	input_xyz_filenames = auxiliary_data["input_xyz_filenames"]
	number_of_xyz_files = len( input_xyz_filenames )
	number_of_timeframes = auxiliary_data["number_of_timeframes"]
	number_of_trajectories = len(labels)
	number_of_beads_per_file = number_of_trajectories // number_of_xyz_files
	molecule_sizes = auxiliary_data["molecule_sizes"]
	molecule_numbers = auxiliary_data["molecule_numbers"]
	if not CM: traj_temp_filename = auxiliary_data["traj_temp_filename"]
	else: traj_temp_filename = auxiliary_data["cm_temp_filename"]

	trajectories = np.memmap( traj_temp_filename, dtype = input_data["float_type"],
						   shape = ( number_of_trajectories, number_of_timeframes, 3 ) )

	### FREUD VERSION START ###

	if FREUD:

		bin_counts = [ np.zeros( input_data["number_of_bins"], dtype = input_data["float_type"] ) for k in range( len(input_data["labels"]) ) ]

		all_bin_counts = np.zeros( input_data["number_of_bins"], dtype = input_data["float_type"] )

		rdfs = [ np.zeros( input_data["number_of_bins"], dtype = input_data["float_type"] ) for k in range( len(input_data["labels"]) ) ]

		coords_all_boxes = [  ]

		all_coords_all_boxes = [  ]

		dr = input_data["box_size"] / 2 / input_data["number_of_bins"]

		for l in range( number_of_xyz_files ):

			coords_this_box = [  ]

			all_coords_this_box = [ [ trajectories[i][j] for i in range( l * number_of_beads_per_file, (l + 1) * number_of_beads_per_file ) ] for j in range( number_of_timeframes ) ]

			for k in range( len(input_data["labels"]) ):

				coords_list = [ [ trajectories[i][j] for i in range( l * number_of_beads_per_file, (l + 1) * number_of_beads_per_file ) if input_data["labels"][k] == labels[i]  ] for j in range( number_of_timeframes ) ]

				coords_this_box.append( np.array( coords_list, dtype = input_data["float_type"] ) )

			coords_all_boxes.append( coords_this_box )

			all_coords_all_boxes.append( all_coords_this_box )

		box = freud.box.Box.cube( input_data["box_size"] )

		for l in range( number_of_xyz_files ):

			for j in range( number_of_timeframes ):

				coords = np.array( [ all_coords_all_boxes[l][j][i] for i in range( number_of_beads_per_file ) ], dtype = input_data["float_type"] )

				link_cell = freud.locality.LinkCell( box=box, points=coords, cell_width = input_data["box_size"]/2 )

				rdf = freud.density.RDF(bins=input_data["number_of_bins"], r_max=input_data["box_size"] / 2)

				a = rdf.compute(system=link_cell, reset=False)

				all_bin_counts += a.bin_counts

			for k in range( len(input_data["labels"]) ):

				for j in range( number_of_timeframes ):

					if not CM: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] * molecule_sizes[ input_data["labels"][k] ] // number_of_xyz_files

					else: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] // number_of_xyz_files

					coords = np.array( [ coords_all_boxes[l][k][j][i] for i in range( number_of_k_particles_per_box ) ], dtype = input_data["float_type"] )

					link_cell = freud.locality.LinkCell( box=box, points=coords, cell_width = input_data["box_size"]/2 )

					rdf = freud.density.RDF(bins=input_data["number_of_bins"], r_max=input_data["box_size"] / 2)

					a = rdf.compute(system=link_cell, reset=False)

					bin_counts[k] += a.bin_counts

					rdfs[k] += a.bin_counts

		bins = a.bin_centers

		for k in range( len(bin_counts) ):

			if not CM: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] * molecule_sizes[ input_data["labels"][k] ] // number_of_xyz_files

			else: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] // number_of_xyz_files

			bin_counts[k] /= ( number_of_timeframes * number_of_xyz_files )

			rdfs[k] /= ( number_of_timeframes * number_of_xyz_files )

			rdfs[k] /= ( 4.0 * np.pi * bins**2 * dr )

			rdfs[k] /= ( number_of_k_particles_per_box / input_data["box_size"]**3 )

			rdfs[k] /= ( number_of_k_particles_per_box - 1 )

		all_bin_counts /= ( number_of_timeframes * number_of_xyz_files )

		distinct_bin_counts = np.copy( all_bin_counts )

		for k in range( len(bin_counts) ):

			distinct_bin_counts -= bin_counts[k]

		denominator = 0.0

		for k in range( len(bin_counts) ):

			if not CM: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] * molecule_sizes[ input_data["labels"][k] ] // number_of_xyz_files

			else: number_of_k_particles_per_box = molecule_numbers[ input_data["labels"][k] ] // number_of_xyz_files

			for l in range( len(bin_counts) ):

				if not CM: number_of_l_particles_per_box = molecule_numbers[ input_data["labels"][l] ] * molecule_sizes[ input_data["labels"][l] ] // number_of_xyz_files

				else: number_of_l_particles_per_box = molecule_numbers[ input_data["labels"][l] ] // number_of_xyz_files

				if k != l:

					denominator += ( number_of_k_particles_per_box * number_of_l_particles_per_box )

		distinct_rdf = distinct_bin_counts / denominator

		distinct_rdf /= ( 4.0 * np.pi * bins**2 * dr )

		distinct_rdf /= ( 1.0 / input_data["box_size"]**3 )

	### FREUD VERSION END ###

	else:

		timestamp('RDF without Freud is not implemented')

	return bins, rdfs, distinct_rdf

#-------------------------------------------------------------------------------

def save_rdfs_to_file( input_data, bins, rdfs, distinct_rdf ):

	output_filename = input_data["input_xyz_template"] + 'rdf.txt'

	with open(output_filename, 'w') as output_file:

		first_line = 'r/AA '
		line = '{} '

		for label in input_data["labels"]:

			first_line += ( label + ' ' )

			line += '{} '

		first_line += 'XXX'

		line += '{}'

		output_file.write(first_line + '\n')

		for i in range( len(bins) ):

			line_values = [ bins[i] ]

			for j in range( len(input_data["labels"]) ):

				line_values.append( rdfs[j][i] )

			line_values.append( distinct_rdf[i] )

			output_file.write( line.format(*line_values) + '\n' )

#===============================================================================

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
					coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ] )
					trajectories[bead_index, timeframe_index - min_time_index] = coords
			counter += 1

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
                if not ( i == 0 and j == 0 and k == 0 ) ] )

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

# TODO: what to do if n = (0,0,0)
def _compute_autocorrelation(orientation, mode = 'direct'):

	for n in orientation:

		if np.any(n != 0.0): n /= np.linalg.norm(n)

	if mode == 'direct':

		return np.array( [1.0] + [ np.dot( orientation[0], orientation[j] ) for j in range(1, len(orientation)) ] )

	elif mode == 'window':

		return np.array( [1.0] + [ np.mean( [ np.dot(orientation[k], orientation[k+j]) for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ] )

#-------------------------------------------------------------------------------

# TODO: what to do if n = (0,0,0)
def _compute_sad(orientation, mode = 'direct'):

	for n in orientation:

		if np.any(n != 0.0): n /= np.linalg.norm(n)

	if mode == 'direct':

		return np.array( [0.0] + [ np.arccos( np.dot( orientation[0], orientation[j] ) )**2 for j in range(1, len(orientation)) ] )

	elif mode == 'window':

		return np.array( [0.0] + [ np.mean( [ np.arccos( np.dot(orientation[k], orientation[k+j]) )**2 for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ] )

#===============================================================================

def _second_Legendre_polynomial(orientation, director):

	return 1.5 * np.dot( orientation, director )**2 - 0.5

#-------------------------------------------------------------------------------

def _number_of_timeframes( input_xyz_filenames, probing_frequency ):

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  	probing_frequency )
							  	for i in range( number_of_xyz_files ) ]

	return min( numbers_of_timeframes )

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

def _find_min_time_index(input_xyz_filenames, min_time, probing_frequency):

	number_of_xyz_files = len( input_xyz_filenames )

	min_time_index = _min_time_index( input_xyz_filenames[0], min_time, probing_frequency )

	for i in range( 1, number_of_xyz_files ):

		assert _min_time_index( input_xyz_filenames[i], min_time, probing_frequency ) == min_time_index

	return min_time_index

#-------------------------------------------------------------------------------

def _count_beads(filename):

	with open(filename, 'r') as file:
		return( int( file.readline() ) )

#-------------------------------------------------------------------------------

def _number_of_beads( input_xyz_filenames ):

	number_of_beads = _count_beads ( input_xyz_filenames[0] )

	for i in range( 1, len(input_xyz_filenames) ):

		assert _count_beads( input_xyz_filenames[i] ) == number_of_beads

	return number_of_beads * len( input_xyz_filenames )

#-------------------------------------------------------------------------------

def _molecule_sizes( input_labels, input_sizes ):

	molecule_sizes = {}

	for label, size  in zip( input_labels, input_sizes ):
		molecule_sizes[label] = size

	return molecule_sizes

#-------------------------------------------------------------------------------

def _molecule_numbers( input_labels, input_sizes, labels, molecule_sizes ):

	molecule_numbers = {}

	for label, size  in zip( input_labels, input_sizes ):
		molecule_numbers[label] = 0

	for label in labels:
		molecule_numbers[label] += 1

	for label in input_labels:
		molecule_numbers[label] = molecule_numbers[label] // molecule_sizes[label]

	return molecule_numbers

#-------------------------------------------------------------------------------

def _number_of_molecules( molecule_numbers, input_labels ):

	number_of_molecules = 0

	for label in input_labels:
		number_of_molecules += molecule_numbers[label]

	return number_of_molecules

#-------------------------------------------------------------------------------

def _keep_bound_beads_in_the_same_box(r, r_ref, box_size):

	for i in range(0, 3):
		if ( r[i] - r_ref[i] >= box_size / 2 ):
			r[i] -= box_size
		elif ( r_ref[i] - r[i] >= box_size / 2 ):
			r[i] += box_size

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

def unify_coordinates(trajectories, box_size):

	for particle in trajectories:
		for i in range(1, len(particle)):
			_coord_unify(particle[i-1], particle[i], box_size)

#-------------------------------------------------------------------------------

def _coord_unify(past, present, box_size):

	for i in range(3):
		while( present[i] - past[i] >= box_size / 2 ):
			present[i] -= box_size
		while( past[i] - present[i] >= box_size / 2 ):
			present[i] += box_size
