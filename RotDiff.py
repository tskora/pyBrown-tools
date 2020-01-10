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
from pyBrown.trajectories import read_trajectories, unify_coordinates, \
								 _count_timeframes
# from pyBrown.trajectories import read_trajectories, read_energies, \
# 								 separate_center_of_mass, \
# 								 compute_msds, compute_menergies, \
# 								 save_msds_to_file, \
# 								 plot_msds, plot_menergies

#-------------------------------------------------------------------------------

def _compute_angular_displacement(vec0, vec1):

	import numpy as np

	print('hey')

	print( np.linalg.norm(vec0) )
	print( np.linalg.norm(vec1) )
	print( np.dot( vec0, vec1 ) )

	if np.linalg.norm(vec0) == 0.0 or np.linalg.norm(vec1) == 0.0:

		return 0.0

	return np.arccos( np.dot( vec0, vec1 ) / np.linalg.norm(vec0) / np.linalg.norm(vec1) )

#-------------------------------------------------------------------------------

def _compute_angle(vec):

	import numpy as np

	return _compute_angular_displacement( vec, np.array([1.0, 0.0, 0.0], float) )

#-------------------------------------------------------------------------------

def _compute_msa(angles, mode):

	import numpy as np

	import freud.box
	import freud.msd

	angles_transformed = [ [ np.array( [ angles[i][j], 0.0, 0.0 ], float ) for i in range( len( angles ) ) ] for j in range( len( angles[0] ) ) ]

	if ( mode == "direct" ):
		msa = freud.msd.MSD( mode = 'direct' )
	elif ( mode == "window" ):
		msa = freud.msd.MSD( mode = 'window' )

	msa.compute( positions = angles_transformed )
	
	return msa.msd

#-------------------------------------------------------------------------------

def _keep_bound_beads_in_the_same_box(r, r_ref, box_size):

	for i in range(0, 3):
		if ( r[i] - r_ref[i] >= box_size / 2 ):
			r[i] -= box_size
		elif ( r_ref[i] - r[i] >= box_size / 2 ):
			r[i] += box_size

#-------------------------------------------------------------------------------

def compute_orientations(input_data, temporary_filename, labels):

	import numpy as np

	import os

	bead_dict = {}
	mol_dict = {}

	which_trajectory = 0
	which_orientation_trajectory = 0

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  input_data["probing_frequency"] )
							  for i in range( number_of_xyz_files ) ]

	number_of_timeframes = min( numbers_of_timeframes )

	number_of_orientation_trajectories = 0

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

		mol_dict[label] = 0

	for label in labels:

		mol_dict[label] += 1

	for label in input_data["labels"]:

		mol_dict[label] = mol_dict[label] // bead_dict[label]

		number_of_orientation_trajectories += mol_dict[label]

	print( 'number of orientation trajectories: {}'.format(number_of_orientation_trajectories) )

	temporary_filename_2 = 'temp2.dat'

	orientation_trajectories = np.memmap( temporary_filename_2, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_orientation_trajectories,
									  number_of_timeframes, 3 ) )

	for i in range(number_of_orientation_trajectories):
		orientation_trajectories[i, :] = np.zeros((number_of_timeframes, 3), dtype = np.float32)

	del orientation_trajectories

	orientation_labels = [ '___' for i in range(number_of_orientation_trajectories) ]

	while( which_trajectory < len( labels ) ):

		multiplicity = bead_dict[ labels[ which_trajectory ] ]

		temp = np.memmap( temporary_filename, dtype = np.float32,
							shape = ( len(labels), number_of_timeframes, 3 ) )

		temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
							shape = ( len(orientation_labels), number_of_timeframes, 3 ) )

		if multiplicity == 1:

			temp2[which_orientation_trajectory, :, :] = temp[which_trajectory, :, :]
			orientation_labels[ which_orientation_trajectory ] = labels[ which_trajectory ]

			which_trajectory += 1
			which_orientation_trajectory += 1

		else:

			for i in range( number_of_timeframes ):

				r_ref = temp[which_trajectory, i, :]

				for j in range( 1, multiplicity ):

					r = temp[which_trajectory + j, i, :]

					_keep_bound_beads_in_the_same_box( r, r_ref, input_data["box_size"] )

					r_ref = r

					temp[which_trajectory + j, i, :] = np.array(r, dtype = np.float32)

			orientations = temp[which_trajectory + multiplicity - 1, :, :] - temp[which_trajectory, :, :]

			orientation_norms = np.linalg.norm(orientations, axis = 1)

			for v, vn in zip(orientations, orientation_norms):
				v /= vn
			
			temp2[which_orientation_trajectory, :] = orientations

			orientation_labels[ which_orientation_trajectory ] = labels[ which_trajectory ]

			which_trajectory += multiplicity
			which_orientation_trajectory += 1

		del temp
		del temp2

	if input_data["verbose"]: print('orientation computation performed')

	os.remove(temporary_filename)

	return temporary_filename_2, orientation_labels

#-------------------------------------------------------------------------------

def _compute_autocorrelation(orientation, mode = 'direct'):

	import numpy as np

	if mode == 'direct':

		return np.array( [ np.dot( orientation[0], orientation[j] ) for j in range(len(orientation)) ], np.float32 )

	elif mode == 'window':

		return np.array( [1.0] + [ np.mean( [ np.dot(orientation[k], orientation[k+j]) for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ], np.float32 )

#-------------------------------------------------------------------------------

def _compute_sad(orientation, mode = 'direct'):

	import numpy as np

	if mode == 'direct':

		return np.array( [ np.arccos( np.dot( orientation[0], orientation[j] ) )**2 for j in range(len(orientation)) ], np.float32 )

	elif mode == 'window':

		return np.array( [0.0] + [ np.mean( [ np.arccos( np.dot(orientation[k], orientation[k+j]) )**2 for k in range(len(orientation)-j) ] ) for j in range(1, len(orientation)) ], np.float32 )

#-------------------------------------------------------------------------------

def compute_average_orientation_autocorrelation(input_data, times, temporary_filename_2, orientation_labels):

	import numpy as np

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	number_of_timeframes = len(times)

	bead_dict = {}

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	temporary_filename_2 = 'temp2.dat'

	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(orientation_labels), len(times), 3 ) )

	temporary_filename_3 = 'temp3.dat'

	temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
							mode = 'w+',
							shape = ( len(orientation_labels),
									  len(times) ) )

	for i in range(len(orientation_labels)):
		temp3[i] = np.zeros(len(times), dtype = np.float32)

	average_autocorrelation = [ np.zeros(len(times)) for i in range(len(input_data["sizes"])) ]
	amounts = np.zeros(len(input_data["sizes"]))

	del temp2

	del temp3

	for i in range( len(orientation_labels) ):

		if input_data["verbose"]: timestamp('computing autocorrelation: {} / {}', i, len(orientation_labels))

		temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(orientation_labels), number_of_timeframes, 3 ) )

		orientation_trajectory = temp2[i]

		del temp2

		autocorrelation = _compute_autocorrelation(orientation_trajectory, mode = 'window')


		temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
					 	   shape = ( len(orientation_labels), number_of_timeframes ) )

		temp3[i] = autocorrelation

		del temp3

	for i in range( len(orientation_labels) ):

		if input_data["verbose"]: timestamp('averaging: {} / {}', i, len(orientation_labels))

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
					 	   shape = ( len(orientation_labels), number_of_timeframes ) )

		autocorrelation = temp3[i]

		average_autocorrelation[counter] += autocorrelation

		del temp3

		amounts[counter] += 1

	for ava, amount in zip( average_autocorrelation, amounts ):

		ava /= amount

	return average_autocorrelation

#-------------------------------------------------------------------------------

def save_aoas_to_file(input_data, times, aoas):

	output_filename = input_data["input_xyz_template"] + 'aoa.txt'

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

				line_values.append( aoas[j][i] )

			output_file.write( line.format(*line_values) + '\n' )

#-------------------------------------------------------------------------------

def save_msads_to_file(input_data, times, msads):

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

def compute_msads(input_data, times, temporary_filename_2, orientation_labels):

	import numpy as np

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	number_of_timeframes = len(times)

	bead_dict = {}

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	temporary_filename_2 = 'temp2.dat'

	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(orientation_labels), len(times), 3 ) )

	temporary_filename_4 = 'temp4.dat'

	temp4 = np.memmap( temporary_filename_4, dtype = np.float32,
							mode = 'w+',
							shape = ( len(orientation_labels),
									  len(times) ) )

	for i in range(len(orientation_labels)):
		temp4[i] = np.zeros(len(times), dtype = np.float32)

	msads = [ np.zeros(len(times)) for i in range(len(input_data["sizes"])) ]
	amounts = np.zeros(len(input_data["sizes"]))

	del temp2

	del temp4

	for i in range( len(orientation_labels) ):

		if input_data["verbose"]: timestamp('computing sad: {} / {}', i, len(orientation_labels))

		temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(orientation_labels), number_of_timeframes, 3 ) )

		orientation_trajectory = temp2[i]

		del temp2

		sad = _compute_sad(orientation_trajectory, mode = 'window')


		temp4 = np.memmap( temporary_filename_4, dtype = np.float32,
					 	   shape = ( len(orientation_labels), number_of_timeframes ) )

		temp4[i] = sad

		del temp4

	for i in range( len(orientation_labels) ):

		if input_data["verbose"]: timestamp('averaging: {} / {}', i, len(orientation_labels))

		counter = 0

		for input_label in input_data["labels"] :

			if orientation_labels[i] == input_label:

				break

			counter += 1

		temp4 = np.memmap( temporary_filename_4, dtype = np.float32,
					 	   shape = ( len(orientation_labels), number_of_timeframes ) )

		sad = temp4[i]

		msads[counter] += sad

		del temp4

		amounts[counter] += 1

	for msad, amount in zip( msads, amounts ):

		msad /= amount

	return msads

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
				"probing_frequency": 1, "min_time": 0.0, "mode": "window"}

	required_keywords = []
	defaults = {}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataRotDiff(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	temporary_filename, times, labels, min_time_index = read_trajectories(i.input_data)
	timestamp( 'Computing the angles' )
	temporary_filename_2, orientation_labels = compute_orientations(i.input_data, temporary_filename, labels)
	del labels

	import matplotlib.pyplot as plt
	import numpy as np

	aoas = compute_average_orientation_autocorrelation(i.input_data, times, temporary_filename_2, orientation_labels)
	timestamp( 'Saving average orientation autocorrelations to a file' )
	save_aoas_to_file(i.input_data, times, aoas)

	msads = compute_msads(i.input_data, times, temporary_filename_2, orientation_labels)
	timestamp( 'Saving mean squared angular displacements to a file' )
	save_msads_to_file(i.input_data, times, msads)

	ys = aoas[1]
	a, b = np.polyfit(times[:len(times)//4], ys[:len(times)//4], 1)
	plt.plot( times, ys, '-')
	plt.plot(times, a * times + b * np.ones(len(times)), ':')
	plt.show()
	plt.close()
	print(-0.5*a)

	zs = msads[1]
	aa, bb = np.polyfit( times[:len(times)//4], zs[:len(times)//4], 1 )
	plt.plot( times, zs, '-' )
	plt.plot( times, aa * times + bb * np.ones(len(times)), ':' )
	plt.show()
	plt.close()
	print(0.25*aa)

	1/0

	timestamp( 'Computing mean square angle displacements' )
	msas = compute_msas(i.input_data, temporary_filename_2, angle_labels, min_time_index)
	del angle_labels
	# timestamp( 'Saving mean square displacements to a file' )
	# save_msds_to_file(i.input_data, times, msds)
	# timestamp( 'Plotting mean square displacements' )
	# plot_msds(i.input_data, times, msds)
	# del times
	# del msds

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
