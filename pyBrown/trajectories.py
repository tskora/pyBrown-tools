# pyBrown is a bundle of tools useful for Brownian and Stokesian dynamics simulations
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

import os
import numpy as np
import matplotlib.pyplot as plt
from pyBrown.plot_config import plot_config
from scipy.constants import Boltzmann

#-------------------------------------------------------------------------------

def read_trajectories(input_data_object):

	input_data = input_data_object.input_data

	if input_data["debug"]:
		print( 'input data from JSON file: {}'.format(input_data) )

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	if input_data["debug"]:
		print( 'input xyz filenames: {}'.format(input_xyz_filenames) )

	number_of_timeframes = _count_timeframes( input_xyz_filenames[0],
											input_data["probing_frequency"] )

	number_of_beads = _count_beads( input_xyz_filenames[0] )

	if input_data["verbose"] or input_data["debug"]:
		print( 'timeframes: {}\nbeads: {}\nfiles: {}'.format( number_of_timeframes,
															number_of_beads,
															len( input_xyz_filenames ) ) )

	# temporary binary file which will contain the trajectories
	temporary_filename = 'temp.dat'

	trajectories = np.memmap( temporary_filename, dtype = np.float32,
							mode = 'w+',
							shape = ( number_of_beads * len(input_xyz_filenames),
									  number_of_timeframes, 3 ) )

	for i in range(number_of_beads * len(input_xyz_filenames)):
		for j in range(number_of_timeframes):
			trajectories[i, j, :] = np.zeros(3, np.float32)

	del trajectories

	times = np.zeros( number_of_timeframes, dtype = np.float32 )
	labels = [ '___' for _ in range( len(input_xyz_filenames) * number_of_beads ) ]

	for i, input_xyz_filename in enumerate(input_xyz_filenames):

		if input_data["verbose"] or input_data["debug"]:
			print( 'xyz file {}/{}'.format( i + 1, len(input_xyz_filenames) ) )

		with open(input_xyz_filename, 'r') as input_xyz_file:

			temp = np.memmap( temporary_filename, dtype = np.float32,
						   shape = ( number_of_beads * len(input_xyz_filenames),
						   number_of_timeframes, 3 ) )

			trajectories = temp[i * number_of_beads : (i + 1) * number_of_beads, :, :]

			_read_trajectories_from_xyz_file(input_xyz_file, trajectories, times,
											 labels, input_data["probing_frequency"])

			del temp

	if input_data["verbose"] or input_data["debug"]: print('xyzs read')

	return temporary_filename, times, labels

#-------------------------------------------------------------------------------

def read_energies(input_data_object):

	input_data = input_data_object.input_data

	input_enr_filenames = [ input_data["input_enr_template"] + str(i) + '.enr'
							for i in range( *input_data["input_enr_range"] ) ]

	energies = []
	times = []

	for i, input_enr_filename in enumerate(input_enr_filenames):

		if input_data["verbose"]: print('enr file {}/{}'.format( i + 1, len(input_enr_filenames) ))

		with open(input_enr_filename, 'r') as input_enr_file:

			_read_energies_from_enr_file(input_enr_file, energies, times)

	if input_data["verbose"]: print('enrs read')

	return np.array(energies), np.array(times, float)

#-------------------------------------------------------------------------------

def separate_center_of_mass(input_data_object, temporary_filename, labels):

	input_data = input_data_object.input_data
	bead_dict = {}

	which_trajectory = 0
	which_cm_trajectory = 0

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_timeframes = _count_timeframes( input_xyz_filenames[0], input_data["probing_frequency"] )

	number_of_cm_trajectories = 0

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	for label in labels: # thad division not works perfect with integers

		number_of_cm_trajectories += 1 / bead_dict[label]

	number_of_cm_trajectories = int( number_of_cm_trajectories )

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

			temp2[which_cm_trajectory, :, :] = temp[which_trajectory, :, :]

			for i in range( 1, multiplicity ):

				temp2[which_cm_trajectory, :, :] += temp[which_trajectory + i, :, :]

			temp2[which_cm_trajectory, :, :] /= multiplicity

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

			for i in range( 1, multiplicity ):

				temp2[which_cm_trajectory, :, :] += temp[which_trajectory + i, :, :]

			temp2[which_cm_trajectory, :, :] /= multiplicity

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

	if input_data["verbose"]: print('cm separation performed')

	os.remove(temporary_filename)

	return temporary_filename_2, cm_labels

#-------------------------------------------------------------------------------

def compute_msds(input_data_object, temporary_filename_2, cm_labels):

	input_data = input_data_object.input_data

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_timeframes = _count_timeframes( input_xyz_filenames[0], input_data["probing_frequency"] )

	bead_dict = {}

	temporary_filename_3 = 'temp3.dat'

	sds = np.memmap( temporary_filename_3, dtype = np.float32,
					 mode = 'w+',
					 shape = ( len(cm_labels), number_of_timeframes ) )

	if input_data["verbose"]: print("filling sds")

	for i in range( len(cm_labels) ):
		sds[i, :] = np.zeros( number_of_timeframes, dtype = np.float32 ) 

	del sds

	if input_data["verbose"]: print("end filling sds")

	msds = np.zeros( ( len( input_data["labels"] ), number_of_timeframes ), float )
	amounts = np.zeros( len( input_data["labels"] ), float )

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	for i in range( len(cm_labels) ):

		if input_data["verbose"]: print('computing sd: {} / {}'.format(i, len(cm_labels)))

		temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(cm_labels), number_of_timeframes, 3 ) )

		cm_trajectory = temp2[i]

		del temp2

		sd = _compute_sd(cm_trajectory, input_data["box_size"])

		temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
					 	   shape = ( len(cm_labels), number_of_timeframes ) )

		temp3[i] = sd

		del temp3

	for i in range( len(cm_labels) ):

		if input_data["verbose"]: print('averaging: {} / {}'.format(i, len(cm_labels)))

		counter = 0

		for input_label in input_data["labels"] :

			if cm_labels[i] == input_label:

				break

			counter += 1

		temp3 = np.memmap( temporary_filename_3, dtype = np.float32,
					 	   shape = ( len(cm_labels), number_of_timeframes ) )

		sds = temp3[i]

		msds[counter] += sds

		del temp3

		amounts[counter] += 1

	for msd, amount in zip( msds, amounts ):

		msd /= amount

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

def plot_msds(input_data_object, times, msds):

	input_data = input_data_object.input_data

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

def plot_menergies(input_data_object, times, menergies):

	input_data = input_data_object.input_data

	plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'E [$\frac{kcal}{mol}$]')

	plt.plot(times / 1000000, menergies, '-')

	plt.savefig(input_data["input_enr_template"] + 'enr.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def _read_trajectories_from_xyz_file(xyz_file, trajectories, times, labels, probing_frequency):

	number_of_beads = int( xyz_file.readline() )

	counter = 0

	label_index = 0
	for label in labels:
		if label != '___':
			label_index += 1
		else:
			break

	for i, line in enumerate(xyz_file):

		if 'xyz' in line.split()[0].split('.'):

			# print('counter: {}; counter // number_of_beads // probing_frequency: {}; len(times): {}'.format(counter, counter // number_of_beads // probing_frequency, len(times)))

			if ( ( counter // number_of_beads ) % probing_frequency ) == 0:

				times[ counter // number_of_beads // probing_frequency ] = float( line.split()[3] )

				# print(times)


		elif len( line.split() ) == 1: continue

		else:

			if counter < number_of_beads:
				labels[label_index + counter] = line.split()[0]
			if ( ( counter // number_of_beads ) % probing_frequency ) == 0:
				coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ], float )
				trajectories[counter % number_of_beads, counter // number_of_beads // probing_frequency] = coords
			counter += 1

#-------------------------------------------------------------------------------

def _read_energies_from_enr_file(enr_file, energies, times):

	energies_from_current_file = []

	if len(times) == 0: add_times = True
	else: add_times = False

	enr_file.readline()

	for i, line in enumerate(enr_file):

		if add_times: times.append( float( line.split()[0] ) )

		else: assert( float( line.split()[0] ) == times[i] )

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

	if frequency == 1:
		return _file_length(filename) // (number_of_beads + 2)
	else:
		return _file_length(filename) // (number_of_beads + 2) // frequency# + 1