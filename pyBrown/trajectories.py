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

import numpy as np
import matplotlib.pyplot as plt
from pyBrown.plot_config import plot_config
from scipy.constants import Boltzmann

#-------------------------------------------------------------------------------

def read_trajectories(input_data_object):

	input_data = input_data_object.input_data

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	trajectories = []
	times = []
	labels = []

	for i, input_xyz_filename in enumerate(input_xyz_filenames):

		print('xyz file {}/{}'.format( i + 1, len(input_xyz_filenames) ))

		with open(input_xyz_filename, 'r') as input_xyz_file:

			_read_trajectories_from_xyz_file(input_xyz_file, trajectories, times, labels, input_data["probing_frequency"])

	print('xyzs read')

	return np.array(trajectories, float), np.array(times, float), labels

#-------------------------------------------------------------------------------

def read_energies(input_data_object):

	input_data = input_data_object.input_data

	input_enr_filenames = [ input_data["input_enr_template"] + str(i) + '.enr'
							for i in range( *input_data["input_enr_range"] ) ]

	energies = []
	times = []

	for i, input_enr_filename in enumerate(input_enr_filenames):

		print('enr file {}/{}'.format( i + 1, len(input_enr_filenames) ))

		with open(input_enr_filename, 'r') as input_enr_file:

			_read_energies_from_enr_file(input_enr_file, energies, times)

	print('enrs read')

	return np.array(energies, float), np.array(times, float)

#-------------------------------------------------------------------------------

def separate_center_of_mass(input_data_object, trajectories, labels):

	input_data = input_data_object.input_data
	bead_dict = {}

	counter = 0
	cm_trajectories = []
	cm_labels = []

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	while( counter < len( labels ) ):

		multiplicity = bead_dict[ labels[ counter ] ]

		if multiplicity == 1:

			cm_trajectories.append( trajectories[ counter ] )
			cm_labels.append( labels[ counter ] )

			counter += 1

		else:

			cm_trajectory = trajectories[ counter ]

			for i in range( counter + 1, counter + multiplicity ):

				cm_trajectory += trajectories[i]

			cm_trajectory /= multiplicity

			cm_trajectories.append( cm_trajectory )
			cm_labels.append( labels[ counter ] )

			counter += multiplicity

	print('cm separation performed')

	return cm_trajectories, cm_labels

#-------------------------------------------------------------------------------

def compute_msds(input_data_object, cm_trajectories, cm_labels):

	input_data = input_data_object.input_data

	bead_dict = {}
	sds = []
	msds = [ np.zeros(len(cm_trajectories[0]), float) for i in range( len( input_data["labels"] ) ) ]
	amounts = [ 0 for i in range( len( input_data["labels"] ) ) ]

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	for cm_trajectory in cm_trajectories:

		sd = _compute_sd(cm_trajectory, input_data["box_size"])

		sds.append( sd )


	for sd, cm_label in zip( sds, cm_labels ):

		counter = 0

		for input_label in input_data["labels"] :

			if cm_label == input_label:

				break

			counter += 1

		msds[counter] += sd

		amounts[counter] += 1

	for msd, amount in zip( msds, amounts ):

		msd /= amount

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
			D = a / 6 # in A**2 / ps
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

	trajectories_from_current_file = []
	counter = 0

	number_of_beads = int( xyz_file.readline() )

	if len(times) == 0: add_times = True
	else: add_times = False

	for i in range( number_of_beads ):

		trajectories_from_current_file.append( [] )

	for i, line in enumerate(xyz_file):

		if (i % 1000000 == 0): print('line: {}'.format(i))

		if 'xyz' in line.split()[0].split('.'):

			if ( ( counter // number_of_beads ) % probing_frequency ) == 0:

				if add_times: times.append( float( line.split()[3] ) )

				else:
					assert( float( line.split()[3] ) == times[ ( counter // number_of_beads ) // probing_frequency ] ) 

		elif len( line.split() ) == 1: continue

		else:

			if counter < number_of_beads: labels.append( line.split()[0] )
			if ( ( counter // number_of_beads ) % probing_frequency ) == 0:
				coords = [ float(line.split()[x]) for x in range(1, 4) ]
				trajectories_from_current_file[counter % number_of_beads].append(coords)
			counter += 1

	for trajectory in trajectories_from_current_file:

		trajectories.append( np.array( trajectory, float ) )

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
        
        if jump >= box_size / 2:
            versors = box_size * np.array( [ [i, j, k]
            	for i in [-1,0,1] for j in [-1,0,1] for k in [-1,0,1]
                            if not (i == 0 and j == 0 and k == 0 ) ], float )

            for versor in versors:
                replica = versor + trajectory[i]
                replica_jump = np.sqrt( np.sum( ( replica - trajectory[i - 1] )**2 ) )
                if replica_jump <= jump:
                    versor_chosen = versor
                    jump = replica_jump

            for j in range(i, len(trajectory)):
                trajectory[j] += versor_chosen
        
        sd[i] = np.sum( ( trajectory[i] - trajectory[0] )**2 )
                          
    return sd

#-------------------------------------------------------------------------------

def file_length(file):

	counter = 0

	for line in file:
		counter += 1

	return counter