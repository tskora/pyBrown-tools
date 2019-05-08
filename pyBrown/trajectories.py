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

	for input_xyz_filename in input_xyz_filenames:

		with open(input_xyz_filename, 'r') as input_xyz_file:

			_read_trajectories_from_xyz_file(input_xyz_file, trajectories, times, labels)

	return trajectories, times, labels

#-------------------------------------------------------------------------------

def read_energies(input_data_object):

	input_data = input_data_object.input_data

	input_enr_filenames = [ input_data["input_enr_template"] + str(i) + '.enr'
							for i in range( *input_data["input_enr_range"] ) ]

	energies = []
	times = []

	for input_enr_filename in input_enr_filenames:

		with open(input_enr_filename, 'r') as input_enr_file:

			_read_energies_from_enr_file(input_enr_file, energies, times)

	return energies, times

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

	return cm_trajectories, cm_labels

#-------------------------------------------------------------------------------

def compute_msds_and_mjumps(input_data_object, cm_trajectories, cm_labels):

	input_data = input_data_object.input_data

	bead_dict = {}
	sds = []
	jumps_ensemble = []
	msds = [ np.zeros(len(cm_trajectories[0]), float) for i in range( len( input_data["labels"] ) ) ]
	mjumps = [ np.zeros(len(cm_trajectories[0]), float) for i in range( len( input_data["labels"] ) ) ]
	amounts = [ 0 for i in range( len( input_data["labels"] ) ) ]

	for label, size  in zip( input_data["labels"], input_data["sizes"] ):

		bead_dict[label] = size

	for cm_trajectory in cm_trajectories:

		sd, jumps = _compute_sd_and_jumps(cm_trajectory, input_data["box_size"])

		sds.append( sd )
		jumps_ensemble.append( jumps )


	for sd, jumps, cm_label in zip( sds, jumps_ensemble, cm_labels ):

		counter = 0

		for input_label in input_data["labels"] :

			if cm_label == input_label:

				break

			counter += 1

		msds[counter] += sd

		mjumps[counter] += jumps

		amounts[counter] += 1

	for msd, mjump, amount in zip( msds, mjumps, amounts ):

		msd /= amount

		mjump /= amount

	return msds, mjumps

#-------------------------------------------------------------------------------

def compute_menergies(energies):

	menergies = np.zeros( len( energies[0] ), float )

	for energies_one_file in energies:

		menergies += energies_one_file

	return menergies / len( energies )

#-------------------------------------------------------------------------------

def plot_msds(input_data_object, times, msds):

	input_data = input_data_object.input_data

	plot_config()

	plt.xlabel('time [ps]')
	plt.ylabel(r'MSD [$\AA ^2$]')

	plt.yticks([])

	for i, msd in enumerate(msds):

		plt.plot(times, msd, '-', label = input_data["labels"][i] )

		a, b = np.polyfit(times, msd, 1)
		D = a / 6
		D_string = '{:7.5f}'.format(D)
		Rh = 10**17 * Boltzmann / (6 * np.pi) * input_data["temperature"] / ( D * 0.1 * input_data["viscosity"] )
		Rh_string = '{:4.2f}'.format(Rh)

		plt.text( 0.0 * times[-1], (0.60 - i * 0.15) * np.ndarray.max( np.array( msds ) ), input_data["labels"][i] )
		plt.text( 0.0 * times[-1], (0.55 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$D =$' + D_string + r' $\frac{\AA ^2}{ps}$')
		plt.text( 0.0 * times[-1], (0.50 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$R_H =$' + Rh_string + ' nm')

		plt.plot(times, a * np.array(times, float) + b * np.ones(len(times)), '--', label = 'linear fit for ' + input_data["labels"][i]  )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'msd.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def plot_mjumps(input_data_object, times, mjumps):

	input_data = input_data_object.input_data

	plot_config()

	plt.xlabel('time [ps]')
	plt.ylabel(r'$\Delta$ D [$\AA$]')

	for i, mjump in enumerate(mjumps):

		plt.plot( times[1:], mjump[1:], '-', label = input_data["labels"][i] )

		plt.plot( times[1:], np.mean(mjump[1:]) * np.ones( len(mjump[1:]) ), '--', label = 'mean ' + input_data["labels"][i] )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'stp.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def plot_menergies(input_data_object, times, menergies):

	input_data = input_data_object.input_data

	plot_config()

	plt.xlabel('time [ps]')
	plt.ylabel(r'E [$\frac{kcal}{mol}$]')

	plt.plot(times, menergies, '-')

	plt.savefig(input_data["input_enr_template"] + 'enr.jpg', dpi = 100)

	plt.close()

#-------------------------------------------------------------------------------

def _read_trajectories_from_xyz_file(xyz_file, trajectories, times, labels):

	trajectories_from_current_file = []
	counter = 0

	number_of_beads = int( xyz_file.readline() )

	if len(times) == 0: add_times = True
	else: add_times = False

	for i in range( number_of_beads ):

		trajectories_from_current_file.append( [] )

	for line in xyz_file:

		if 'xyz' in line.split()[0].split('.'):

			if add_times: times.append( float( line.split()[3] ) )

		elif len( line.split() ) == 1: continue

		else:

			if counter < number_of_beads: labels.append( line.split()[0] )
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

	for line in enr_file:

		if add_times: times.append( float( line.split()[0] ) )

		energies_from_current_file.append( float( line.split()[1] ) )

	energies.append( np.array( energies_from_current_file, float ) )

#-------------------------------------------------------------------------------

def _compute_sd_and_jumps(trajectory, box_size):

    sd = [0.0]
    jumps = [0.0]
    
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
        
        sd.append( np.sum( ( trajectory[i] - trajectory[0] )**2 ) )
        jumps.append( jump )
                          
    return sd, jumps