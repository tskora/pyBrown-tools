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

#-------------------------------------------------------------------------------

# WARNING: NOT EVERYTHING WORKS SMOOTHLY

def read_trajectory(filename):

	time = []
	trajectory = []

	with open(filename, 'r') as file:

		counter = 0
		number_of_beads = int( file.readline() )

		for i in range(number_of_beads):

			trajectory.append([])

		for line in file:

			if 'xyz' in line.split()[0].split('.'):
				time.append( float( line.split()[3] ) )

			#elif str(number_of_beads) in line.split():
			elif len( line.split() ) == 1:
				continue

			else:
				coords = [ float(line.split()[x]) for x in range(1, 4) ]
				trajectory[counter % number_of_beads].append(coords)
				counter += 1

	return time, trajectory

#-------------------------------------------------------------------------------

def msd(trajectory):
    
	msds = []
        
	for bead in trajectory:

		bead_msd = []

		for trajectory_point in bead:

			_msd = ( ( np.array( trajectory_point ) - np.array( bead[0] ) )**2 ).sum()
			bead_msd.append( _msd )

		msds.append(bead_msd)
    
	return msds

#-------------------------------------------------------------------------------

def msd_cm(trajectory):
    
	trajectory_cm = []
    
	for i in range( len( trajectory[0] ) ):

		cm = np.array([0.0, 0.0, 0.0], float)

		for j in range(len(trajectory)):

			cm += np.array(trajectory[j][i]) / len(trajectory)
        
		trajectory_cm.append(cm)
            
	return msd( [trajectory_cm] )

#-------------------------------------------------------------------------------

def save_msd_to_file(time, m, filename):

	with open(filename, 'w') as output_file:

		for i in range(len(time)):

			m_values = ""

			for bead in m:

				m_values += " " + str( bead[i] )

			output_file.write( str( time[i] ) + m_values + '\n' )





			# output_file.write( str( time[i] ) + ' ' + str( m[0][i] ) + '\n' )

#-------------------------------------------------------------------------------

def read_msd_from_file(filename):

	time = []
	m = [  ]

	with open(filename, 'r') as input_file:

		for n, line in enumerate(input_file):

			if n == 0:

				for i in range( 1, len( line.split() ) ):

					m.append( [] )

			time.append(float( line.split()[0] ) )

			for i in range( 1, len( line.split() ) ):

				m[i].append( float( line.split()[i] ) )

			# m[0].append(float(line.split()[1]))

	return time, m

	# I HAVE REFACTORIZED UP TO THIS POINT BUT NOT YET TESTED

#-------------------------------------------------------------------------------

def plot_msd(time, m, filename):

	plt.title('Diffusion of DNA in dilute regime')
	plt.xlabel('time [ps]')
	plt.ylabel(r'MSD [$\AA ^2$]')

	plt.plot(time, m[0], '-', label = 'mean', lw = 1)

	a, b = np.polyfit(time, m[0], 1)
	diff = '{:7.6f}'.format(a / 6)
	rh = '{:7.6f}'.format(6 * 0.0214 / a)

	# plt.text(len(time) / 2, 2000, r'$D =$' + diff + r' $\frac{\AA ^2}{ps}$')
	# plt.text(len(time) / 2, 20, r'$R_H =$' + rh + ' nm')
	plt.text(0.5 * len(time), 0.2 * np.ndarray.max( m ), r'$D =$' + diff + r' $\frac{\AA ^2}{ps}$')
	plt.text(0.5 * len(time), 0.1 * np.ndarray.max( m ), r'$R_H =$' + rh + ' nm')

	plt.plot(time, a * np.array(time, float) + b * np.ones(len(time)), '--', label = 'linear fit', lw = 1)

	plt.legend()

	plt.savefig(filename, dpi = 300)

	plt.close()

#-------------------------------------------------------------------------------

def plot_d(time, m, filename):

	plt.title('Diffusion of DNA in dilute regime')
	plt.xlabel('time [ps]')
	plt.ylabel(r'D [$\frac{\AA ^2}{ps}$]')

	plt.plot(time[1:], m[0][1:] / np.array(time[1:]) / 6, '-', label = 'mean', lw = 1)

	plt.savefig(filename, dpi = 300)

	plt.close()

#-------------------------------------------------------------------------------

def prepare_time_and_msd_cm(start, end, steps, filename_template):

	### SOMETHING DOES NOT WORK HERE WHEN STEPS < STEPS IN FILE

	number_of_samples = end - start + 1
	
	m = [ np.zeros(steps, float) ]
	counter = 0
	
	for i in range(start, end + 1):
		time, trajectory = read_trajectory(filename_template + str(i) + '.xyz')
		m += np.array( msd_cm(trajectory), float ) / number_of_samples
		counter += 1
		print(str(counter) + '/' + str(number_of_samples))
	
	print('---')

	return time, m

#-------------------------------------------------------------------------------

########## DNA BONDS

def elongation_of_dna_wrapped():

	original_length = 22.8
	filename = 'fix_occ_05_comp.xyz'
	time, trajectory = read_trajectory(filename)
	av_bonds = []
	lengths = []
	for i, t in enumerate(time):

		av_bond = 0.0
		length = original_length

		for j in range( len(trajectory) - 1 ):

			diff = np.array( trajectory[j][i], float ) - \
				   np.array( trajectory[j + 1][i], float )

			av_bond += np.sqrt( np.transpose( diff ) @ diff ) / ( len(trajectory) - 1 )
			length += np.sqrt( np.transpose( diff ) @ diff )

		av_bonds.append(av_bond)
		lengths.append(length)

	plt.title('Bond length of DNA in dilute regime')
	plt.xlabel('time [ps]')
	plt.ylabel(r'R [$\AA$]')
	plt.plot(time, av_bonds, '--', color = 'red')
	plt.hlines(22.8, time[0], time[-1], color = 'black')
	plt.savefig('bond_0_986589.jpg', dpi = 300)
	# plt.xlim((0.0, 1000.0))
	# plt.show()
	plt.close()

	plt.title('DNA strand length elongation in dilute regime')
	plt.xlabel('time [ps]')
	plt.ylabel(r'R [$\AA$]')
	plt.plot(time, np.array(lengths) - np.ones( len(lengths) ) * 8 * original_length, '-', color = 'red')
	plt.savefig('strand_elongation_0_986589.jpg', dpi = 300)
	plt.close()

#-------------------------------------------------------------------------------

########## FOR DNA HYDRODYNAMIC RADIUS

def hydrodynamic_radius_of_dna_wrapped():

	start = 1
	end = 25
	steps = 500000
	
	filename_template = '/Users/tomaszskora/Documents/plgtskora/long_check/dna_986589_'
	filename_template_d = '/Users/tomaszskora/Documents/plgtskora/long_check/dna_986589_d_av'
	filename_template_msd = '/Users/tomaszskora/Documents/plgtskora/long_check/dna_986589_msd_av'
	
	# PREPARE TIME AND MSD
	time, m = prepare_time_and_msd_cm(start, end, steps, filename_template)
	
	# number_of = end - start + 1
	# m = [ np.zeros(steps, float) ]
	# counter = 0
	# for i in range(start, end + 1):
	# 	time, trajectory = read_trajectory('/Users/tomaszskora/Documents/plgtskora/long_check/dna_' + str(i) + '.xyz')
	# 	m += np.array( msd_cm(trajectory), float ) / number_of
	# 	counter += 1
	# 	print(str(counter) + '/' + str(number_of))
	# print('---')
	
	# PLOT D
	plot_d(time, m, filename_template_d + '.jpg')
	
	# plt.xlabel('time [ps]')
	# plt.ylabel(r'D [$\frac{\AA ^2}{ps}$]')
	# plt.plot(time[1:], m[0][1:] / np.array(time[1:]) / 6, '-', label = 'mean', lw = 1)
	# plt.savefig('/Users/tomaszskora/Documents/plgtskora/long_check/dna_d_av.jpg', dpi = 300)
	# plt.close()
	
	# PLOT MSD
	plot_msd(time, m, filename_template_msd + '.jpg')
	
	# plt.xlabel('time [ps]')
	# plt.ylabel(r'MSD [$\AA ^2$]')
	# plt.plot(time, m[0], '-', label = 'mean', lw = 1)
	# a, b = np.polyfit(time, m[0], 1)
	# diff = '{:7.6f}'.format(a / 6)
	# rh = '{:7.6f}'.format(6 * 0.0214 / a)
	# plt.text(len(time) / 2, 2000, r'$D =$' + diff + r' $\frac{\AA ^2}{ps}$')
	# plt.text(len(time) / 2, 20, r'$R_H =$' + rh + ' nm')
	# plt.plot(time, a * np.array(time, float) + b * np.ones(len(time)), '--', label = 'linear fit', lw = 1)
	# plt.legend()
	# plt.savefig('/Users/tomaszskora/Documents/plgtskora/long_check/dna_msd_av.jpg', dpi = 300)
	# plt.close()
	
	# SAVE MSD TO FILE
	save_msd_to_file(time, m, filename_template_msd + '.txt')
	
	# with open('/Users/tomaszskora/Documents/plgtskora/long_check/dna_traj_av.txt', 'w') as output_file:
	# 	for i in range(len(time)):
	# 		output_file.write(str(time[i])+' '+str(m[0][i])+'\n')
	
	###
	
	# time, m = read_msd_from_file('/Users/tomaszskora/Documents/plgtskora/long_check/dna_msd_av.txt')
	# plot_msd(time, m, '/Users/tomaszskora/Documents/plgtskora/long_check/dna_msd_av_corr.jpg')

if __name__ == '__main__':

	# elongation_of_dna_wrapped()
	hydrodynamic_radius_of_dna_wrapped()