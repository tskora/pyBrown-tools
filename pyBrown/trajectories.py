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

from scipy.constants import Boltzmann

from pyBrown.input import print_input
from pyBrown.messaging import timestamp
from pyBrown.plot_config import plot_config

#-------------------------------------------------------------------------------

def read_trajectories(input_data):

	input_xyz_filenames = [ input_data["input_xyz_template"] + str(i) + '.xyz' 
							for i in range( *input_data["input_xyz_range"] ) ]

	number_of_xyz_files = len( input_xyz_filenames )

	numbers_of_timeframes = [ _count_timeframes( input_xyz_filenames[i],
											  input_data["probing_frequency"] )
							  for i in range( number_of_xyz_files ) ]

	number_of_timeframes = min( numbers_of_timeframes )

	number_of_beads = _count_beads( input_xyz_filenames[0] )

	for i in range( 1, number_of_xyz_files ):

		assert _count_beads( input_xyz_filenames[i] ) == number_of_beads

	number_of_xyz_files = len( input_xyz_filenames )

	if input_data["debug"]:
		print_input( input_data )
		print( 'Input xyz filenames: {}'.format(input_xyz_filenames) )
	if input_data["debug"] or input_data["verbose"]:
		print( 'Number of timeframes: {}'.format(number_of_timeframes) )
		print( 'Number of beads: {}'.format(number_of_beads) )
		print( 'Number of files: {}'.format(number_of_xyz_files) )

	# temporary binary file which will contain the trajectories
	temporary_filename = 'temp.dat'

	trajectories = np.memmap( temporary_filename, dtype = np.float32,
							  mode = 'w+',
							  shape = ( number_of_beads * number_of_xyz_files,
									    number_of_timeframes, 3 ) )

	for i in range( number_of_beads * number_of_xyz_files ):
		for j in range(number_of_timeframes):
			trajectories[i, j, :] = np.zeros(3, np.float32)

	del trajectories

	times = np.zeros( number_of_timeframes, dtype = np.float32 )
	labels = [ '___' for _ in range( number_of_xyz_files * number_of_beads ) ]

	for i, input_xyz_filename in enumerate(input_xyz_filenames):

		if input_data["verbose"] or input_data["debug"]:
			timestamp( 'XYZ file {} / {}', i + 1, number_of_xyz_files )

		with open(input_xyz_filename, 'r') as input_xyz_file:

			temp = np.memmap( temporary_filename, dtype = np.float32,
						   	  shape = ( number_of_beads * len(input_xyz_filenames),
						      			number_of_timeframes, 3 ) )

			trajectories = temp[i * number_of_beads : (i + 1) * number_of_beads, :, :]

			_read_trajectories_from_xyz_file(input_xyz_file, trajectories, times,
											 labels, input_data["probing_frequency"])

			del temp

	return temporary_filename, times, labels

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

def compute_msds(input_data, temporary_filename_2, cm_labels):

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

	### FREUD VERSION START

	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
						   shape = ( len(cm_labels), number_of_timeframes, 3 ) )

	unify_coordinates(temp2, input_data["box_size"])

	trajs_all = []

	for k in range(len(input_data["labels"])):

		trajs_list = [ [ temp2[i][j] for i in range( len(cm_labels) ) if input_data["labels"][k] == cm_labels[i]  ] for j in range( number_of_timeframes ) ]

		trajs = np.array( trajs_list )

		trajs_all.append(trajs)

	import freud.box
	import freud.msd

	msds = []

	for k in range(len(input_data["labels"])):

		box = freud.box.Box.cube(750.0)
		box2 = freud.box.Box.cube(750.0)

		msd = freud.msd.MSD(box, 'direct')
		msd2 = freud.msd.MSD(box2, 'window')
		msd.compute( positions = trajs_all[k] )
		msd2.compute( positions = trajs_all[k] )

		if (input_data["mode"] == "direct"):
			msds.append(msd.msd)
		elif (input_data["mode"] == "window"):
			msds.append(msd2.msd)

		# plt.plot(msd.msd, 'o', label = 'direct')
		# plt.plot(msd2.msd, '+', label = 'window')

	# plt.title('Mean Squared Displacement')
	# plt.xlabel('$t$')
	# plt.ylabel('MSD$(t)$')
	# plt.legend()
	# plt.savefig('test.jpg', dpi=300)

	# print(msd.msd)
	# print(msd2.msd)

	# plt.show()

	# 1/0

	### FREUD VERSION END

	### OLD VERSION START

	# del trajs

	# del temp2

	# ###

	# for i in range( len(cm_labels) ):

	# 	if input_data["verbose"]: print('computing sd: {} / {}'.format(i, len(cm_labels)))

	# 	temp2 = np.memmap( temporary_filename_2, dtype = np.float32,
	# 					   shape = ( len(cm_labels), number_of_timeframes, 3 ) )

	# 	cm_trajectory = temp2[i]

	# 	del temp2

	# 	sd = _compute_sd(cm_trajectory, input_data["box_size"])

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
		###
		# M1 = [0., 1246.0022, 2532.1602, 3825.5596, 5116.626, 6366.3667,\
		# 7579.986,    8753.173,   10035.295,   11044.939,   12464.446,   13869.375,\
		# 14804.468,   16271.135,   17103.334,   18242.799,   19509.38 ,   20975.828,\
		# 22015.893,   23295.648,   24582.496,   25975.95 ,   27107.523,   28601.045,\
		# 30230.736,   31809.893,   32862.04 ,   34215.918,   35210.367,   36058.445,\
		# 37397.77 ,   38873.633,   39922.15 ,   41460.73 ,   43001.266,   44665.52,\
		# 45899.973,   46985.992,   48113.758,   49475.652,   50679.477,   52023.867,\
		# 52793.035,   53951.918,   55073.332,   56086.375,   57149.48 ,   58714.406,\
		# 59968.043,   61164.332,   62932.805,   64342.75 ,   65396.984,   66588.27,\
		# 67911.53 ,   69166.49 ,   70285.03 ,   71828.59 ,   73010.58 ,   74359.26,\
		# 75756.12 ,   77339.64 ,   78576.42 ,   79267.33 ,   80428.04 ,   81985.61,\
		# 82749.34 ,   84365.22 ,   86237.586,   87111.055,   88242.44 ,   89471.195,\
		# 90779.516,   91899.15 ,   93661.97 ,   94915.21 ,   96089.89 ,   97442.16,\
		# 98179.03 ,   99122.52 ,  100565.31 ,  101971.63 ,  103023.305,  104673.18,\
		# 106190.29 ,  107579.914,  109123.305,  110563.41 ,  111790.234,  112988.41,\
		# 114616.664,  115803.555,  117083.35 ,  119206.336,  120394.05 ,  121738.15,\
		# 123370.03 ,  124204.65 ,  125741.41 ,  126307.66 ,  127838.66 ,  129446.96,\
		# 130932.586,  132259.05 ,  134096.33 ,  135243.62 ,  136783.42 ,  137894.98,\
		# 138973.8  ,  139796.28 ,  140380.38 ,  141418.81 ,  142595.86 ,  143493.11,\
		# 144297.69 ,  145539.69 ,  146826.34 ,  147807.3  ,  149334.88 ,  149926.06,\
		# 151321.75 ,  152607.62 ,  153540.27 ,  154351.86 ,  155831.47 ,  156805.73,\
		# 158000.52 ,  159072.8  ,  160179.06 ,  161219.52 ,  162102.42 ,  163154.27,\
		# 164814.83 ,  165668.67 ,  166246.75 ,  167451.73 ,  168869.36 ,  170067.42,\
		# 171586.14 ,  172582.97 ,  173771.06 ,  175669.31 ,  176895.98 ,  177741.56,\
		# 178144.52 ,  179595.98 ,  180744.84 ,  181897.12 ,  182716.72 ,  183575.6,\
		# 184376.16 ,  185454.88 ,  186773.03 ,  188430.53 ,  190201.78 ,  191279.28,\
		# 191763.48 ,  193033.17 ,  194146.39 ,  195841.53 ,  197517.42 ,  198998.27,\
		# 200101.69 ,  201258.53 ,  202317.86 ,  203344.94 ,  204459.84 ,  205859.28,\
		# 206596.73 ,  207479.39 ,  209463.44 ,  211008.44 ,  212376.2  ,  213445.52,\
		# 214535.34 ,  216172.03 ,  217436.5  ,  218389.33 ,  219291.28 ,  219466.,\
		# 220881.38 ,  222059.48 ,  224078.02 ,  225536.17 ,  226562.23 ,  228279.14,\
		# 229326.67 ,  230481.06 ,  231739.9  ,  233774.78 ,  234496.25 ,  235599.14,\
		# 236731.61 ,  238121.38 ,  239102.88 ,  241044.58 ,  242547.7  ,  243218.31,\
		# 244784.36 ,  246451.78 ]
		# M2 = [2.74313486e-03, 1.24715484e+03, 2.48381490e+03, 3.71929854e+03,\
		# 4.94982572e+03, 6.17564631e+03, 7.40140089e+03, 8.62852459e+03,\
		# 9.85367837e+03, 1.10815858e+04, 1.23085221e+04, 1.35349160e+04,\
		# 1.47539785e+04, 1.59735761e+04, 1.71964367e+04, 1.84195795e+04,\
		# 1.96440735e+04, 2.08688894e+04, 2.20913731e+04, 2.33195665e+04,\
		# 2.45516099e+04, 2.57782713e+04, 2.70070729e+04, 2.82378716e+04,\
		# 2.94681800e+04, 3.06929469e+04, 3.19172546e+04, 3.31444588e+04,\
		# 3.43756874e+04, 3.56045305e+04, 3.68393909e+04, 3.80758210e+04,\
		# 3.93092208e+04, 4.05474825e+04, 4.17872714e+04, 4.30287591e+04,\
		# 4.42737299e+04, 4.55167034e+04, 4.67585105e+04, 4.80010430e+04,\
		# 4.92391994e+04, 5.04757705e+04, 5.17111132e+04, 5.29515880e+04,\
		# 5.41952216e+04, 5.54369951e+04, 5.66793486e+04, 5.79234098e+04,\
		# 5.91677697e+04, 6.04097508e+04, 6.16465879e+04, 6.28799587e+04,\
		# 6.41117968e+04, 6.53463693e+04, 6.65832058e+04, 6.78238650e+04,\
		# 6.90692785e+04, 7.03111085e+04, 7.15583058e+04, 7.28082151e+04,\
		# 7.40571983e+04, 7.53050311e+04, 7.65480983e+04, 7.77939268e+04,\
		# 7.90432228e+04, 8.02889234e+04, 8.15343149e+04, 8.27822228e+04,\
		# 8.40331640e+04, 8.52741149e+04, 8.65201266e+04, 8.77627843e+04,\
		# 8.90062372e+04, 9.02518280e+04, 9.15064211e+04, 9.27545597e+04,\
		# 9.40051356e+04, 9.52622255e+04, 9.65152662e+04, 9.77750487e+04,\
		# 9.90414639e+04, 1.00309860e+05, 1.01580695e+05, 1.02845179e+05,\
		# 1.04109397e+05, 1.05376984e+05, 1.06640826e+05, 1.07901690e+05,\
		# 1.09150722e+05, 1.10410195e+05, 1.11668894e+05, 1.12923430e+05,\
		# 1.14168104e+05, 1.15409423e+05, 1.16650125e+05, 1.17891038e+05,\
		# 1.19130629e+05, 1.20360977e+05, 1.21593274e+05, 1.22815158e+05,\
		# 1.24040651e+05, 1.25260096e+05, 1.26476270e+05, 1.27680110e+05,\
		# 1.28886567e+05, 1.30090506e+05, 1.31302698e+05, 1.32513323e+05,\
		# 1.33722888e+05, 1.34937113e+05, 1.36161214e+05, 1.37380175e+05,\
		# 1.38594444e+05, 1.39807493e+05, 1.41021094e+05, 1.42242689e+05,\
		# 1.43461699e+05, 1.44686444e+05, 1.45920825e+05, 1.47163487e+05,\
		# 1.48417598e+05, 1.49662162e+05, 1.50909296e+05, 1.52155617e+05,\
		# 1.53408091e+05, 1.54653759e+05, 1.55892954e+05, 1.57138773e+05,\
		# 1.58398073e+05, 1.59648252e+05, 1.60902932e+05, 1.62149561e+05,\
		# 1.63395737e+05, 1.64630774e+05, 1.65868960e+05, 1.67123626e+05,\
		# 1.68390011e+05, 1.69654942e+05, 1.70922672e+05, 1.72185897e+05,\
		# 1.73448029e+05, 1.74703904e+05, 1.75938296e+05, 1.77158045e+05,\
		# 1.78390794e+05, 1.79633904e+05, 1.80871902e+05, 1.82112928e+05,\
		# 1.83351231e+05, 1.84586591e+05, 1.85837536e+05, 1.87087811e+05,\
		# 1.88346347e+05, 1.89602637e+05, 1.90853864e+05, 1.92095503e+05,\
		# 1.93339658e+05, 1.94599944e+05, 1.95857271e+05, 1.97144913e+05,\
		# 1.98422388e+05, 1.99657681e+05, 2.00874319e+05, 2.02071932e+05,\
		# 2.03269076e+05, 2.04458361e+05, 2.05667875e+05, 2.06885203e+05,\
		# 2.08111733e+05, 2.09341987e+05, 2.10557514e+05, 2.11735853e+05,\
		# 2.12899396e+05, 2.14046228e+05, 2.15225112e+05, 2.16431356e+05,\
		# 2.17655900e+05, 2.18882797e+05, 2.20113370e+05, 2.21365518e+05,\
		# 2.22676132e+05, 2.23961250e+05, 2.25254075e+05, 2.26522305e+05,\
		# 2.27784577e+05, 2.29046768e+05, 2.30258695e+05, 2.31494138e+05,\
		# 2.32753685e+05, 2.34068071e+05, 2.35265377e+05, 2.36406683e+05,\
		# 2.37530669e+05, 2.38721050e+05, 2.39945570e+05, 2.41259878e+05,\
		# 2.42600813e+05, 2.43652594e+05, 2.44875710e+05, 2.46450412e+05]
		# plt.plot( times / 1000000, M1, 'o', label = input_data["labels"][i], color = 'black' )
		# plt.plot( times / 1000000, M2, '+', label = input_data["labels"][i], color = 'blue' )
		###

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

		if counter // number_of_beads // probing_frequency >= len(times): break

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

	if frequency == 1:
		return _file_length(filename) // (number_of_beads + 2)
	else:
		return _file_length(filename) // (number_of_beads + 2) // frequency# + 1

#-------------------------------------------------------------------------------

def unify_coordinates(trajectories, box_length):

	for particle in trajectories:
		for i in range(1, len(particle)):
			_coord_unify(particle[i-1], particle[i], box_length)

#-------------------------------------------------------------------------------

def _coord_unify(past, present, box_length):

	for i in range(3):
		while( present[i] - past[i] >= box_length * np.sqrt(3) / 2 ):
			present[i] -= box_length
		while( past[i] - present[i] >= box_length * np.sqrt(3) / 2 ):
			present[i] += box_length