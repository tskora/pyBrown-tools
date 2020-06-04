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

from pyBrown.input import InputData
from pyBrown.parse import parse_input_filename
from pyBrown.monte_carlo import MonteCarlo, place_crowders_linearly, place_tracers_linearly, place_crowders_xyz
from pyBrown.sphere import Sphere, overlap

import numpy as np

from tqdm import tqdm


def _read_snapshot_from_xyz_file(xyz_file, time):

	snapshot = [  ]
	labels = [  ]

	number_of_beads_per_file = int( xyz_file.readline() )

	this_one = False

	for i, line in enumerate(xyz_file):

		is_line_unimportant = len( line.split() ) == 1

		if 'xyz' in line.split()[0].split('.'):

			time_of_snapshot = float( line.split()[3] )

			if time == time_of_snapshot:

				this_one = True

			else:

				this_one = False

		elif is_line_unimportant: continue

		else:

			if this_one:

				labels.append( line.split()[0] )

				coords = np.array( [ float(line.split()[x]) for x in range(1, 4) ], float )
						
				snapshot.append( coords )

	return labels, snapshot

#-------------------------------------------------------------------------------

def estimate_excluded_volume(tfs, input_labels, input_radii, r_tracer, number_of_trials, box_size, to_be_withdrawn):

	result = []

	for time_filename in tfs:

		time, xyz_filename = time_filename

		with open(xyz_filename, 'r') as xyz_file:

			labels, snapshot = _read_snapshot_from_xyz_file(xyz_file, time)

		import random

		# print('before\n{}\n'.format(labels))

		while len(to_be_withdrawn) > 0:

			random_index = random.randint( 0, len(labels) - 1 )

			if labels[random_index] == to_be_withdrawn[0]:

				del labels[random_index]

				del snapshot[random_index]

				del to_be_withdrawn[0]

		# print('after\n{}\n'.format(labels))

		r_crowders = []

		for label in labels:

			for i, input_label in enumerate( input_labels ):

				if label == input_label:

					r_crowders.append( input_radii[i] )

		crowders = place_crowders_xyz(r_crowders, snapshot)

		# crowders = crowders[:-1]
		# crowders = crowders[1:]
	
		count = 0

		# for crowder in crowders:
		# 	crowder.r -= 1.5

		for i in tqdm( range(number_of_trials) ):

			tracers = place_tracers_linearly(r_tracer, box_size)

			# for tracer in tracers:
			# 	tracer.r -= 1.5

			if overlap(crowders, tracers, 0.0):
				count += 1

			else:
				versors = [ np.array([nx * box_size,
									  ny * box_size,
									  nz * box_size])
							for nx in np.arange(-1, 2, 1)
							for ny in np.arange(-1, 2, 1)
							for nz in np.arange(-1, 2, 1) ]

				if_overlap = False

				for versor in versors:

					for tracer in tracers:
						tracer.translate( versor )

					if overlap(crowders, tracers, 0.0):
						if_overlap = True

					for tracer in tracers:
						tracer.translate( -versor )

					if if_overlap:
						count += 1
						break

		ex_vol = count / number_of_trials

		result.append( [time, xyz_filename, ex_vol] )

	return result