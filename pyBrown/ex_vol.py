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
import random

from tqdm import tqdm

def read_radii_from_str_file(input_str_filename, mode):

	with open(input_str_filename) as str_file:
		labels = []
		radii = []
		for line in str_file:
			if line.split()[0] == 'sub':
				label = line.split()[1]
				hydrodynamic_radius = float(line.split()[6])
				lennard_jones_radius = float(line.split()[8]) / 2

				if label not in labels:
					labels.append(label)
					if mode == 'hydrodynamic': radii.append( hydrodynamic_radius )
					elif mode == 'lennard-jones': radii.append( lennard_jones_radius )
					else: print('radii mode unknown')

	return labels, radii

#-------------------------------------------------------------------------------

def _read_snapshot_from_xyz_file(xyz_filename, snapshot_time):

	snapshot = [ ]
	labels = [ ]

	start_snapshot = False

	with open(xyz_filename, 'r') as xyz_file:
		for line in xyz_file:
			if start_snapshot:
				# 'time' here means it is already the next snapshot, 
                # so we can break
				if 'time' in line.split():
					break
				elif len( line.split() ) == 4:
					label = line.split()[0]
					labels.append(label)
					coords = [ line.split()[i] for i in range(1,4) ]
					snapshot.append( coords )
			if 'time' in line.split():
				time = float(line.split()[-1])
				if time >= snapshot_time:
					# print('Saving snapshot for time %f (requested: %f)' % (time, snapshot_time))
					start_snapshot = True

	return labels, snapshot

#-------------------------------------------------------------------------------

def estimate_excluded_volume(tfs, input_labels, input_radii, r_tracer, number_of_trials, box_size, to_be_withdrawn):

	result = []

	for time_filename in tfs:

		time, xyz_filename = time_filename

		labels, snapshot = _read_snapshot_from_xyz_file(xyz_filename, time)

		# print('before\n{}\n'.format(labels))

		tbw = to_be_withdrawn.copy()

		while len(tbw) > 0:

			random_index = random.randint( 0, len(labels) - 1 )

			if labels[random_index] == tbw[0]:

				del labels[random_index]

				del snapshot[random_index]

				del tbw[0]

		# print('after\n{}\n'.format(labels))

		r_crowders = []

		for label in labels:

			for i, input_label in enumerate( input_labels ):

				if label == input_label:

					r_crowders.append( input_radii[i] )

		crowders = place_crowders_xyz(r_crowders, snapshot)
	
		count = 0

		for i in tqdm( range(number_of_trials) ):

			tracers = place_tracers_linearly(r_tracer, box_size)

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