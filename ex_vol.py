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

# TODO: Refactor that code and separate into functions
# TODO: Sphere volume/surface may be thrown into Sphere class
if __name__ == '__main__':

	# input_filename = parse_input_filename()
	# i = InputData(input_filename, [])
	# print(i.input_data)

	# box_size = i.input_data["box_size"]
	# r_crowders = i.input_data["crowders_radii"]
	# r_tracer = i.input_data["tracers_radii"]
	# number_of_trials = i.input_data["number_of_draws"]

	box_size = 750.0
	# r_tracer = [11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4, 11.4]
	# r_tracer = 39.2
	r_tracer = 51.0

	# number_of_trials = 100000
	number_of_trials = 100000

	input_labels = ["FIC", "DNA", "DNS"]
	input_radii = [51.0, 11.4, 39.2]

	xyz_filenames = []
	for which in range(1, 4):
		xyz_filenames += [ 'ficoll_39_DNA_39_'+str(which)+'.xyz' ]
		# xyz_filenames += [ 'ficoll_42_'+str(which)+'.xyz', 'ficoll_41_DNA_14_'+str(which)+'.xyz', 'ficoll_39_DNA_39_'+str(which)+'.xyz', 'ficoll_35_DNA_104_'+str(which)+'.xyz' ]
	# xyz_filenames = [ 'ficoll_41_DNA_14_'+str(which)+'.xyz', 'ficoll_39_DNA_39_'+str(which)+'.xyz', 'ficoll_35_DNA_104_'+str(which)+'.xyz' ]
	# xyz_filenames = [ 'ficoll_18_DNS_54_'+str(which)+'.xyz' ]
	# xyz_filenames = [ 'ficoll_18_DNS_54_'+str(which)+'.xyz' ]
	print( xyz_filenames )

	times = [0.0, 1000000.0, 2000000.0, 3000000.0, 4000000.0, 4500000.0]
	# times = [2000000.0, 3000000.0, 4000000.0]

	for xyz_filename in xyz_filenames:

		print(xyz_filename)

		for time in times:

			with open(xyz_filename, 'r') as xyz_file:

				labels, snapshot = _read_snapshot_from_xyz_file(xyz_file, time)

			r_crowders = []

			for label in labels:

				for i, input_label in enumerate( input_labels ):

					if label == input_label:

						r_crowders.append( input_radii[i] )

			crowders = place_crowders_xyz(r_crowders, snapshot)

			crowders = crowders[:-1]
	
			count = 0

			for crowder in crowders:
				crowder.r -= 1.5

			for i in tqdm( range(number_of_trials) ):
				tracers = place_tracers_linearly(r_tracer, box_size)
				for tracer in tracers:
					tracer.r -= 1.5
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
	
			ex_vol = count / number_of_trials * box_size**3
			print(ex_vol / box_size**3)