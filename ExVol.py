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

from pyBrown.input_ExVol import InputDataExVol
from pyBrown.ex_vol import estimate_excluded_volume, read_radii_from_str_file
from pyBrown.messaging import timestamp

import multiprocessing
from multiprocessing import Pool
from functools import partial

import numpy as np

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["box_size", "tracer_radii", "number_of_trials",
						 "xyz_templates", "xyz_range", "times",
						 "input_str_filename", "radii_mode"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"bond_lengths":'hydrodynamic_radii', "withdraw":[], "float_type": 32}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataExVol(input_filename, required_keywords, defaults)
	i.input_data["labels"], i.input_data["crowder_radii"] = read_radii_from_str_file(
														i.input_data["input_str_filename"],
											 			i.input_data["radii_mode"] )

	timestamp( 'Input data:\n{}', i )

	nproc = multiprocessing.cpu_count()
	print('You have {0:1d} CPUs'.format(nproc))

	xyz_filenames = []
	for template in i.input_data["xyz_templates"]:
		for which in range(*i.input_data["xyz_range"]):
			xyz_filenames += [ template + str(which) + '.xyz' ]

	times_filenames = [ [time, filename] for time in i.input_data["times"] for filename in xyz_filenames ]

	tf_per_proc = len(times_filenames) // nproc

	tf_for_proc = [ times_filenames[m * tf_per_proc:( m + 1 ) * tf_per_proc] for m in range(nproc - 1) ] \
					+ [ times_filenames[ ( nproc - 1 ) * tf_per_proc:] ]

	pool = Pool(processes=nproc)

	estimate_excluded_volume_partial = partial( estimate_excluded_volume, input_labels = i.input_data["labels"], input_radii = i.input_data["crowder_radii"], r_tracer = i.input_data["tracer_radii"], number_of_trials = i.input_data["number_of_trials"], box_size = i.input_data["box_size"], to_be_withdrawn = i.input_data["withdraw"] )

	excluded_volume = pool.map( estimate_excluded_volume_partial, tf_for_proc )

	excluded_volume_sorted = [ [] for i in range(len(i.input_data["xyz_templates"])) ]

	for element in excluded_volume:
		for subelement in element:
			words = subelement[1][:-4].split('_')[:-1]
			template = ''
			for word in words:
				template += word + '_'
			counter = 0
			for xyz_template in i.input_data["xyz_templates"]:
				if template == xyz_template:
					excluded_volume_sorted[counter].append( subelement[2] )
					break
				else:
					counter += 1

	for template, exvol in zip(i.input_data["xyz_templates"], excluded_volume_sorted):
		print(exvol)
		print('{}: fex = {} +/- {}'.format(template[:-1], np.mean(exvol), np.std(exvol, ddof=1)))



#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
