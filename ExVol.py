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
						 "input_xyz_template", "input_xyz_range", "times",
						 "input_str_filename", "radii_mode"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"bond_lengths":'hydrodynamic_radii', "withdraw":[], "float_type": 32,
				"scan_mode":False, "scan_density":2, "verbose":False, "debug":False,
				"omp_cores": 0}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataExVol(input_filename, required_keywords, defaults)
	i.input_data["labels"], i.input_data["crowder_radii"] = read_radii_from_str_file(
														i.input_data["input_str_filename"],
											 			i.input_data["radii_mode"] )

	timestamp( 'Input data:\n{}', i )

	if i.input_data["omp_cores"] == 0:
		nproc = multiprocessing.cpu_count()
	else:
		nproc = i.input_data["omp_cores"]
	print('You have {0:1d} CPUs'.format(nproc))

	results = []

	if i.input_data["scan_mode"]: scales = np.linspace(0.0, 1.0, i.input_data["scan_density"])
	else: scales = [1.0]

	for scale in scales:

		xyz_filenames = [ i.input_data["input_xyz_template"] + str(which) + '.xyz'
						  for which in range(*i.input_data["input_xyz_range"]) ]

		tracer_radii = list( scale * np.array( i.input_data["tracer_radii"] ) )
	
		if i.input_data["verbose"]: timestamp('input xyz filenames: {}', xyz_filenames)
	
		times_filenames = [ [time, filename] for time in i.input_data["times"] for filename in xyz_filenames ]
	
		tf_per_proc = len(times_filenames) // nproc
	
		tf_for_proc = [ times_filenames[m * tf_per_proc:( m + 1 ) * tf_per_proc] for m in range(nproc - 1) ] \
						+ [ times_filenames[ ( nproc - 1 ) * tf_per_proc:] ]
	
		pool = Pool(processes=nproc)
	
		estimate_excluded_volume_partial = partial( estimate_excluded_volume, input_labels = i.input_data["labels"], input_radii = i.input_data["crowder_radii"], r_tracer = tracer_radii, number_of_trials = i.input_data["number_of_trials"], box_size = i.input_data["box_size"], to_be_withdrawn = i.input_data["withdraw"] )
	
		excluded_volume = pool.map( estimate_excluded_volume_partial, tf_for_proc )
	
		if i.input_data["verbose"]: timestamp('single results: {}', excluded_volume)
	
		exvol = np.array([ excluded_volume[i][j][-1] for i in range( len(excluded_volume) ) for j in range( len(excluded_volume[i]) ) ])
	
		if i.input_data["debug"]: timestamp('single results: {}', exvol)
	
		print('r = {}; fex = {} +/- {}'.format(tracer_radii, np.mean(exvol), np.std(exvol, ddof=1)))
		print('R = {}'.format((np.mean(exvol)*i.input_data["box_size"]**3/(4/3*np.pi))**(1/3)))

		results.append( [ tracer_radii[0], np.mean(exvol), np.std(exvol, ddof=1) ] )

	with open(i.input_data["input_xyz_template"]+'fex.txt', 'w') as output_file:
		output_file.write('R_tr fex dfex\n')
		for result in results:
			output_file.write( '{} {} {}\n'.format(*result) )

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
