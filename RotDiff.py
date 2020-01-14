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

import numpy as np

from pyBrown.input_RotDiff import InputDataRotDiff
from pyBrown.messaging import timestamp
from pyBrown.trajectories import read_trajectories, add_auxiliary_data_multibeads, \
								 compute_orientations, compute_mean_orientation_autocorrelation, \
								 save_mean_orientation_autocorrelation_to_file, \
								 compute_mean_squared_angular_displacements, \
								 save_mean_squared_angular_displacements_to_file, \
								 compute_nematic_order, compute_mean_director

#-------------------------------------------------------------------------------

@click.command()
@click.argument('input_filename',
				type = click.Path( exists = True ))
def main(input_filename):

	# here the list of keywords that are required for program to work is provided
	required_keywords = ["labels", "sizes", "box_size", "temperature", "viscosity",
						 "input_xyz_template", "input_xyz_range"]

	# here the dict of keywords:default values is provided
	# if given keyword is absent in JSON, it is added with respective default value
	defaults = {"debug": False, "verbose": False, "fit_MSD": False,
				"probing_frequency": 1, "min_time": 0.0, "mode": "window"}

	required_keywords = []
	defaults = {}

	timestamp( 'Reading input from {} file', input_filename )
	i = InputDataRotDiff(input_filename, required_keywords, defaults)
	timestamp( 'Input data:\n{}', i )

	timestamp( 'Reading trajectories' )
	times, labels, auxiliary_data = read_trajectories(i.input_data)
	add_auxiliary_data_multibeads( i.input_data, labels, auxiliary_data )
	timestamp( 'Separating the rotational movement' )
	orientation_labels = compute_orientations( i.input_data, labels, auxiliary_data )
	del labels
	timestamp( 'Computing director' )
	mean_director =  compute_mean_director( i.input_data, orientation_labels, auxiliary_data )
	timestamp( 'Computing nematic order parameters' )
	a = compute_nematic_order( i.input_data, orientation_labels, auxiliary_data, director = np.array([1.0, 0.0, 0.0], dtype=np.float32) )
	b = compute_nematic_order( i.input_data, orientation_labels, auxiliary_data, director = np.array([0.0, 1.0, 0.0], dtype=np.float32) )
	c = compute_nematic_order( i.input_data, orientation_labels, auxiliary_data, director = np.array([0.0, 0.0, 1.0], dtype=np.float32) )
	d = compute_nematic_order( i.input_data, orientation_labels, auxiliary_data, director = mean_director[1] )
	print( mean_director )
	print(a)
	print(b)
	print(c)
	print(d)
	1/0
	print( compute_nematic_order( i.input_data, orientation_labels, auxiliary_data, director = mean_director[1] ) )
	1/0

	timestamp( 'Computing mean orientation autocorrelations' )
	moas = compute_mean_orientation_autocorrelation( i.input_data, orientation_labels, auxiliary_data )
	timestamp( 'Saving mean orientation autocorrelations to a file' )
	save_mean_orientation_autocorrelation_to_file(i.input_data, times, moas)
	timestamp( 'Computing mean squared angular displacements' )
	msads = compute_mean_squared_angular_displacements( i.input_data, orientation_labels, auxiliary_data )
	timestamp( 'Saving mean squared angular displacements to a file' )
	save_mean_squared_angular_displacements_to_file(i.input_data, times, msads)
	del orientation_labels

	import matplotlib.pyplot as plt

	divider = 10

	ys = moas[1]
	a, b = np.polyfit(times[:len(times)//divider], np.log( ys[:len(times)//divider] ), 1)
	plt.plot( times, np.log( ys ), '-')
	plt.plot(times, a * times + b * np.ones(len(times)), ':')
	plt.show()
	plt.close()
	print(-0.5*a)

	zs = msads[1]
	aa, bb = np.polyfit( times[:len(times)//divider], zs[:len(times)//divider], 1 )
	plt.plot( times, zs, '-' )
	plt.plot( times, aa * times + bb * np.ones(len(times)), ':' )
	plt.show()
	plt.close()
	print(0.25*aa)

	1/0

	# timestamp( 'Plotting mean square displacements' )
	# plot_msds(i.input_data, times, msds)
	del times
	del moas
	del msads

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	
	main()
