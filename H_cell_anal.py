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

import click
import numpy as np

import matplotlib.pyplot as plt

from pyBrown.h_cell import read_cells_from_files,\
						   perform_differential_analysis_for_2_substances,\
						   perform_analysis_for_many_substances

#-------------------------------------------------------------------------------

@click.command()
@click.option(
	'--inp', '-i',
	help = 'template of input filename (output filenames from ' +
		   'H cell sim is organized in a following way:\n' +
		   '"template"_"INT"_snapshot_"FLOAT".txt where in place ' +
		   'of "template" there is chosen job name, "INT" is used to ' +
		   'enumerate substances/diffusion coefficients and "FLOAT" ' +
		   'is the value of x in given snapshot'
	)
@click.option(
	'--snapshots', '-s',
	help = 'x values for concentration profiles to be plotted ' +
		   '(you may provide a range with syntax a:b:c where a ' +
		   'is the minimal value, b is the maximal value and c ' +
		   'is the number of equally distributed values)'
	)
@click.option(
	'--num-subs', '-n',
	default = 2,
	help = 'number of substances (different diffusion constants)')
@click.option(
	'--num-outlets', '-x',
	default = 2,
	help = 'number of H-cell outlets')
@click.option(
	'--out-idx', '-d',
	default = 1,
	help = 'index of outlet to be followed (indexing starts from 1)')
def main(inp, snapshots, num_subs, num_outlets, out_idx):

	ss = parse_snapshots(snapshots)

	hs, Ds = read_cells_from_files(inp, ss, num_subs)

	if num_subs == 2:

		perform_differential_analysis_for_2_substances(hs, Ds, ss,
			num_outlets, inp)

	else:

		perform_analysis_for_many_substances(hs, Ds, ss,
			num_outlets, out_idx, inp)

#-------------------------------------------------------------------------------

def parse_snapshots(snapshots):

	if not ( ":" in snapshots ):

		return [ float( element ) for element in snapshots.split() ]

	else:

		return list( np.linspace( *[ float( snapshots.split(':')[i] )
			for i in range(3) ] ) )

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
