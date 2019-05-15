# pyVrown is a bound of tools useful for Brownian and Stokesian dynamics simulations
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

from pyBrown.input import InputData
from pyBrown.trajectories import read_trajectories, read_energies, \
								 separate_center_of_mass, \
								 compute_msds, compute_menergies, \
								 plot_msds, plot_menergies
from pyBrown.parse import parse_input_filename

#-------------------------------------------------------------------------------

if __name__ == '__main__':

	required_keywords = []

	input_filename = parse_input_filename()
	i = InputData(input_filename, required_keywords)

	# energies, times = read_energies(i)
	# menergies = compute_menergies(energies)
	# del energies
	# plot_menergies(i, times, menergies)
	# del times
	# del menergies

	temporary_filename, times, labels = read_trajectories(i)
	temporary_filename_2, cm_labels = separate_center_of_mass(i, temporary_filename, labels)
	del labels
	msds = compute_msds(i, temporary_filename_2, cm_labels)
	del cm_labels
	plot_msds(i, times, msds)
	del times
	del msds
