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

import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import Boltzmann
from scipy.optimize import curve_fit

from pyBrown_tools.plot_config import plot_config

#-------------------------------------------------------------------------------

def MSD_fit_func(msd, D):

	return msd/6/D

#-------------------------------------------------------------------------------

def plot_msds(input_data, times, msds):

	colors, _ = plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'MSD [$\AA ^2$]')

	for i, msd in enumerate(msds):

		plt.plot( times / 1000000, msd, '-', label = input_data["labels"][i], color = colors[i] )

		if input_data["fit_MSD"]:
			# a, b = np.polyfit(times, msd, 1)
			D = curve_fit(MSD_fit_func, msd, times)[0][0]
			# D = a / 6 # in A**2 / ps
			D_string = '{:7.5f}'.format(D * 1000) # in AA**2 / ns
			Rh = 10**17 * Boltzmann / (6 * np.pi) * input_data["temperature"] / ( D * 0.1 * input_data["viscosity"] ) # in nm
			Rh_string = '{:4.2f}'.format(Rh) # in nm

			plt.text( 0.0 * times[-1], (0.70 - i * 0.15) * np.ndarray.max( np.array( msds ) ), input_data["labels"][i] )
			plt.text( 0.0 * times[-1], (0.65 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$D =$' + D_string + r' $\frac{\AA ^2}{ns}$')
			plt.text( 0.0 * times[-1], (0.60 - i * 0.15) * np.ndarray.max( np.array( msds ) ), r'$R_H =$' + Rh_string + ' nm')

			plt.plot(times / 1000000, 6 * D * times, '--', label = 'linear fit for ' + input_data["labels"][i], color = colors[i]  )
			# plt.plot(times / 1000000, a * np.array(times, float) + b * np.ones(len(times)), '--', label = 'linear fit for ' + input_data["labels"][i], color = colors[i]  )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'msd.pdf')

	plt.close()

#-------------------------------------------------------------------------------

def plot_msads(input_data, times, msads):

	colors, _ = plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'MSAD [$rad ^2$]')

	for i, msad in enumerate(msads):

		plt.plot( times / 1000000, msad, '-', label = input_data["labels"][i], color = colors[i%5] )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'msad.pdf')

	plt.close()

#-------------------------------------------------------------------------------

def plot_moas(input_data, times, moas):

	colors, _ = plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'MOA')

	for i, moa in enumerate(moas):

		plt.plot( times / 1000000, moa, '-', label = input_data["labels"][i], color = colors[i%5] )

	plt.legend()

	plt.savefig(input_data["input_xyz_template"] + 'moa.pdf')

	plt.close()

#-------------------------------------------------------------------------------

def plot_menergies(input_data, times, menergies):

	plot_config()

	plt.xlabel(r't [$\mu s$]')
	plt.ylabel(r'E [$\frac{kcal}{mol}$]')

	plt.plot(times / 1000000, menergies, '-')

	plt.savefig(input_data["input_enr_template"] + 'enr.pdf')

	plt.close()

#-------------------------------------------------------------------------------

def plot_pore_size_distribution(input_data, bins, psd):

	plot_config()

	plt.xlabel(r'd [$\AA$]')
	plt.ylabel(r'PSD')

	bins += 0.5*(bins[1] - bins[0])

	plt.plot(bins[:-1], psd, '-')

	plt.savefig(input_data["input_xyz_template"] + 'psd.pdf')

	plt.close()
