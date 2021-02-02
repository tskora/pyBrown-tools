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

from pyBrown.sphere import _distance_pbc

def hs_potential(distance, d):

	if distance <= d: return np.inf

	else: return 0

#-------------------------------------------------------------------------------

def cesp_hs_potential(distance, dc, de, Ur, d0):

	if distance <= dc: return np.inf

	if de - dc != 0: return Ur/2 * ( 1 - np.tanh(dc/d0/(de-dc)*(distance-(de+dc)/2)) )

	else: return 0

#-------------------------------------------------------------------------------

def cesp_macrolj_pointlj_potential(distance, dc, de, E, Ur, sigma, d0, eps):

	b = dc - sigma

	if distance <= b: return np.inf

	if distance <= b + eps: distance = b + eps

	Vhard = 4*E*sigma**12 * 4*np.pi/45*b**3 * (15*distance**6 + 63*distance**4*b**2 + 45*distance**2*b**4 + 5*b**6) / (distance**2 - b**2)**9

	if de-dc != 0: Vsoft = Ur/2 * ( 1 - np.tanh(dc/d0/(de-dc)*(distance-(de+dc)/2)) )

	else: Vsoft = 0

	return Vhard + Vsoft

#-------------------------------------------------------------------------------

# def cesp_macrolj_macrolj_potential(distance, dc1, dc2, de, E, Ur, sigma, d0, eps):

# 	b = dc1 + dc2 - sigma

# 	dc = dc1 + dc2

# 	if distance <= b: return np.inf

# 	if distance <= b + eps: distance = b + eps

# 	Vhard = E*np.pi**2/315.0*dc1*dc2/(dc)*sigma**6/(distance-dc)**7

# 	if de-dc != 0: Vsoft = Ur/2 * ( 1 - np.tanh(dc/d0/(de-dc)*(distance-(de+dc)/2)) )

# 	else: Vsoft = 0

# 	return Vhard + Vsoft

#-------------------------------------------------------------------------------

def compute_potential(tracers, crowders, box_size, potential_name, params):

	# import matplotlib.pyplot as plt

	# dc1 = 51
	# dc2 = 0.1
	# de = 51
	# E = 0.39
	# Ur = 0.67
	# eps = 0.8
	# sigma = 1.5
	# d0 = 10.0

	# rs = np.linspace(49, 150, 1000)
	# u1 = np.array([cesp_macrolj_pointlj_potential(r, dc1, de, E, Ur, sigma, d0, eps) for r in rs])
	# u2 = np.array([cesp_macrolj_macrolj_potential(r, dc1, dc2, de, E, Ur, sigma, d0, eps) for r in rs])
	# _ = plt.plot(rs, u1, '-', color = 'red', label = 'point')
	# _ = plt.plot(rs, u2, '-.', color = 'blue', label = 'limit')

	# _ = plt.legend()

	# _ = plt.show()

	# 1/0

	if potential_name == "cesp-hs":

		potential = cesp_hs_potential

	elif potential_name == "cesp-macrolj" and params["dc"]["tracer"] == 0:

		potential = cesp_macrolj_pointlj_potential

	elif potential_name == "cesp-macrolj" and params["dc"]["tracer"] != 0:

		potential = cesp_macrolj_macrolj_potential

	else:

		return None

	V = 0

	for tracer in tracers:

		for crowder in crowders:

			if potential_name == "cesp-hs":
				dc = params["dc"][crowder.label] + params["dc"]["tracer"]
				de = params["de"][crowder.label] + params["de"]["tracer"]
				Ur = np.sqrt( params["Ur"][crowder.label] * params["Ur"]["tracer"] )
				d0 = params["d0"]
				pot_params = [dc, de, Ur, d0]

			dist = _distance_pbc(tracer, crowder, box_size)

			V += potential(dist, *pot_params)

	return V