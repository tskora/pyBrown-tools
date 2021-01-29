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

# def cesp_macrolj_pointlj_potential(distance, dc, de, E, Ur, sigma, d0 = 10.0, eps = 0.8):

# 	b = dc - sigma

# 	if distance <= b: return np.inf

# 	if distance <= b + eps: distance = b + eps

# 	Vhard = 4*E*sigma**12 * 4*np.pi/45*b**3 * (15*distance**6 + 63*distance**4*b**2 + 45*distance**2*b**4 + 5*b**6) / (distance**2 - b**2)**9

# 	if de-dc != 0: Vsoft = Ur/2 * ( 1 - np.tanh(dc/d0/(de-dc)*(distance-(de+dc)/2)) )

# 	else: Vsoft = 0

# 	return Vhard + Vsoft

#-------------------------------------------------------------------------------

def compute_potential(tracers, crowders, box_size, potential, params):

	V = 0

	for tracer in tracers:

		for crowder in crowders:

			dc = params["dc"][crowder.label] + params["dc"]["tracer"]
			de = params["de"][crowder.label] + params["de"]["tracer"]
			Ur = np.sqrt( params["Ur"][crowder.label] * params["Ur"]["tracer"] )
			d0 = params["d0"]

			dist = _distance_pbc(tracer, crowder, box_size)

			V += potential(dist, dc = dc, de = de, Ur = Ur, d0 = d0)

	return V
