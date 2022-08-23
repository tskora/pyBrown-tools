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

def compute_flux_single(r0, r1, plane_normal_vector, plane_point, box_size):

	for i in range(0, 3):
		if ( r1[i] - r0[i] >= box_size / 2 ) or ( r1[i] - r0[i] <= -box_size / 2 ):
			return 0

	f0 = np.dot( plane_normal_vector, (r0 - plane_point) ) > 0.0
	f1 = np.dot( plane_normal_vector, (r1 - plane_point) ) > 0.0

	if f0 == f1: return 0
	elif f0: return -1
	else: return 1