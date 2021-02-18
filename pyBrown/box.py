# pyBrown is a bound of tools useful for Brownian and Stokesian dynamics simulations
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

from pyBrown.hydrodynamics import M_rpy#, M_rpy_smith, R_lub_corr

import math
import numpy as np
from scipy.constants import Boltzmann

from pyBrown.bead import overlap_pbc, distance_pbc

class Box():

	def __init__(self, beads, box_length, T, viscosity):

		self.beads = beads
		self.box_length = box_length
		self.T = T
		self.viscosity = viscosity

	def propagate(self, dt, build_D = True, cholesky = True, overlaps = True):

		if build_D:
			self.compute_rijmatrix()
			self.compute_Dmatrix()
		if cholesky:
			self.decompose_Dmatrix()
		
		while True:

			BX = self.B @ np.random.normal(0.0, 1.0, 3 * len(self.beads)) * math.sqrt(2 * dt)

			for i, bead in enumerate( self.beads ):
				bead.translate( BX[3 * i: 3 * (i + 1)] )

			if self.check_overlaps():
				for i, bead in enumerate( self.beads ):
					bead.translate( -BX[3 * i: 3 * (i + 1)] )
			else:
				for i, bead in enumerate( self.beads ):
					bead.keep_in_box(self.box_length)
				break

	def check_overlaps(self):

		overlaps = False
		for i in range(len(self.beads)-1):
			for j in range(i+1, len(self.beads)):
				pointer = self.beads[i].r - self.beads[j].r
				radii_sum = self.beads[i].a + self.beads[j].a
				radii_sum_pbc = self.box_length - radii_sum
				if ( pointer[0] > radii_sum and pointer[0] < radii_sum_pbc ) or ( pointer[0] < -radii_sum and pointer[0] > -radii_sum_pbc ):
					continue
				elif ( pointer[1] > radii_sum and pointer[1] < radii_sum_pbc ) or ( pointer[1] < -radii_sum and pointer[1] > -radii_sum_pbc ):
					continue
				elif ( pointer[2] > radii_sum and pointer[2] < radii_sum_pbc ) or ( pointer[2] < -radii_sum and pointer[2] > -radii_sum_pbc ):
					continue
				else:
					if overlap_pbc(self.beads[i], self.beads[j], self.box_length): overlaps = True
		return overlaps

	def compute_rijmatrix(self, nearest = True):

		self.rij = np.zeros((len(self.beads), len(self.beads), 3))

		for i in range(1, len(self.beads)):
			for j in range(0, i):
				self.rij[i][j] = self._pointer(self.beads[i], self.beads[j], nearest)
				self.rij[j][i] = -self.rij[i][j]

	def compute_Dmatrix(self):

		self.D = Boltzmann * self.T * 10**19 * M_rpy(self.beads, self.rij) / self.viscosity

		# self.D = Boltzmann * self.T * 10**19 * M_rpy_smith(self.beads, self.rij) / self.viscosity

	def decompose_Dmatrix(self):

		self.B = np.linalg.cholesky(self.D)

	def _pointer(self, p1, p2, nearest = True):

		point0 = p2.r - p1.r

		dist0 = math.sqrt( point0[0]**2 + point0[1]**2 + point0[2]**2 )

		if nearest and dist0 > self.box_length / 2:

			versor = np.array( [0.0, 0.0, 0.0] )

			for i in range(3):

				if point0[i] > self.box_length / 2:
					versor[i] += self.box_length
				elif point0[i] < -self.box_length / 2:
					versor[i] -= self.box_length

			p1.translate(versor)

			point1 = p2.r - p1.r

			dist1 = math.sqrt( point1[0]**2 + point1[1]**2 + point1[2]**2 )

			p1.translate(-versor)

			if dist1 <= dist0:

				return point1

			else:

				print('ALERT')

				1/0

		return point0

	def __str__(self):

		return 'a'