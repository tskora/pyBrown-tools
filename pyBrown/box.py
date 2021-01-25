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

from pyBrown.hydrodynamics import M_rpy, M_rpy_smith, R_lub_corr

import math
import numpy as np
from scipy.constants import Boltzmann

class Box():

	def __init__(self, beads, box_length, T, viscosity):

		self.beads = beads
		self.box_length = box_length
		self.T = T
		self.viscosity = viscosity

	def propagate(self, dt, build_D = True, cholesky = True):

		if build_D: self.compute_Dmatrix()
		# if cholesky: self.decompose_Dmatrix()
		# BX = self.B @ np.random.normal(0.0, 1.0, 3 * len(self.beads)) * math.sqrt(2 * dt)
		X = np.random.normal(0.0, 1.0, 3 * len(self.beads))
		for i in range(len(self.beads)): X[3*i: 3*(i+1)] /= np.linalg.norm(X[3*i: 3*(i+1)])
		X *= math.sqrt(6 * dt)

		for i, bead in enumerate( self.beads ):
			# bead.translate( BX[3 * i: 3 * (i + 1)] )
			bead.translate( math.sqrt(self.D[0][0]) * X[3*i: 3*(i+1)])
			bead.keep_in_box(self.box_length)

	def compute_Dmatrix(self):

		# self.compute_pointers()

		self.D = np.identity(3*len(self.beads)) * 10**19 * Boltzmann / (6 * np.pi) * self.T / ( self.beads[0].a * self.viscosity )

		# self.D = np.asfortranarray( Boltzmann * self.T * 10**19 * M_rpy(self.beads, self.pointers) / self.viscosity )

		# print(self.D)

		# self.D = np.asfortranarray( Boltzmann * self.T * 10**19 * M_rpy_smith(self.beads, self.box_length, math.sqrt(np.pi), 3, 3) / self.viscosity )

		# print(self.D)

		# self.D = np.asfortranarray( Boltzmann * self.T * 10**19 * np.linalg.inv( ( np.linalg.inv( M_rpy(self.beads, self.pointers) ) + R_lub_corr(self.beads) ) ) / self.viscosity  )

	def decompose_Dmatrix(self):

		self.B = np.linalg.cholesky(self.D)

	def compute_pointers(self, nearest = True):

		self.pointers = [ [ self._pointer(self.beads[i], self.beads[j], nearest) for i in range(len(self.beads)) ] for j in range(len(self.beads)) ]

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