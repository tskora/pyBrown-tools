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

import numpy as np

class Bead():

	def __init__(self, coords, radius, label):

		self.r = np.array(coords)
		self.a = radius
		self.label = label

	def translate(self, vector):

		self.r += vector

	def keep_in_box(self, box_length):

		for i in range(3):
			while self.r[i] < -box_length / 2:
				self.r[i] += box_length
			while self.r[i] >= box_length / 2:
				self.r[i] -= box_length

	def __str__(self):

		return "{}, radius = {}".format(self.r, self.a)

	def __repr__(self):

		return self.__str__()

	def __eq__(self, p):

		if isinstance( p, Bead ):
			return ( np.all( self.r == p.r ) ) and ( self.a == p.a )
		return False

def distance_pbc(bead1, bead2, box_size):

	pointing = bead1.r - bead2.r	

	for i in range(3):
		while pointing[i] >= box_size/2:
			pointing[i] -= box_size
		while pointing[i] <= -box_size/2:
			pointing[i] += box_size

	return np.sqrt( np.sum( pointing**2 ) )

def overlap_pbc(bead1, bead2, box_size):

	dist = distance_pbc(bead1, bead2, box_size)

	return dist <= bead1.a + bead2.a