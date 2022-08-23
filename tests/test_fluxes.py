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

import unittest
import numpy as np

import os.path
import sys
sys.path.insert(0, os.path.abspath( os.path.join(os.path.dirname(__file__), '..') ))

from pyBrown_tools.fluxes import compute_flux_single

#-------------------------------------------------------------------------------

class TestFluxes(unittest.TestCase):

	def setUp(self):

		self.box_size = 200.0
		
		return 0

	#---------------------------------------------------------------------------

	def test_nonmoving_particle(self):

		r0 = np.array([20.0, 10.0, 20.0])
		r1 = r0
		plane_normal_vectors = [ np.array([i, j, k]) for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1] ]
		plane_points = [ r0 ] + [ i for i in plane_normal_vectors ]

		for point in plane_points:
			for vector in plane_normal_vectors:
				self.assertEqual( compute_flux_single(r0, r1, vector, point, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_forward(self):

		r0 = np.array([20.0, 10.0, -20.0])
		r1 = np.array([20.0, 10.0, 20.0])
		plane_normal_vector = np.array([0.0, 0.0, 1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_forward_reversed_normal(self):

		r0 = np.array([20.0, 10.0, -20.0])
		r1 = np.array([20.0, 10.0, 20.0])
		plane_normal_vector = np.array([0.0, 0.0, -1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_backward(self):

		r0 = np.array([20.0, 10.0, 20.0])
		r1 = np.array([20.0, 10.0, -20.0])
		plane_normal_vector = np.array([0.0, 0.0, 1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), -1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_backward_reversed_normal(self):

		r0 = np.array([20.0, 10.0, 20.0])
		r1 = np.array([20.0, 10.0, -20.0])
		plane_normal_vector = np.array([0.0, 0.0, -1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), 1 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_perpendicular(self):

		r0 = np.array([20.0, 10.0, 20.0])
		r1 = np.array([20.0, 10.0, -20.0])
		plane_normal_vector = np.array([-1.0, 1.0, 0.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_through_pbc(self):

		r0 = np.array([20.0, 10.0, -90.0])
		r1 = np.array([20.0, 10.0, 90.0])
		plane_normal_vector = np.array([0.0, 0.0, 1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

	#---------------------------------------------------------------------------

	def test_moving_through_pbc_reversed(self):

		r0 = np.array([20.0, 10.0, 90.0])
		r1 = np.array([20.0, 10.0, -90.0])
		plane_normal_vector = np.array([0.0, 0.0, 1.0])
		plane_point_1 = np.array([0.0, 0.0, -25.0])
		plane_point_2 = np.array([0.0, 0.0, -15.0])
		plane_point_3 = np.array([0.0, 0.0, -5.0])
		plane_point_4 = np.array([0.0, 0.0, 5.0])
		plane_point_5 = np.array([0.0, 0.0, 15.0])
		plane_point_6 = np.array([0.0, 0.0, 25.0])

		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_1, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_2, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_3, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_4, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_5, self.box_size), 0 )
		self.assertEqual( compute_flux_single(r0, r1, plane_normal_vector, plane_point_6, self.box_size), 0 )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()