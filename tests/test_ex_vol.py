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

import sys
sys.path.insert(0, '../pyBrown')
import numpy as np
import copy as cp

from pyBrown.ex_vol import _compute_pore_radius
from pyBrown.sphere import Sphere

#-------------------------------------------------------------------------------

class TestExVol(unittest.TestCase):

	def setUp(self):
		
		self.tracer = Sphere([0.0, 0.0, 0.0], 0.0)

		self.box_size = 10.0

	#---------------------------------------------------------------------------

	def test_in_void(self):

		crowders = [  ]

		r_max = 5.0

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, r_max, 0.001)

		self.assertEqual( pore_radius, r_max )

	#---------------------------------------------------------------------------

	def test_in_cube(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-2,2] for j in [-2,2] for k in [-2,2] ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 2*np.sqrt(3)-1, places=2 )

	#---------------------------------------------------------------------------

	def test_in_cube_2(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-2.3,1.7] for j in [-1.7,2.3] for k in [-2.3,1.7] ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 2*np.sqrt(3)-1, places=2 )

	#---------------------------------------------------------------------------

	def test_in_cube_3(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [0,4] for j in [0,4] for k in [-2,2] ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 1, places=2 )

	#---------------------------------------------------------------------------

	def test_in_cube_fcc(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-2,2] for j in [-2,2] for k in [-2,2] ]

		crowders += [ Sphere([i, j, k], 1.0) for i in [-2,0,2] for j in [-2,0,2] for k in [-2,0,2] if i**2 + j**2 + k**2 == 4 ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 1, places=2 )

	#---------------------------------------------------------------------------

	def test_in_cube_bcc(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-6,-2,2,6] for j in [-6,-2,2,6] for k in [-6,-2,2,6] ]

		crowders += [ Sphere([0.0, 0.0, 0.0], 1.0) ]

		for crowder in crowders: crowder.translate(-np.array([1,1,1]))

		print(crowders)

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, np.sqrt(3)-1, places=2 )

	#---------------------------------------------------------------------------

	def test_in_cube_bcc_2(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-6,-2,2,6] for j in [-6,-2,2,6] for k in [-6,-2,2,6] ]

		crowders += [ Sphere([0.0, 0.0, 0.0], 1.0) ]

		for crowder in crowders: crowder.translate(-np.array([1.1,0.9,0.95]))

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertTrue( pore_radius < 2*np.sqrt(2)-1 )
		self.assertTrue( pore_radius > 1)

	#---------------------------------------------------------------------------

	def test_in_cuboid(self):

		crowders = [ Sphere([i, j, k], 1.0) for i in [-2,2] for j in [-3,3] for k in [-4,4] ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, np.sqrt(29)-1, places=2 )

	#---------------------------------------------------------------------------

	def test_overlapping(self):

		crowders = [ Sphere([0.0, 0.0, 0.0], 1.0) ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 0.0, places=2 )

	#---------------------------------------------------------------------------

	def test_overlapping_pbc(self):

		crowders = [ Sphere([50.0, 50.0, 50.0], 1.0) ]

		pore_radius = _compute_pore_radius(self.tracer, crowders, self.box_size, 50.0, 0.001)

		self.assertAlmostEqual( pore_radius, 0.0, places=2 )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
