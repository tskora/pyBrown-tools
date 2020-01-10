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

from pytest import approx

import sys
sys.path.insert(0, '../pyBrown')
import numpy as np
import copy as cp

from RotDiff import _compute_angular_displacement, _compute_msa

#-------------------------------------------------------------------------------

class TestTrajectory(unittest.TestCase):

	def setUp(self):
		
		return 0

	#---------------------------------------------------------------------------

	def test_compute_angular_displacement(self):

		vec0 = np.array( [0.0, 0.0, 0.0], float )

		vec1 = np.array( [1.0, 0.0, 0.0], float )
		vec2 = np.array( [0.0, 1.0, 0.0], float )
		vec3 = np.array( [0.0, 0.0, 1.0], float )

		self.assertEqual( _compute_angular_displacement( vec1, vec1 ), 0.0 )
		self.assertEqual( _compute_angular_displacement( vec1, -vec1 ), np.pi )

		self.assertEqual( _compute_angular_displacement( vec1, vec2 ), np.pi / 2.0 )
		self.assertEqual( _compute_angular_displacement( vec1, vec3 ), np.pi / 2.0 )
		self.assertEqual( _compute_angular_displacement( -vec1, vec2 ), np.pi / 2.0 )
		self.assertEqual( _compute_angular_displacement( -vec1, vec3 ), np.pi / 2.0 )

		self.assertEqual( _compute_angular_displacement( vec1, 3.0 * vec1 ), 0.0 )
		self.assertEqual( _compute_angular_displacement( vec1, -2.0 * vec2 ), np.pi / 2.0 )
		self.assertEqual( _compute_angular_displacement( vec1, 14.0 * vec3 ), np.pi / 2.0 )

		self.assertEqual( _compute_angular_displacement( vec0, vec1 ), 0.0 )
		self.assertEqual( _compute_angular_displacement( vec0, vec2 ), 0.0 )
		self.assertEqual( _compute_angular_displacement( vec0, vec3 ), 0.0 )

	#---------------------------------------------------------------------------

	def test_compute_msa(self):

		angles = [ [0.0, 0.0, 0.0], [1.0, 1.0, 1.0] ]

		self.assertSequenceEqual( list( _compute_msa(angles, mode = 'window') ), [ 0.0, 0.0, 0.0 ] )
		self.assertSequenceEqual( list( _compute_msa(angles, mode = 'direct') ), [ 0.0, 0.0, 0.0 ] )

		angles2 = [ [ 0.0, np.pi / 2.0, 0.0 ] ]

		self.assertTrue( approx( list( _compute_msa(angles2, 'window') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )
		self.assertTrue( approx( list( _compute_msa(angles2, 'direct') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )

		angles3 = [ [ 0.0, np.pi / 2.0, 0.0 ], [ 0.0, np.pi / 2.0, 0.0 ] ]

		self.assertTrue( approx( list( _compute_msa(angles3, 'window') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )
		self.assertTrue( approx( list( _compute_msa(angles3, 'direct') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )

		angles4 = [ [ 0.0, np.pi / 2.0, 0.0 ], [ np.pi / 2.0, 0.0, np.pi / 2.0 ] ]

		self.assertTrue( approx( list( _compute_msa(angles4, 'window') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )
		self.assertTrue( approx( list( _compute_msa(angles4, 'direct') ), [ 0.0, np.pi**2 / 4, 0.0 ] ) )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
