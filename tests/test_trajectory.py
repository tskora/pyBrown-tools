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

from pyBrown.trajectories import _compute_sd
from pyBrown.trajectories import _coord_unify

#-------------------------------------------------------------------------------

class TestTrajectory(unittest.TestCase):

	def setUp(self):
		
		return 0

	#---------------------------------------------------------------------------

	def test_compute_sd(self):

		t = np.array([ [ 25.0, 25.0, 25.0 ], [ 26.0, 25.0, 25.0 ], [27.0, 25.0, 25.0] ],
			float)

		s1 = _compute_sd(t, 50.0)
		s2 = _compute_sd(t, 24.0)
		s3 = _compute_sd(t, 25.0)
		s4 = _compute_sd(t, 26.5)

		self.assertSequenceEqual( list( s1 ), [0.0, 1.0, 4.0] )
		self.assertSequenceEqual( list( s2 ), [0.0, 1.0, 4.0] )
		self.assertSequenceEqual( list( s3 ), [0.0, 1.0, 4.0] )
		self.assertSequenceEqual( list( s4 ), [0.0, 1.0, 4.0] )

		t = np.array([ [ 25.0, 24.0, 27.0 ], [ 26.0, 25.0, 26.0 ], [25.0, 24.0, 30.0] ],
			float)

		s1 = _compute_sd(t, 50.0)
		s2 = _compute_sd(t, 24.0)
		s3 = _compute_sd(t, 25.0)
		s4 = _compute_sd(t, 26.5)

		init = t[0]

		ans = [ np.sum( ( t[i] - t[0] )**2 ) for i in range( len(t) ) ]

		self.assertSequenceEqual( list( s1 ), ans )
		self.assertSequenceEqual( list( s2 ), ans )
		self.assertSequenceEqual( list( s3 ), ans )
		self.assertSequenceEqual( list( s4 ), ans )

	#---------------------------------------------------------------------------

	def test_unify_coordinates(self):

		past = [0,0,0]
		present = [1,1,1]
		bs = 100
		_coord_unify( past, present, bs )

		self.assertSequenceEqual( past, [0,0,0] )
		self.assertSequenceEqual( present, [1,1,1] )

		past = [0,0,0]
		present = [1,1,1]
		bs = 0.15
		_coord_unify( past, present, bs )

		self.assertSequenceEqual( past, [0,0,0] )
		answer = [-0.05,-0.05,-0.05]
		for coord1, coord2 in zip(present, answer):
			self.assertAlmostEqual( coord1, coord2, places = 7 )

		past = [0.49, 0.49, 0.49]
		present = [-0.49, -0.49, -0.49]
		bs = 1.0
		_coord_unify( past, present, bs )

		self.assertSequenceEqual( past, [0.49,0.49,0.49] )
		self.assertSequenceEqual( present, [0.51,0.51,0.51] )

		present = [0.49, 0.49, 0.49]
		past = [-0.49, -0.49, -0.49]
		bs = 1.0
		_coord_unify( past, present, bs )

		self.assertSequenceEqual( past, [-0.49,-0.49,-0.49] )
		self.assertSequenceEqual( present, [-0.51,-0.51,-0.51] )


#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
