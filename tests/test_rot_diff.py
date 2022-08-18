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
import copy as cp

from pytest import approx

from pyBrown_tools.trajectories import _compute_sad, _compute_autocorrelation

#-------------------------------------------------------------------------------

class TestTrajectory(unittest.TestCase):

	def setUp(self):
		
		return 0

	#---------------------------------------------------------------------------

	def test_compute_sad(self):

		vec0 = np.array( [0.0, 0.0, 0.0], float )

		vec1 = np.array( [1.0, 0.0, 0.0], float )
		vec2 = np.array( [0.0, 1.0, 0.0], float )
		vec3 = np.array( [0.0, 0.0, 1.0], float )

		self.assertSequenceEqual( list( _compute_sad( [vec1, vec1], 'direct' ) ), [0.0, 0.0] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, -vec1], 'direct' ) ), [0.0, np.pi**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, vec2], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, vec3], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec2, vec3], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [-vec1, vec2 ], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [-vec1, vec3 ], 'direct' ) ), [0.0, (np.pi/2)**2] )

		self.assertSequenceEqual( list( _compute_sad( [vec1, vec1], 'window' ) ), [0.0, 0.0] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, -vec1], 'window' ) ), [0.0, np.pi**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, vec2], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, vec3], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec2, vec3], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [-vec1, vec2 ], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [-vec1, vec3 ], 'window' ) ), [0.0, (np.pi/2)**2] )

		self.assertSequenceEqual( list( _compute_sad( [vec1, 3*vec1], 'direct' ) ), [0.0, 0.0] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, -2*vec2], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, 14*vec3], 'direct' ) ), [0.0, (np.pi/2)**2] )

		self.assertSequenceEqual( list( _compute_sad( [vec1, 3*vec1], 'window' ) ), [0.0, 0.0] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, -2*vec2], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec1, 14*vec3], 'window' ) ), [0.0, (np.pi/2)**2] )

		# TODO: what to do if n = (0,0,0)
		self.assertSequenceEqual( list( _compute_sad( [vec0, vec1], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec0, vec2], 'direct' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec0, vec3], 'direct' ) ), [0.0, (np.pi/2)**2] )

		self.assertSequenceEqual( list( _compute_sad( [vec0, vec1], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec0, vec2], 'window' ) ), [0.0, (np.pi/2)**2] )
		self.assertSequenceEqual( list( _compute_sad( [vec0, vec3], 'window' ) ), [0.0, (np.pi/2)**2] )

	#---------------------------------------------------------------------------
	
	def test_compute_autocorrelation(self):

		vec0 = np.array( [0.0, 0.0, 0.0], float )

		vec1 = np.array( [1.0, 0.0, 0.0], float )
		vec2 = np.array( [0.0, 1.0, 0.0], float )
		vec3 = np.array( [0.0, 0.0, 1.0], float )

		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, vec1], 'direct' ) ), [1.0, 1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec2, vec2], 'direct' ) ), [1.0, 1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec3, vec3], 'direct' ) ), [1.0, 1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, -vec1], 'direct' ) ), [1.0, -1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec2, -vec2], 'direct' ) ), [1.0, -1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec3, -vec3], 'direct' ) ), [1.0, -1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, vec2], 'direct' ) ), [1.0, 0.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec2, vec3], 'direct' ) ), [1.0, 0.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec3, vec1], 'direct' ) ), [1.0, 0.0] )

		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, 3*vec1], 'direct' ) ), [1.0, 1.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, -2*vec2], 'direct' ) ), [1.0, 0.0] )
		self.assertSequenceEqual( list( _compute_autocorrelation( [vec1, 14*vec3], 'direct' ) ), [1.0, 0.0] )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
