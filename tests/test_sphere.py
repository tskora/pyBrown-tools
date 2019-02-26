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

from pyBrown.sphere import Sphere, _overlap

#-------------------------------------------------------------------------------

class TestSphere(unittest.TestCase):

	def setUp(self):
		
		self.s = Sphere([0.0, 0.0, 0.0], 1.0)

		self.s1 = Sphere([1.0, 0.0, 0.0], 1.0)

		self.s2 = Sphere([0.0, 1.0, 0.0], 1.0)

		self.s3 = Sphere([0.0, 0.0, 1.0], 1.0)

	#---------------------------------------------------------------------------

	def test_translate_unit(self):

		vector_of_translation = [1.0, 1.0, 1.0]

		self.s.translate( vector_of_translation )

		self.s1.translate( vector_of_translation )

		self.s2.translate( vector_of_translation )

		self.s3.translate( vector_of_translation )

		self.assertSequenceEqual( list( self.s.coords ), [1.0, 1.0, 1.0] )

		self.assertSequenceEqual( list( self.s1.coords ), [2.0, 1.0, 1.0] )

		self.assertSequenceEqual( list( self.s2.coords ), [1.0, 2.0, 1.0] )

		self.assertSequenceEqual( list( self.s3.coords ), [1.0, 1.0, 2.0] )

	#---------------------------------------------------------------------------

	def test_translate_null(self):

		self.s.translate( [0.0, 0.0, 0.0] )

		self.s1.translate( [0.0, 0.0, 0.0] )

		self.s2.translate( [0.0, 0.0, 0.0] )

		self.s3.translate( [0.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s.coords ), [0.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s1.coords ), [1.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s2.coords ), [0.0, 1.0, 0.0] )

		self.assertSequenceEqual( list( self.s3.coords ), [0.0, 0.0, 1.0] )

	#---------------------------------------------------------------------------

	def test_rotate_null(self):

		self.s.rotate( [0.0, 0.0, 0.0] )

		self.s1.rotate( [0.0, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s.coords ), [0.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s1.coords ), [1.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s2.coords ), [0.0, 1.0, 0.0] )

		self.assertSequenceEqual( list( self.s3.coords ), [0.0, 0.0, 1.0] )

	#---------------------------------------------------------------------------

	def test_rotate_unit_0(self):

		self.s1.rotate( [0.0, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s1.coords ), [1.0, 0.0, 0.0] )

		self.assertSequenceEqual( list( self.s2.coords ), [0.0, 1.0, 0.0] )

		self.assertSequenceEqual( list( self.s3.coords ), [0.0, 0.0, 1.0] )

	#---------------------------------------------------------------------------

	def test_rotate_unit_90(self):

		cases = []
		answers = []

		self.s1.rotate( [0.0, 0.0, np.pi / 2] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s1.rotate( [0.0, 0.0, -np.pi / 2] )

		self.s1.rotate( [0.0, np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [0.0, 0.0, -1.0] )
		self.s1.rotate( [0.0, -np.pi / 2, 0.0] )

		self.s1.rotate( [np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [-np.pi / 2, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, np.pi / 2] )
		cases.append( cp.copy( self.s2) )
		answers.append( [-1.0, 0.0, 0.0] )
		self.s2.rotate( [0.0, 0.0, -np.pi / 2] )

		self.s2.rotate( [0.0, np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [0.0, -np.pi / 2, 0.0] )

		self.s2.rotate( [np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s2.rotate( [np.pi / 2, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, np.pi / 2] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [0.0, 0.0, -np.pi / 2] )

		self.s3.rotate( [0.0, np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s3.rotate( [0.0, -np.pi / 2, 0.0] )
		
		self.s3.rotate( [np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, -1.0, 0.0] )
		self.s3.rotate( [-np.pi / 2, 0.0, 0.0] )

		for case, answer in zip( cases, answers ):

			for coord, answer_coord in zip( case.coords, answer ):

				self.assertAlmostEqual( coord, answer_coord, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_180(self):

		# IN PROGRESS

		cases = []
		answers = []

		self.s1.rotate( [0.0, 0.0, np.pi] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [-1.0, 0.0, 0.0] )
		self.s1.rotate( [0.0, 0.0, -np.pi] )

		self.s1.rotate( [0.0, np.pi, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [-1.0, 0.0, 0.0] )
		self.s1.rotate( [0.0, -np.pi, 0.0] )

		self.s1.rotate( [np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [-np.pi, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, np.pi] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, -1.0, 0.0] )
		self.s2.rotate( [0.0, 0.0, -np.pi] )

		self.s2.rotate( [0.0, np.pi, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [0.0, -np.pi, 0.0] )

		self.s2.rotate( [np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, -1.0, 0.0] )
		self.s2.rotate( [np.pi, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, np.pi] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [0.0, 0.0, -np.pi] )

		self.s3.rotate( [0.0, np.pi, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, -1.0] )
		self.s3.rotate( [0.0, -np.pi, 0.0] )
		
		self.s3.rotate( [np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, -1.0] )
		self.s3.rotate( [-np.pi, 0.0, 0.0] )

		for case, answer in zip( cases, answers ):

			for coord, answer_coord in zip( case.coords, answer ):

				self.assertAlmostEqual( coord, answer_coord, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_270(self):

		cases = []
		answers = []

		self.s1.rotate( [0.0, 0.0, 3 * np.pi / 2] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [0.0, -1.0, 0.0] )
		self.s1.rotate( [0.0, 0.0, -3 * np.pi / 2] )

		self.s1.rotate( [0.0, 3 * np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s1.rotate( [0.0, -3 * np.pi / 2, 0.0] )

		self.s1.rotate( [3 * np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [-3 * np.pi / 2, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, 3 * np.pi / 2] )
		cases.append( cp.copy( self.s2) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s2.rotate( [0.0, 0.0, -3 * np.pi / 2] )

		self.s2.rotate( [0.0, 3 * np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [0.0, -3 * np.pi / 2, 0.0] )

		self.s2.rotate( [3 * np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 0.0, -1.0] )
		self.s2.rotate( [-3 * np.pi / 2, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, 3 * np.pi / 2] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [0.0, 0.0, -3 * np.pi / 2] )

		self.s3.rotate( [0.0, 3 * np.pi / 2, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [-1.0, 0.0, 0.0] )
		self.s3.rotate( [0.0, -3 * np.pi / 2, 0.0] )
		
		self.s3.rotate( [3 * np.pi / 2, 0.0, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s3.rotate( [-3 * np.pi / 2, 0.0, 0.0] )

		for case, answer in zip( cases, answers ):

			for coord, answer_coord in zip( case.coords, answer ):

				self.assertAlmostEqual( coord, answer_coord, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_360(self):

		cases = []
		answers = []

		self.s1.rotate( [0.0, 0.0, 2 * np.pi] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [0.0, 0.0, -2 * np.pi] )

		self.s1.rotate( [0.0, 2 * np.pi, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [0.0, -2 * np.pi, 0.0] )

		self.s1.rotate( [2 * np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s1 ) )
		answers.append( [1.0, 0.0, 0.0] )
		self.s1.rotate( [-2 * np.pi, 0.0, 0.0] )

		self.s2.rotate( [0.0, 0.0, 2 * np.pi] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [0.0, 0.0, -2 * np.pi] )

		self.s2.rotate( [0.0, 2 * np.pi, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [0.0, -2 * np.pi, 0.0] )

		self.s2.rotate( [2 * np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s2) )
		answers.append( [0.0, 1.0, 0.0] )
		self.s2.rotate( [-2 * np.pi, 0.0, 0.0] )

		self.s3.rotate( [0.0, 0.0, 2 * np.pi] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [0.0, 0.0, -2 * np.pi] )

		self.s3.rotate( [0.0, 2 * np.pi, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [0.0, -2 * np.pi, 0.0] )
		
		self.s3.rotate( [2 * np.pi, 0.0, 0.0] )
		cases.append( cp.copy( self.s3 ) )
		answers.append( [0.0, 0.0, 1.0] )
		self.s3.rotate( [-2 * np.pi, 0.0, 0.0] )

		for case, answer in zip( cases, answers ):

			for coord, answer_coord in zip( case.coords, answer ):

				self.assertAlmostEqual( coord, answer_coord, places = 7 )

	#---------------------------------------------------------------------------

	def test_overlap(self):

		self.assertTrue(_overlap(self.s, self.s1, 0.0))
		self.assertTrue(_overlap(self.s, self.s2, 0.0))
		self.assertTrue(_overlap(self.s, self.s3, 0.0))
		self.assertTrue(_overlap(self.s1, self.s2, 0.0))
		self.assertTrue(_overlap(self.s1, self.s3, 0.0))
		self.assertTrue(_overlap(self.s2, self.s3, 0.0))

		self.s1.translate([1.0, 0.0, 0.0])
		self.s2.translate([0.0, 1.0, 0.0])
		self.s3.translate([0.0, 0.0, 1.0])

		self.assertFalse(_overlap(self.s, self.s1, 0.0))
		self.assertFalse(_overlap(self.s, self.s2, 0.0))
		self.assertFalse(_overlap(self.s, self.s3, 0.0))
		self.assertFalse(_overlap(self.s1, self.s2, 0.0))
		self.assertFalse(_overlap(self.s1, self.s3, 0.0))
		self.assertFalse(_overlap(self.s2, self.s3, 0.0))

		self.s.r = 2.0

		self.assertTrue(_overlap(self.s, self.s1, 0.0))
		self.assertTrue(_overlap(self.s, self.s2, 0.0))
		self.assertTrue(_overlap(self.s, self.s3, 0.0))

		self.s_copy = cp.copy( self.s )
		self.s.r = 0.0
		self.assertTrue(_overlap(self.s, self.s_copy, 0.0))

	def test_overlap_lists(self):

		# TODO: further testing

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
