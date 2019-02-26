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

from pyBrown.sphere import Sphere
import numpy as np

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

	def test_rotate_null_center(self):

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

		self.assertSequenceEqual( list( self.s1.coords ), [1.0, 0.0, 0.0] )

	#---------------------------------------------------------------------------

	def test_rotate_unit_90(self):

		self.s1.rotate( [0.0, 0.0, np.pi / 2] )
		answers = [0.0, 1.0, 0.0]

		for coord, answer in zip( self.s1.coords, answers ):
			
			self.assertAlmostEqual( coord, answer, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_180(self):

		self.s1.rotate( [0.0, 0.0, np.pi] )
		answers = [-1.0, 0.0, 0.0]

		for coord, answer in zip( self.s1.coords, answers ):
			
			self.assertAlmostEqual( coord, answer, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_270(self):

		self.s1.rotate( [0.0, 0.0, 3 * np.pi / 2] )
		answers = [0.0, -1.0, 0.0]

		for coord, answer in zip( self.s1.coords, answers ):

			self.assertAlmostEqual( coord, answer, places = 7 )

	#---------------------------------------------------------------------------

	def test_rotate_unit_360(self):

		self.s1.rotate( [0.0, 0.0, 2 * np.pi] )
		answers = [1.0, 0.0, 0.0]

		for coord, answer in zip( self.s1.coords, answers ):

			self.assertAlmostEqual( coord, answer, places = 7 )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
