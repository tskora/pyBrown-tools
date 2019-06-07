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

from pyBrown.h_cell import H_cell

#-------------------------------------------------------------------------------

class TestH_cell(unittest.TestCase):

	def setUp(self):

		self.grid_sizes = [10, 25, 99]

		self.hs = [ H_cell(a = 1.0, b = 1.0, diff_coef = 1.0,
			grid_points = g, v = 1) for g in self.grid_sizes ]

	#---------------------------------------------------------------------------

	def test_ygrid(self):

		for g, h in zip( self.grid_sizes, self.hs ):
			# checking correct y grid initialization
			self.assertEqual( h.grid_points, g )
			# checking correct cs and ys build
			self.assertEqual( len( h.cs ), len( h.ys ) )
			# checking ys comprising the grid
			self.assertEqual( 2 * h.a / ( g + 1 ), h.dy )

		# checking initial concentration profile
		# BEWARE, IT IS NOT FLEXIBLE SO CHANGING SETUP MAY DISRUPT IT
		self.assertSequenceEqual( list( self.hs[0].cs ), [ 0 ]*5 + [ 1 ]*5 )
		self.assertSequenceEqual( list( self.hs[2].cs ), [ 0 ]*49 + [ 1 ]*50 )
		self.assertSequenceEqual( list( self.hs[1].cs ), [ 0 ]*12 + [ 1 ]*13 )

	#---------------------------------------------------------------------------


#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
