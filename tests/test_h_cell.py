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

#-------------------------------------------------------------------------------

class TestSphere(unittest.TestCase):

	def setUp(self):
		
		return 0

	#---------------------------------------------------------------------------

	def test_translate_unit(self):

		return 0		

		self.assertSequenceEqual( list( self.s3.coords ), [1.0, 1.0, 2.0] )

	#---------------------------------------------------------------------------


#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()
