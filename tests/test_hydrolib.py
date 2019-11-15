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

from pyBrown.hydrolib import Particle, XA11

#-------------------------------------------------------------------------------

class TestHydrolib(unittest.TestCase):

	def setup(self):

		# TODO
		return Particle([0.0, 0.0, 0.0], 1.0)

	#---------------------------------------------------------------------------

	def test_jeffrey(self):

		self.assertAlmostEqual( XA11(2.0, 1.0), 0.9954, 3 )

#-------------------------------------------------------------------------------

if __name__ == '__main__':
	unittest.main()