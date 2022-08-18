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

import json

from pyBrown_tools.input import InputData

#-------------------------------------------------------------------------------

class InputDataPores(InputData):

    def __init__(self, input_filename, obligatory_keywords = [], defaults = {}):

        super().__init__(input_filename, obligatory_keywords, defaults)

        # if self.input_data["bond_lengths"] == 'hydrodynamic_radii':
        #     self.bond_lengths_if_only_radii_given()
    
    #---------------------------------------------------------------------------

    # def bond_lengths_if_only_radii_given(self):

    #     bond_lengths = [ ]

    #     for j, bfc in enumerate( self.input_data["bond_force_constants"] ):
    #         bond_lengths.append( [ ] )
    #         for k in range( len( bfc ) ):
    #             bond_lengths[j].append( self.input_data["hydrodynamic_radii"][j][k] +
    #                                     self.input_data["hydrodynamic_radii"][j][k+1] )