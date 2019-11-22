# pyBrown is a bundle of tools useful for Brownian and Stokesian dynamics
# simulations. Copyright (C) 2018  Tomasz Skora (tskora@ichf.edu.pl)
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

from pyBrown.input import InputData

#-------------------------------------------------------------------------------

class InputDataPack(InputData):

    def __init__(self, input_filename, obligatory_keywords = [], defaults = {}):

        super().__init__(input_filename, obligatory_keywords, defaults)

        if not isinstance( self.input_data["max_bond_lengths"], list ):
            self.max_bond_lengths_if_only_one_given()

        if self.input_data["bond_lengths"] == 'hydrodynamic_radii':
            self.bond_lengths_if_only_radii_given()
    
    #---------------------------------------------------------------------------

    def max_bond_lengths_if_only_one_given(self):

        max_bond_lengths = [ ]
    
        for bfc in self.input_data["bond_force_constants"]:
            max_bond_lengths.append( [ self.input_data["max_bond_lengths"]
                                       for j in range( len( bfc ) ) ] )
        
        self.input_data["max_bond_lengths"] = max_bond_lengths

    #---------------------------------------------------------------------------

    def bond_lengths_if_only_radii_given(self):

        bond_lengths = [ ]

        for j, bfc in enumerate( self.input_data["bond_force_constants"] ):
            bond_lengths.append( [ ] )
            for k in range( len( bfc ) ):
                bond_lengths[j].append( self.input_data["hydrodynamic_radii"][j][k] +
                                        self.input_data["hydrodynamic_radii"][j][k+1] )

        self.input_data["bond_lengths"] = bond_lengths
