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

import numpy as np
from pyBrown.sphere import Sphere

#-------------------------------------------------------------------------------

class MonteCarlo:

    def __init__(self, box_size, spherical = True):

        self.spherical = spherical
        self.box_size = box_size

    #---------------------------------------------------------------------------

    def get_values(self):

        if self.spherical: dimension = 3
        else: dimension = 6

        draw = np.random.rand(dimension)
        draw[0:3] = draw[0:3] * self.box_size - self.box_size / 2

        if not self.spherical:
            draw[3] = 2 * np.pi * draw[3] - np.pi
            draw[4] = np.arcsin( 2 * draw[4] - 1 )
            draw[5] = 2 * np.pi * draw[5] # ad hoc

        return draw

    #---------------------------------------------------------------------------

def place_crowders_linearly(radii, box_size):

    if not isinstance(radii, list): radii = [radii]

    crowders = []
    length = 0.0

    for radius in radii:
        length += radius
        crowders.append( Sphere([0.0, 0.0, 0.0 + length], radius) )
        length += radius

    for crowder in crowders:
        crowder.translate([0.0, 0.0, -length / 2])

    return crowders

#-------------------------------------------------------------------------------

def place_crowders_xyz(radii, positions):

    crowders = []

    for radius, position in zip( radii, positions ):
        crowders.append( Sphere(position, radius) )

    return crowders

#-------------------------------------------------------------------------------

def place_tracers_linearly(radii, box_size, bond_lengths = None):

    if not isinstance(radii, list):
        spherical = True
        radii = [radii]
        mc = MonteCarlo(box_size)
    else:
        spherical = False
        mc = MonteCarlo(box_size, spherical)

    tracers = []
    length = 0.0

    if bond_lengths == None:

        for radius in radii:
            length += radius
            tracers.append( Sphere([0.0, 0.0, 0.0 + length], radius) )
            length += radius

        # print(radii)

        # assert length == 2 * np.sum( radii ) - radii[0] - radii[-1]

    else:

        for i, radius in enumerate(radii):
            tracers.append( Sphere([0.0, 0.0, 0.0 + length], radius) )
            if i < len(bond_lengths): length += bond_lengths[i]

        length += radii[0] + radii[-1]

        # assert length == np.sum( bond_lengths ) + radii[0] + radii[-1]

    for tracer in tracers:
            tracer.translate([0.0, 0.0, -length / 2])

    draw = mc.get_values()

    for tracer in tracers:
        if not spherical:
            angles = list(draw[3:5]) + [0.0]
            tracer.rotate(angles)

        tracer.translate(draw[0:3])

    return tracers



