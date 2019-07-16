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

import random
import numpy as np
from pyBrown.monte_carlo import place_tracers_linearly
from pyBrown.sphere import overlap, fit

from pyBrown.sphere import distance_matrix
from pyBrown.monte_carlo import MonteCarlo

#-------------------------------------------------------------------------------

def pack_molecules(input_data_object):

    input_data = input_data_object.input_data

    if input_data['packing_mode'] == 'regular':

        grid = _generate_grid(input_data_object)

        return _populate_grid(grid, input_data_object)

    elif input_data['packing_mode'] == 'monte_carlo' or \
         input_data['packing_mode'] == 'monte_carlo_fluct':

        return _populate_monte_carlo(input_data_object)

#-------------------------------------------------------------------------------

# TODO: how to use grid for mixtures?

# def _generate_grid(input_data_object):

#     input_data = input_data_object.input_data

#     box_size = input_data["box_size"]
#     molecule_radius = input_data["hydrodynamic_radii"][0][0]
#     min_dist_between_surfaces = input_data["minimal_distance_between_surfaces"]

#     unit_dist = 2 * molecule_radius + min_dist_between_surfaces

#     num_of_unit_x = int(box_size[0] / unit_dist)
#     num_of_unit_y = int(box_size[1] / unit_dist)
#     num_of_unit_z = int(box_size[2] / unit_dist)

#     grid_points = num_of_unit_x * num_of_unit_y * num_of_unit_z

#     print('Number of available grid points is: {}'.format( grid_points ))

#     return [[(0.5 + i) * unit_dist, (0.5 + j) * unit_dist, (0.5 + k) * unit_dist]
#             for i in range(num_of_unit_x)
#             for j in range(num_of_unit_y)
#             for k in range(num_of_unit_z)]

#-------------------------------------------------------------------------------

# TODO: how to use grid for mixtures?

# def _populate_grid(grid, input_data_object):

#     input_data = input_data_object.input_data

#     number_of_molecules = input_data["number_of_molecules"]

#     populated_grid = []

#     assert number_of_molecules <= len( grid ), \
#     'Too many molecules for that box, assuming given minimal distance'

#     while len(populated_grid) < number_of_molecules:
#         index_to_be_included = random.randint(0, len(grid)-1)
#         if grid[index_to_be_included] not in populated_grid:
#             populated_grid.append(grid[index_to_be_included])

#     return populated_grid

#-------------------------------------------------------------------------------

def _populate_monte_carlo(input_data_object):

    input_data = input_data_object.input_data

    numbers_of_molecules = input_data["numbers_of_molecules"]
    box_size = input_data["box_size"]
    radii = input_data["hydrodynamic_radii"]
    min_dist_between_surfaces = input_data["minimal_distance_between_surfaces"]

    populated_box = []

    for i, n_mol in enumerate(numbers_of_molecules):

        thrown = 0

        while thrown < n_mol:

            # works correctly only for cubic boxes
            print('{} / {}'.format(thrown, n_mol))

            tracers = place_tracers_linearly(radii[i], box_size[0])

            print( distance_matrix( tracers ) )
            
            if input_data["packing_mode"] == 'monte_carlo_fluct':

                mc = MonteCarlo( input_data["fluct_length"] )
                
                for tracer in tracers:
                    tracer.translate( mc.get_values() )

            print( distance_matrix( tracers ) )

            if overlap(tracers, populated_box, min_dist_between_surfaces):
                continue

            elif not fit(tracers, box_size[0]):
                continue

            else:
                versors = [ np.array([nx * box_size[0],
                                      ny * box_size[1],
                                      nz * box_size[2]])
                            for nx in np.arange(-1, 2, 1)
                            for ny in np.arange(-1, 2, 1)
                            for nz in np.arange(-1, 2, 1) ]

                if_overlap = False

                for versor in versors:

                    for tracer in tracers:
                        tracer.translate( versor )

                    if overlap(tracers, populated_box, min_dist_between_surfaces):
                        if_overlap = True

                    for tracer in tracers:
                        tracer.translate( -versor )

                if not if_overlap:
                    populated_box += tracers
                    thrown += 1

    return [[populated_box[i].x, populated_box[i].y, populated_box[i].z]
            for i in range( len(populated_box) )]
