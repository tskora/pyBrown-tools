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

import scipy.interpolate as interpolate
from scipy.integrate import quad

from pyBrown.monte_carlo import place_tracers_linearly
from pyBrown.sphere import overlap, fit

from pyBrown.sphere import distance_matrix
from pyBrown.monte_carlo import MonteCarlo

#-------------------------------------------------------------------------------

def pack_molecules(input_data):

    if input_data['packing_mode'] == 'regular':

        # that option does not work

        grid = _generate_grid(input_data)

        return _populate_grid(grid, input_data)

    elif input_data['packing_mode'] == 'monte_carlo':

        return _populate_monte_carlo(input_data)

    elif input_data['packing_mode'] == 'monte_carlo_fluct':

        return _populate_monte_carlo_fluct(input_data)

    else:

        return None

#-------------------------------------------------------------------------------

# TODO: how to use grid for mixtures?
def _generate_grid(input_data):

    return None

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
def _populate_grid(grid, input_data):

    return None

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

# TODO: it should work also for noncubic boxes
# TODO: bond lengths not utilised
def _populate_monte_carlo(input_data):

    numbers_of_molecules = input_data["numbers_of_molecules"]
    box_size = input_data["box_size"]

    radii = input_data["hydrodynamic_radii"]
    # bond_lengths = input_data["bond_lengths"]

    min_dist_between_surfaces = input_data["minimal_distance_between_surfaces"]

    populated_box = []

    for i, n_mol in enumerate(numbers_of_molecules):

        thrown = 0

        while thrown < n_mol:

            print('{} / {}'.format(thrown, n_mol))

            tracers = place_tracers_linearly(radii[i], box_size[0])

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

#-------------------------------------------------------------------------------

# TODO: ad hoc solution -- drawing predefined amount of random numbers
# TODO: bond lengths and open/close radii are redundant info, bond lengths should be removed
def _populate_monte_carlo_fluct(input_data, n_buff = 100000):

    numbers_of_molecules = input_data["numbers_of_molecules"]
    box_size = input_data["box_size"]
    temperature = input_data["temperature"]

    radii = input_data["hydrodynamic_radii"]
    bond_potential = input_data["bond_potential"]
    bond_force_constants = input_data["bond_force_constants"]
    open_radii = input_data["open_radii"]
    close_radii = input_data["close_radii"]

    min_dist_between_surfaces = input_data["minimal_distance_between_surfaces"]

    # ad hoc solution
    random_bond_lengths = [ draw_bond_lengths(open_radii[n], close_radii[n], bond_force_constants[n], bond_potential, temperature, n_samples = n_buff) for n in range( len( numbers_of_molecules) ) ]

    print( random_bond_lengths )

    populated_box = []

    for i, n_mol in enumerate(numbers_of_molecules):

        thrown = 0

        while thrown < n_mol:

            # works correctly only for cubic boxes
            print('{} / {}'.format(thrown, n_mol))

            rbl = [ random_bond_lengths[i][counter][thrown] for counter in range( len(bond_force_constants[i]) )]

            print( 'random bond lengths: {}'.format(random_bond_lengths) )

            tracers = place_tracers_linearly(radii[i], box_size[0], rbl)

            print('coordinates: {}'.format(tracers) )

            print('distance matrix: {}'.format( distance_matrix(tracers) ) )

            if overlap(tracers, populated_box, min_dist_between_surfaces):
                print('overlap')
                continue

            elif not fit(tracers, box_size[0]):
                print('not fit')
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

#-------------------------------------------------------------------------------

def draw_bond_lengths(open_radii, close_radii, bond_force_constants,
                      bond_potential, temperature, n_samples = 1):

    random_bond_lengths = [ ]

    for n, bond_force_constant in enumerate( bond_force_constants ):

        open_length = open_radii[n] + open_radii[n + 1]
        close_length = close_radii[n] + close_radii[n + 1]

        random_bond_lengths.append(
            _draw_bond_length( open_length, close_length, bond_force_constant,
                               bond_potential, temperature, n_samples )
            )

    return random_bond_lengths

#-------------------------------------------------------------------------------

def _double_potential(length, open_length, close_length, bond_force_constant):

    mult = 16 * bond_force_constant / ( open_length - close_length )**4

    return mult * ( length - open_length )**2 * ( length - close_length )**2

#-------------------------------------------------------------------------------

def _close_potential(length, open_length, close_length, bond_force_constant):

    mult = 16 * bond_force_constant / ( open_length - close_length )**4 * \
           ( length - close_length )**2

    if length < close_length:
        return mult * ( length - open_length )**2
    else:
        return mult * ( length - 2 * close_length + open_length )**2

#-------------------------------------------------------------------------------

def _open_potential(length, open_length, close_length, bond_force_constant):

    mult = 16 * bond_force_constant / ( open_length - close_length )**4 * \
           ( length - open_length )**2

    if length > open_length:
        return mult * ( length - close_length )**2
    else:
        return mult * ( length - 2 * open_length + close_length )**2

#-------------------------------------------------------------------------------

# TODO: +10 and -10 as tresholds in lengths is an ad-hoc solution
def _draw_bond_length(open_length, close_length, bond_force_constant,
                      bond_potential, temperature, n_samples = 1):

    lengths = np.linspace(close_length - np.min([10.0, close_length]), open_length + 10.0, 10000)

    kb_T = ( temperature / 298.15 ) * 0.5924812013256604

    if bond_potential == 'double': potential = _double_potential
    elif bond_potential == 'open': potential = _open_potential
    elif bond_potential == 'close': potential = _close_potential

    boltz_factor = lambda l: l**2 * \
     np.exp( -potential(l, open_length, close_length, bond_force_constant) / kb_T )

    I= []
    norm, _ = quad(boltz_factor, lengths[0], lengths[-1])
    for i in range( len( lengths ) ):
        integral, _ = quad(boltz_factor, lengths[0], lengths[i])
        I.append( integral / norm )

    inv_cdf = interpolate.interp1d(I, lengths)

    r = np.random.rand(n_samples)
    return inv_cdf(r)