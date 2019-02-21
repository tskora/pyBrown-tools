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

from pyBrown.input import InputData
from pyBrown.parse import parse_input_filename
from pyBrown.monte_carlo import MonteCarlo, place_crowders_linearly, place_tracers_linearly
from pyBrown.sphere import Sphere, overlap

import numpy as np

input_filename = parse_input_filename()
i = InputData(input_filename, [])
# print(i.input_data)

box_size = i.input_data["box_size"]
r_crowders = i.input_data["crowders_radii"]
r_tracer = i.input_data["tracers_radii"]
number_of_trials = i.input_data["number_of_draws"]

def sphere_volume(radius):
    return 4 / 3 * np.pi * radius**3

def sphere_surface(radius):
    return 4 * np.pi * radius**2

def cylinder_volume(radius, length):
    return np.pi * radius**2 * length

def excluded_volume_sphere_Nspheres(radius_tracer, radius, N):
    return N * sphere_volume(radius) + sphere_volume(radius_tracer) + N * sphere_surface(radius) * radius_tracer + (N + 1) * cylinder_volume(radius_tracer, 2 * radius)

count = 0

crowders = place_crowders_linearly(r_tracer, box_size)
# print(crowders)
for i in range(number_of_trials):
    tracers = place_tracers_linearly(r_crowders, box_size)
    if overlap(crowders, tracers):
        count += 1
    print( '{}\t{}'.format(i, count / (i + 1) * box_size**3) )


ex_vol = count / number_of_trials * box_size**3

# print(ex_vol)
# print( excluded_volume_sphere_Nspheres(2.289, 1.14, 8) )

### SOME TEST CASES

# s_crowders = [Sphere([0.0, 0.0, -r_crowder], r_crowder), Sphere([0.0, 0.0, r_crowder], r_crowder)]

# for i in range(number_of_trials):

#     with Sphere(mc.get_values(), r_tracer) as s_tracer:
#         if any([overlap(s_crowder, s_tracer) for s_crowder in s_crowders]):
#             count += 1

# ex_vol = count / number_of_trials * box_size**3
# expected = 2 * 4 / 3 * np.pi * (r_tracer + r_crowder)**3 - 4 / 3 * np.pi * r_tracer**3 - 2 * np.pi * r_crowder * r_tracer**2

# print('Computed excluded volume:\t{}'.format(ex_vol))
# print('Expected excluded volume:\t{}'.format(expected))

# ###

# box_size = 10.0
# r_crowder = 2.0
# r_tracer = 2.0
# number_of_trials = 1000000

# mc = MonteCarlo(False, box_size)
# count = 0

# s_crowder = Sphere([0.0, 0.0, 0.0], r_crowder)

# for i in range(number_of_trials):

#     shot = mc.get_values()

#     with Sphere(shot[0:3], r_tracer) as s1, Sphere(shot[0:3], r_tracer) as s2:

#         s_tracers = [s1, s2]
#         s1.translate([0.0, 0.0, r_crowder])
#         s2.translate([0.0, 0.0, -r_crowder])

#         for s in s_tracers:
#             angles = list(shot[3:5]) + [0.0]
#             s.rotate(angles)

#         # if any([overlap(s_crowder, s_tracer) for s_tracer in s_tracers]):
#         #     count += 1

#         if overlap(s_crowder, s_tracers):
#             count += 1

# ex_vol = count / number_of_trials * box_size**3
# expected = 2 * 4 / 3 * np.pi * (r_tracer + r_crowder)**3 - 4 / 3 * np.pi * r_tracer**3 - 2 * np.pi * r_crowder * r_tracer**2

# print('Computed excluded volume:\t{}'.format(ex_vol))
# print('Expected excluded volume:\t{}'.format(expected))