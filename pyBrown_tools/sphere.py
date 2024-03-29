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

#-------------------------------------------------------------------------------

class Sphere:

    def __init__(self, coords, radius):
        self.coords = np.array(coords, float)
        self.x, self.y, self.z = self.coords[0], self.coords[1], self.coords[2]
        self.r = radius

    #---------------------------------------------------------------------------

    def __str__(self):
        return "{} r = {}".format(self.coords, self.r)

    #---------------------------------------------------------------------------

    def __repr__(self):
        return self.__str__()

    #---------------------------------------------------------------------------

    def __enter__(self):
        return self

    #---------------------------------------------------------------------------

    def __exit__(self, type, value, traceback):
        del self

    #---------------------------------------------------------------------------

    def translate(self, vector):
        vector = np.array(vector, float)
        self.coords += vector
        self.x, self.y, self.z = self.coords[0], self.coords[1], self.coords[2]

    #---------------------------------------------------------------------------

    def rotate(self, angles):
        vector = np.array(angles, float)
        rotation1 = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(angles[0]), -np.sin(angles[0])], [0.0, np.sin(angles[0]), np.cos(angles[0])]])
        rotation2 = np.array([[np.cos(angles[1]), 0.0, np.sin(angles[1])], [0.0, 1.0, 0.0], [-np.sin(angles[1]), 0.0, np.cos(angles[1])]])
        rotation3 = np.array([[np.cos(angles[2]), -np.sin(angles[2]), 0.0], [np.sin(angles[2]), np.cos(angles[2]), 0.0], [0.0, 0.0, 1.0]])
        self.coords = rotation1 @ rotation2 @ rotation3 @ self.coords
        self.x, self.y, self.z = self.coords[0], self.coords[1], self.coords[2]

    #---------------------------------------------------------------------------

    def volume(self):
        return 4.0 * np.pi * self.r**3 / 3.0

    #---------------------------------------------------------------------------

    def surface(self):
        return 4.0 * np.pi * self.r**2

#-------------------------------------------------------------------------------

def _overlap(sphere1, sphere2, minimal_distance):

    return (sphere1.x-sphere2.x)**2 + (sphere1.y-sphere2.y)**2 + (sphere1.z-sphere2.z)**2 < ( sphere1.r + sphere2.r + minimal_distance )**2

#-------------------------------------------------------------------------------

def overlap(sphere1, sphere2, minimal_distance):
    if not isinstance(sphere1, list): sphere1 = [sphere1]
    if not isinstance(sphere2, list): sphere2 = [sphere2]

    return any([_overlap(s1, s2, minimal_distance)
                for s1 in sphere1 for s2 in sphere2])

#-------------------------------------------------------------------------------

def overlap_pbc(sphere1, sphere2, minimal_distance, box_size):
    if not isinstance(sphere1, list): sphere1 = [sphere1]
    if not isinstance(sphere2, list): sphere2 = [sphere2]

    for s1 in sphere1:
        for s2 in sphere2:

            pointing = np.array( [s1.x - s2.x,
                                  s1.y - s2.y,
                                  s1.z - s2.z] )

            radii_sum = s1.r + s2.r + minimal_distance
            radii_sum_pbc = box_size - radii_sum

            if ( pointing[0] > radii_sum and pointing[0] < radii_sum_pbc ) or ( pointing[0] < -radii_sum and pointing[0] > -radii_sum_pbc ):
                continue
            elif ( pointing[1] > radii_sum and pointing[1] < radii_sum_pbc ) or ( pointing[1] < -radii_sum and pointing[1] > -radii_sum_pbc ):
                continue
            elif ( pointing[2] > radii_sum and pointing[2] < radii_sum_pbc ) or ( pointing[2] < -radii_sum and pointing[2] > -radii_sum_pbc ):
                continue
            else:
                for i in range(3):
                    while pointing[i] >= box_size/2:
                        pointing[i] -= box_size
                    while pointing[i] <= -box_size/2:
                        pointing[i] += box_size

                if np.sum( pointing**2 ) <= radii_sum**2 + minimal_distance: return True

    return False

#-------------------------------------------------------------------------------

def _distance(sphere1, sphere2):

    return np.sqrt( (sphere1.x - sphere2.x)**2 +
                    (sphere1.y - sphere2.y)**2 +
                    (sphere1.z - sphere2.z)**2 )

#-------------------------------------------------------------------------------

def distance_matrix(spheres):

    return np.array([ [ _distance(s1, s2) for s2 in spheres ]
                      for s1 in spheres ], float)

#-------------------------------------------------------------------------------

def _fit(sphere, box_size):

    return all( [ ( sphere.coords[i] < box_size / 2 and
                    sphere.coords[i] > -box_size / 2 ) for i in range(3) ] )

#-------------------------------------------------------------------------------

def fit(spheres, box_size):

    return all([_fit(s, box_size) for s in spheres])
