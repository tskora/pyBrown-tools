# pyVrown is a bound of tools useful for Brownian and Stokesian dynamics simultions
# Copyright (C) 2018  Tomasz Skóra (tskora@ichf.edu.pl)
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

def write_structure(input_data_object, coordinates):

    input_data = input_data_object.input_data

    output_str_filename = input_data["output_structure_filename"]

    with open(output_str_filename, "w") as write_file:
        for i, molecule in enumerate(coordinates):
            line = _line_pattern(input_data, i, molecule)
            write_file.write(line+"\n")

def _line_pattern(input_data, index, coords):

    return "sub ATM {} {:.1f} {:.1f} {:.1f} {} {} {} {} {}".format(
        index+1,
        coords[0],
        coords[1],
        coords[2],
        input_data["hydrodynamic_radius"],
        input_data["charge"],
        2*input_data["lennard-jones_radius"],
        input_data["lennard-jones_energy"],
        input_data["mass"]
    )


