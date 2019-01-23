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

from pyBrown.input import InputData
from pyBrown.grid import pack_molecules
from pyBrown.write import write_structure
from pyBrown.parse import parse_input_filename

input_filename = parse_input_filename()
i = InputData(input_filename, [])
# print(i.input_data)
coords = pack_molecules(i)
# print(coords)
write_structure(i, coords)