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

import json

class InputData:

    def __init__(self, input_filename, keywords):

        self.read_input_file(input_filename)

        self.check_for_missing_keywords(keywords)

    def read_input_file(self, input_filename):

        with open(input_filename, "r") as read_file:
            self.input_data = json.load(read_file)

    def check_for_missing_keywords(self, keywords):

        for keyword in keywords:
            assert keyword in self.input_data.keys(),\
                'Missing {} keyword in input JSON file.'.format(keyword)


