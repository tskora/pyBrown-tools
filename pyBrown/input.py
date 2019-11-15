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

#-------------------------------------------------------------------------------

class InputData:

    def __init__(self, input_filename, obligatory_keywords = [], defaults = {}):

        self._read_input_file(input_filename)

        self._complete_with_defaults(defaults)

        self._check_for_missing_keywords(obligatory_keywords)

    #---------------------------------------------------------------------------

    def __str__(self):

        string_representation = ''

        for keyword in self.input_data.keys():
            string_representation += '{}: {}\n'.format(keyword, self.input_data[keyword])

        return string_representation

    #---------------------------------------------------------------------------

    def __repr__(self):

        return self.__str__()

    #---------------------------------------------------------------------------

    def _read_input_file(self, input_filename):

        with open(input_filename, "r") as read_file:
            self.input_data = json.load(read_file)

    #---------------------------------------------------------------------------

    def _complete_with_defaults(self, defaults):

        for default_keyword in defaults.keys():

            if default_keyword not in self.input_data.keys():

                self.input_data[default_keyword] = defaults[default_keyword]

    #---------------------------------------------------------------------------

    def _check_for_missing_keywords(self, obligatory_keywords):

        for keyword in obligatory_keywords:
            assert keyword in self.input_data.keys(),\
                'Missing {} keyword in input JSON file.'.format(keyword)