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


