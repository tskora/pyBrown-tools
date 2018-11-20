import argparse

def parse_input_filename():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_filename")
    args = parser.parse_args()
    return args.input_filename