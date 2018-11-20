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


