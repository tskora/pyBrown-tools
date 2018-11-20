import random

def pack_molecules(input_data_object):

    grid = _generate_grid(input_data_object)

    return _populate_grid(grid, input_data_object)

def _generate_grid(input_data_object):

    input_data = input_data_object.input_data

    box_size = input_data["box_size"]
    molecule_radius = input_data["hydrodynamic_radius"]
    min_dist_between_surfaces = input_data["minimal_distance_between_surfaces"]

    unit_dist = 2 * molecule_radius + min_dist_between_surfaces

    num_of_unit_x = int(box_size[0] / unit_dist)
    num_of_unit_y = int(box_size[1] / unit_dist)
    num_of_unit_z = int(box_size[2] / unit_dist)

    return [[(0.5 + i) * unit_dist, (0.5 + j) * unit_dist, (0.5 + k) * unit_dist]
            for i in range(num_of_unit_x)
            for j in range(num_of_unit_y)
            for k in range(num_of_unit_z)]

def _populate_grid(grid, input_data_object):

    input_data = input_data_object.input_data

    number_of_molecules = input_data["number_of_molecules"]

    populated_grid = []

    while len(populated_grid) < number_of_molecules:
        index_to_be_included = random.randint(0, len(grid)-1)
        if grid[index_to_be_included] not in populated_grid:
            populated_grid.append(grid[index_to_be_included])

    return populated_grid
