import random

def generate_grid(box_size, molecule_radius, min_dist_between_surfaces):

    unit_dist = 2 * molecule_radius + min_dist_between_surfaces

    num_of_unit_x = int(box_size[0] / unit_dist)
    num_of_unit_y = int(box_size[1] / unit_dist)
    num_of_unit_z = int(box_size[2] / unit_dist)

    return [[(0.5 + i) * unit_dist, (0.5 + j) * unit_dist, (0.5 + k) * unit_dist]
            for i in range(num_of_unit_x)
            for j in range(num_of_unit_y)
            for k in range(num_of_unit_z)]

def populate_grid(grid, number_of_molecules):

    populated_grid = []

    while len(populated_grid) < number_of_molecules:
        index_to_be_included = random.randint(0, len(grid)-1)
        if grid[index_to_be_included] not in populated_grid:
            populated_grid.append(grid[index_to_be_included])

    return populated_grid



box_x = 500.0
box_y = 500.0
box_z = 500.0

box_size = [box_x, box_y, box_z]

molecule_radius = 5.1
number_of_molecules = 2
min_dist_surf = 15.0

grid = generate_grid(box_size, molecule_radius, min_dist_surf)
print(grid)
pop_grid = populate_grid(grid, number_of_molecules)
print(pop_grid)
