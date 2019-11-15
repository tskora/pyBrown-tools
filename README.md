# pyBrown

`pyBrown` is a bundle of tools useful for **Brownian** and **Stokesian** dynamics
simulations and squbsequent analysis of the results.

Copyright ©2018-2019 Tomasz Skóra [tskora@ichf.edu.pl](mailto:tskora@ichf.edu.pl)

## Features

- [x] computing Mean Squared Displacement
- [x] structure (`.str`) files generation
- [ ] diffusion/mobility/resistance matrix analysis
- [ ] computing Radial Distribution Function
- [ ] computing dynamic geometry features

## Table of Contents

1. [Generating structure files](#strs)
* [Keywords](#strs.keywords)
2. [Trajectory analysis](#traj)
* [Keywords](#traj.keywords)
3. [H Cell simulator](#hcells)

<a name="strs"></a>
## Generating structure files
<a name="strs.keywords"></a>
### Keywords

<a name="traj"></a>
## Trajectory analysis
<a name="traj.keywords"></a>
### Keywords
**Required keywords:**

* `"labels": [string, ...]` -- bead labels in input XYZ file
* `"sizes": [integer, ...]` -- numbers of bead representing individual entities
* `"box_size": float` -- size of simulation (cubic) box (*Å*)
* `"temperature": float` -- temperature (*K*)
* `"viscosity": float` -- dynamic viscosity (*P*)

*pyBrown demands from input `xyz` files a following naming scheme:
..., `TEMPLATE_NUMBER.xyz`, 'TEMPLATE_NUMBER.xyz', ...
(where TEMPLATE is a string variable defined with the keyword `"input_xyz_template"` and NUMBER is an integer from range defined with the keyword `"input_xyz_range"`)*

* `"input_xyz_template": string` -- template of input xyz filenames.
* `"input_xyz_range": [integer, integer]` -- the number range defining input xyz filenames.

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns 1, 2 and 3 (without 4!).)*

* `"debug": boolean` -- (default: `false`)
* `"verbose": boolean` -- (default: `false`)
* `"fit_MSD": boolean` -- (default: `false`)
* `"probing_frequency: integer"` -- (default: `1`)
* `"min_time: float"` -- (default: `0.0`)
* `"mode": option` -- (options: `direct`/`window`, default: `window`)

**Optional keywords:**

*pyBrown demands from input `enr` files a following naming scheme:
..., `TEMPLATE_NUMBER.enr`, 'TEMPLATE_NUMBER.enr', ...
(where TEMPLATE is a string variable defined with the keyword `"input_enr_template"` and NUMBER is an integer from range defined with the keyword `"input_enr_range"`)*

* `"input_enr_template": string` -- template of input enr filenames.
* `"input_enr_range": [integer, integer]` -- the number range defining input enr filenames.

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns 1, 2 and 3 (without 4!).)*

<a name="hcells"></a>
## H Cell simulator

`H_cell_sim.py` enables one to compute the concentration profiles emerging in the H-cell experiments.

How to use:

1. `make`
2. `make test`
3. `python H_cell_sim.py --help`
4. `python H_cell_sim.py -a 1 -b 1 -d "1 1.3" -n 100 -x 0.000001 -m 1 -o hcellsim -s "0.001 0.35"`
