# pyBrown

[![Build Status](https://travis-ci.com/tskora/pyBrown.svg?branch=develop)](https://travis-ci.com/tskora/pyBrown)

`pyBrown` is a bundle of tools useful for **Brownian** and **Stokesian** dynamics
simulations and squbsequent analysis of the results.

Copyright &copy;2018-2019 Tomasz Skóra [tskora@ichf.edu.pl](mailto:tskora@ichf.edu.pl)

## Features

- [x] computing Mean Squared Displacement
- [x] computing Radial Distribution Function
- [x] structure (`.str`) files generation
- [ ] diffusion/mobility/resistance matrix analysis

## First steps

Type those commands in a terminal:

`make`

`make test`

## Table of Contents

1. [Generating structure files](#strs)
* [Keywords](#strs.keywords)
2. [Trajectory analysis](#traj)
* [Keywords](#traj.keywords)
* [Example input file](#traj.example)
* [Usage](#traj.usage)
* [Output files](#traj.output)
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

* `"labels": [string, ...]` &mdash; bead labels in input XYZ file
* `"sizes": [integer, ...]` &mdash; numbers of bead representing individual entities
* `"box_size": float` &mdash; size of simulation (cubic) box (*Å*)
* `"temperature": float` &mdash; temperature (*K*)
* `"viscosity": float` &mdash; dynamic viscosity (*P*)

*pyBrown demands from input `xyz` files a following naming scheme:
..., `(TEMPLATE)(NUMBER).xyz`, `(TEMPLATE)(NUMBER).xyz`, ...
(where TEMPLATE is a string variable defined with the keyword `"input_xyz_template"` and NUMBER is an integer from range defined with the keyword `"input_xyz_range"`)*

* `"input_xyz_template": string` &mdash; template of input xyz filenames.
* `"input_xyz_range": [integer, integer]` &mdash; the number range defining input xyz filenames.

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns 1, 2 and 3 (without 4!).)*

* `"debug": boolean` &mdash; print extra information useful for debugging purposes (default: `false`)
* `"verbose": boolean` &mdash; print extra information (default: `false`)
* `"fit_MSD": boolean` &mdash; fit linear functions to the computed MSDs and plot them in the output figure (default: `false`)
* `"probing_frequency: integer"` &mdash; read every *N*-th geometry (default: `1`)
* `"min_time: float"` &mdash; not include snapshots with time smaller than `"min_time"` (default: `0.0`)
* `"mode": option` &mdash; learn more [here](https://freud.readthedocs.io/en/v2.0.1/modules/msd.html) (options: `direct`/`window`, default: `window`)

**Optional keywords:**

*pyBrown demands from input `enr` files a following naming scheme:
..., `(TEMPLATE)(NUMBER).enr`, `(TEMPLATE)(NUMBER).enr`, ...
(where TEMPLATE is a string variable defined with the keyword `"input_enr_template"` and NUMBER is an integer from range defined with the keyword `"input_enr_range"`)*

* `"input_enr_template": string` &mdash; template of input enr filenames.
* `"input_enr_range": [integer, integer]` &mdash; the number range defining input enr filenames.

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns 1, 2 and 3 (without 4!).)*

<a name="traj.example"></a>
### Example input file

```json
{
  "labels": ["FIC", "DNA"],
  "sizes": [1, 8],
  "box_size": 750.0,
  "temperature": 293.15,
  "viscosity": 0.01005,
  "input_xyz_template": "ficoll_41_DNA_14_",
  "input_enr_template": "ficoll_41_DNA_14_",
  "input_xyz_range": [1, 94],
  "input_enr_range": [1, 94],
  "fit_MSD": true,
  "verbose": true,
  "debug": true,
  "probing_frequency": 10,
  "mode": "window"
}
```

<a name="traj.usage"></a>
### Usage
If you have already prepared an input JSON file (using keywords introduced above), you can run the `MSD.py` program using following command:

`python MSD.py input.json`

<a name="traj.output"></a>
### Output files

Successful computations should produce:
* `(TEMPLATE)msd.txt` data file with MSD as a function of time,
* `(TEMPLATE)msd.jpg` image file with MSD as a function of time (and optionally, linear fit),
* `(TEMPLATE)enr.jpg` image file with total energy as a function of time, if `"input_enr_template"` and `"input_enr_range"` are present in the input JSON file.

<a name="hcells"></a>
## H Cell simulator

`H_cell_sim.py` enables one to compute the concentration profiles emerging in the H-cell experiments.

How to use:

1. `make`
2. `make test`
3. `python H_cell_sim.py --help`
4. `python H_cell_sim.py -a 1 -b 1 -d "1 1.3" -n 100 -x 0.000001 -m 1 -o hcellsim -s "0.001 0.35"`
