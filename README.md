# pyBrown

[![Build Status](https://travis-ci.com/tskora/pyBrown.svg?branch=master)](https://travis-ci.com/tskora/pyBrown)

`pyBrown` is a bundle of tools useful for **Brownian** and **Stokesian** dynamics
simulations and subsequent analysis of the results.

Copyright &copy;2018- Tomasz Skóra [tskora@ichf.edu.pl](mailto:tskora@ichf.edu.pl)

## Features

- [x] computing Mean Square Displacement
- [x] computing Mean Square Angular Displacement and Mean Orientation Autocorrelation
- [x] computing Radial Distribution Function
- [x] structure (`.str`) files generation
- [x] Monte Carlo Excluded Volume computation
- [ ] diffusion/mobility/resistance matrix analysis

## First steps

Type following commands in a terminal:

`make`

`make test`

## Table of Contents

0. [General remarks](#gen)
    * [Units](#gen.units)
1. [Generating structure files](#strs)
    * [Keywords](#strs.keywords)
2. [Trajectory analysis](#traj)
    * [`MSD.py`](#traj.msd)
        * [Keywords](#traj.msd.keywords)
        * [Example input file](#traj.msd.example)
        * [Usage](#traj.msd.usage)
        * [Output files](#traj.msd.output)
    * [`Energy.py`](#traj.enr)
        * [Keywords](#traj.enr.keywords)
        * [Example input file](#traj.enr.example)
        * [Usage](#traj.enr.usage)
        * [Output files](#traj.enr.output)
3. [Monte Carlo Excluded Volume](#mcev)
    * [`ExVol.py`](#mcev.ev)
        * [Keywords](#mcev.ev.keywords)
        * [Example input files](#mcev.ev.example)
        * [Usage](#mcev.ev.usage)
        * [Output](#mcev.ev.output)
4. [Snapshot voxelization](#vox)
    * [`Voxels.py`](#vox.vox)
        * [Keywords](#vox.vox.keywords)
        * [Example input file](#vox.vox.example)
        * [Usage](#vox.vox.usage)
        * [Output files](#vox.vox.output)
    * [`PlotVoxels.py`](#vox.plt)

<a name="gen"></a>
## General remarks
<a name="gen.units"></a>
### Units
| Physical property | Units |
|---|---|
| Temperature | kelvin (*K*) |
| Viscosity | poise (*P*) |
| Time | picosecond (*ps*) |
| Distance | Angstrom (*Å*) |

<a name="strs"></a>
## Generating structure files
<a name="strs.keywords"></a>
### Keywords

<a name="traj"></a>
## Trajectory analysis
<a name="traj.msd"></a>
### `MSD.py`
<a name="traj.msd.keywords"></a>
#### Keywords
**Required keywords:**

* `"labels": [string, ...]` &mdash; bead labels in input XYZ file,
* `"sizes": [integer, ...]` &mdash; numbers of beads representing individual entities,
* `"box_size": float` &mdash; size of simulation (cubic) box (*Å*),
* `"temperature": float` &mdash; temperature (*K*),
* `"viscosity": float` &mdash; dynamic viscosity (*P*),
* `"input_xyz_template": string` &mdash; template of input xyz filenames,
* `"input_xyz_range": [integer, integer]` &mdash; the number range defining input xyz filenames.

*pyBrown expects input `xyz` files to follow a specific naming scheme:
..., `(TEMPLATE)(NUMBER).xyz`, `(TEMPLATE)(NUMBER).xyz`, ...
(where TEMPLATE is a string variable defined with the keyword `"input_xyz_template"` and NUMBER is an integer from range defined with the keyword `"input_xyz_range"`)*

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns [1, 2, 3]  (without 4!).)*

**Optional keywords:**

* `"debug": boolean` &mdash; print extra information useful for debugging purposes (default: `false`)
* `"verbose": boolean` &mdash; print extra information (default: `false`)
* `"fit_MSD": boolean` &mdash; fit linear functions to the computed MSDs and plot them in the output figure (default: `false`)
* `"probing_frequency": integer` &mdash; read every *N*-th geometry (default: `1`)
* `"min_time:" float` &mdash; not include snapshots with time smaller than `"min_time"` (default: `0.0`)
* `"mode": option` &mdash; learn more [here](https://freud.readthedocs.io/en/v2.0.1/modules/msd.html) (options: `"direct"`/`"window" `, default: `window`)
* `"float_type": option` &mdash; number of bits per float number (options: `32`/`64`, default: `32`)

<a name="traj.msd.example"></a>
#### Example input file

```json
{
  "labels": ["FIC", "DNA"],
  "sizes": [1, 8],
  "box_size": 750.0,
  "temperature": 293.15,
  "viscosity": 0.01005,
  "input_xyz_template": "ficoll_41_DNA_14_",
  "input_xyz_range": [1, 94],
  "fit_MSD": true,
  "verbose": true,
  "debug": true,
  "probing_frequency": 10,
  "mode": "window"
}
```

<a name="traj.msd.usage"></a>
#### Usage
If you have already prepared an input JSON file (using keywords introduced above), you can run the `MSD.py` program using following command:

`python MSD.py input.json`

<a name="traj.msd.output"></a>
#### Output files

Successful computations should produce:
* `(TEMPLATE)msd.txt` data file with MSD as a function of time,
* `(TEMPLATE)msd.pdf` image file with MSD as a function of time (and optionally, linear fit),
<!-- * `(TEMPLATE)enr.jpg` image file with total energy as a function of time, if `"input_enr_template"` and `"input_enr_range"` are present in the input JSON file. -->

<a name="traj.enr"></a>
### `Energy.py`
<a name="traj.enr.keywords"></a>
#### Keywords
**Required keywords:**

* `"input_enr_template": string` &mdash; template of input enr filenames,
* `"input_enr_range": [integer, integer]` &mdash; the number range defining input enr filenames.

*pyBrown expects input `enr` files to follow a specific naming scheme:
..., `(TEMPLATE)(NUMBER).enr`, `(TEMPLATE)(NUMBER).enr`, ...
(where TEMPLATE is a string variable defined with the keyword `"input_enr_template"` and NUMBER is an integer from range defined with the keyword `"input_enr_range"`)*

*(Have in mind, that ranges in python are defined in such a way that the upper limit is not contained in a range. For example, range(1,4) returns [1, 2, 3]  (without 4!).)*

**Optional keywords:**

* `"debug": boolean` &mdash; print extra information useful for debugging purposes (default: `false`)
* `"verbose": boolean` &mdash; print extra information (default: `false`)
* `"float_type": option` &mdash; number of bits per float number (options: `32`/`64`, default: `32`)

<a name="traj.enr.example"></a>
#### Example input file

```json
{
  "input_enr_template": "enzymes_double_1_",
  "input_enr_range": [1, 11],
  "verbose": false,
  "debug": false
}
```

<a name="traj.enr.usage"></a>
#### Usage
If you have already prepared an input JSON file (using keywords introduced above), you can run the `Energy.py` program using following command:

`python Energy.py input.json`

<a name="traj.enr.output"></a>
#### Output files

Successful computations should produce:
* `(TEMPLATE)enr.txt` data file with mean energy as a function of time,
* `(TEMPLATE)enr.pdf` image file with mean energy as a function of time (and optionally, linear fit)

<a name="mcev"></a>
## Monte Carlo Excluded Volume
<a name="mcev.ev"></a>
### `ExVol.py`
<a name="mcev.ev.keywords"></a>
#### Keywords
**Required keywords:**

* `"box_size": float` &mdash; size of simulation (cubic) box (*Å*),
* `"tracer_radii": [float, ...]` &mdash; radii of tracers randomly inserted into the box: if list contains one value, spherical particle is being inserted; else, linear particle composed of beads of given radii is being inserted (*Å*),
* `"number_of_trials": integer` &mdash; number of Monte Carlo insertions per snapshot,
* `"input_str_filename": string` &mdash; input str file, from which crowder sizes are loaded,
* `"radii_mode": option` &mdash; bead radius definition (options: `"hydrodynamic"`/`"lennard-jones"`)
* `"times": [float, ...]` &mdash; times for which configurations are loaded from `.xyz` files,
* `"input_xyz_template": string` &mdash; template of input xyz filenames,
* `"input_xyz_range": [integer, integer]` &mdash; the number range defining input xyz filenames.

*pyBrown expects input `xyz` files to follow a specific naming scheme:
..., `(TEMPLATE)(NUMBER).xyz`, `(TEMPLATE)(NUMBER).xyz`, ...
(where TEMPLATE is a string variable defined with the keyword `"input_xyz_template"` and NUMBER is an integer from range defined with the keyword `"input_xyz_range"`)*

**Optional keywords:**

* `"debug": boolean` &mdash; print extra information useful for debugging purposes (default: `false`)
* `"verbose": boolean` &mdash; print extra information (default: `false`)
* `"bond_lengths": option` &mdash; way of defining separations between beads in case of nonspherical tracer (default (*for now only option*): `"hydrodynamic_radii"`)
* `"withdraw": [string, ...]` &mdash; list of labels representing particles that are randomly withdrawn before every insertion (default: `[]`),
* `"scan_mode": boolean` &mdash; if `false`: excluded volume is computed for tracer defined by `"tracer_radii"`; if `true`: excluded volume is computed for tracers of size from `0.0` to `"tracer_radii"` (default: `false`)
* `"scan_density": integer` &mdash; number of tracer radii scanned in `"scan_mode"` (default: `2`)
* `"float_type": option` &mdash; number of bits per float number (options: `32`/`64`, default: `32`)

<a name="mcev.ev.example"></a>
#### Example input files

```json
{
  "box_size": 750.0,
  "tracer_radii": [0.0],
  "number_of_trials": 1000,
  "input_xyz_template": "ficoll_42_",
  "input_xyz_range": [1, 4],
  "input_str_filename": "ficoll_42_1.str",
  "radii_mode": "hydrodynamic",
  "times": [0.0, 1000000.0, 2000000.0, 3000000.0, 4000000.0, 4500000.0],
  "withdraw": ["FIC", "FIC"]
}
```
```json
{
  "box_size": 750.0,
  "tracer_radii": [51],
  "number_of_trials": 1000,
  "input_xyz_template": "ficoll_42_",
  "input_xyz_range": [1, 4],
  "input_str_filename": "ficoll_42_1.str",
  "radii_mode": "hydrodynamic",
  "times": [0.0, 1000000.0, 2000000.0, 3000000.0, 4000000.0, 4500000.0],
  "withdraw": [],
  "scan_mode": true,
  "scan_density": 3
}
```

<a name="mcev.ev.usage"></a>
#### Usage
If you have already prepared an input JSON file (using keywords introduced above), you can run the `ExVol.py` program using following command:

`python ExVol.py input.json`

<a name="mcev.ev.output"></a>
#### Output files

Successful computations should produce:
* `*fex.txt` data file with excluded volume fractions as a function of tracer radii,

<a name="vox"></a>
## Snapshot voxelization
<a name="vox.vox"></a>
### `Voxels.py`
<a name="vox.vox.keywords"></a>
#### Keywords
**Required keywords:**

* `"box_size": float` &mdash; size of simulation (cubic) box (*Å*),
* `"input_xyz_filename": string` &mdash; input xyz filename,
* `"input_str_filename": string` &mdash; input str filename,
* `"snapshot_time": float` &mdash; time defining single snapshot from `.xyz` file (*ps*),
* `"grid_density": integer` &mdash; number of voxels per box side,
* `"radii_mode": option` &mdash; bead radius definition (options: `"hydrodynamic"`/`"lennard-jones"`)

**Optional keywords:**

* `"debug": boolean` &mdash; print extra information useful for debugging purposes (default: `false`)
* `"verbose": boolean` &mdash; print extra information (default: `false`)
* `"float_type": option` &mdash; number of bits per float number (options: `32`/`64`, default: `32`)

<a name="vox.vox.example"></a>
#### Example input file

```json
{
  "box_size": 750.0,
  "input_xyz_filename": "cytoplasm_250_1.xyz",
  "input_str_filename": "cytoplasm_250_1.str",
  "verbose": true,
  "debug": true,
  "grid_density": 10,
  "snapshot_time": 1629699.184723,
  "radii_mode": "lennard-jones"
}
```

<a name="vox.vox.usage"></a>
#### Usage
If you have already prepared an input JSON file (using keywords introduced above), you can run the `Voxels.py` program using following command:

`python Voxels.py input.json`

<a name="vox.vox.output"></a>
#### Output files

Successful computations should produce:
* `*voxels.txt` data file with digitized snapshot,

<a name="vox.plt"></a>
### `PlotVoxels.py`

If you have already voxelized snapshot you should have `*voxels.txt` file. You can run the `PlotVoxels.py` program using following command:

`python PlotVoxels.py *voxels.txt`