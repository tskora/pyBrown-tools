# pyBrown

`pyBrown` is a bundle of tools useful for **Brownian** and **Stokesian** dynamics
simulations and squbsequent analysis of the results.

Copyright ©2018-2019 Tomasz Skóra [tskora@ichf.edu.pl](mailto:tskora@ichf.edu.pl)

## Table of contents
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

<a name="hcells"></a>
## H Cell simulator

`H_cell_sim.py` enables one to compute the concentration profiles emerging in the H-cell experiments.

How to use:

1. `make`
2. `make test`
3. `python H_cell_sim.py --help`
4. `python H_cell_sim.py -a 1 -b 1 -d "1 1.3" -n 100 -x 0.000001 -m 1 -o hcellsim -s "0.001 0.35"`
