# LJ SIMULATION

A program to simulate particles under the Lennard-Jones potential utilizing
periodic boundary conditions.

author: -redacted-
-------------------------------------------------------------------------------
## Requirements:
- numpy library

## Use:
The main program <lj_sim.py> must be run from terminal along with two auxiliary
files that set the parameters of the simulation, and optionally a third parameter
to specify the name of the output file (defaults to traj.xyz if not specified)

For example, when simulating the gas phase:

> py lj_sim.py setup_gas.dat data_gas.dat traj_gas.xyz

Note that the order of the files matters.

## Reduced Units:
All inputs, outputs and variables within the code are in the following 
reduced units:
- length: r* = r / σ
- energy: E* = E / ε
- mass: m* = m / m0  (for particles of mass m0)
- time: t* = t / τ where τ = σ*sqrt{m/ε}

Where σ and ε are parameters in the LJ potential.
(these depend on what is being simulated)

## Input files:

lj_sim.py requires two existing files, which contain the parameters for the simulation.
Using a generalized form for naming convention in the previous section, we will call these 
<setup.dat> and <data.dat> respectively.
One can duplicate and modify the existing files <setup.dat> and <data.dat>,
as long as the formatting is not changed, as the variables are set by looking at hardcoded
line numbers within their respective files.

Every odd line has a comment that specifies what parameter the number on the line below refers to.
The value of each parameter should be the only thing in its respective line.
In the case of disabling and enabling observables use 0 and 1 respectively.

### setup.dat
This file contains the variables to set the properties of the material that is simulated.
It contains the temperature, mass and density. 
Since the program works in reduced units the mass should be always set to 1.

### data.dat
This file contains the variables to set the properties of the simulation.
It contains the number of simulated particles, the total simulation time, the simulation 
time step, and which observables to log.

Observables: potential energy, kinetic energy, total energy
(total energy will override potential and kinetic if enabled), MSD and RDF.

Furthermore there are also variables to choose every how many time steps we want the MSD
and RDF to be calculated, and how many columns the RDF histogram should have.

## Outputs

All outputs will be inside the <output> subdirectory. The simulation will always at least
create the trajectory file, which by default is <traj.xyz>, any of the other enabled observables
will also be saved here. The format for the trajectory file is XYZ compliant in order for it to be
visualized in a program such as VMD.

The files for the observables (except RDF) have the following structure, at each line:

> time value

The RDF instead has:

> distance value

## Changelog: (from original project document)

31/1
- changed mdutilities to return box_length instead of box_size
- changed variables in data.dat to use 0/1 instead of False/True 
- new_particle now only takes in N and m and returns dummy list
- changed position and velocity in mdutilities to pos, vel.
- changed 'new_particle' to 'new_particles'
- spacing on print statements for cleaner output
- added mass to setup.dat

6/2
- changed p3d instance label from i to p{i} for better readability in traj.xyz
- changed p3d instance label to start at 0 for better working with arrays

10/2
- changed force_i to new_force_i 

12/2
- made pbc functions input and output np.arrays instead of list
- removed unecessary parameters in setup.dat

17/2
- made output files save to a directory 'output'
- obs.msd now takes in initial_positions and p3d_list

18/2
- added net_forces() to calculate the net forces list 

19/2
- added rdf and msd parameters in data.dat
- added plotting module to plot observables

25/2
- changed how the coefficient is calculated for RDF (now done in rdf_norm)
- changed how forces are calculated by using [N,N,3] matrices
