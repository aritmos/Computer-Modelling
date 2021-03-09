# LJ SIMULATION

A program to simulate particles under the Lennard-Jones potential utilizing
periodic boundary conditions.

author: -redacted-
-------------------------------------------------------------------------------
## Requirements:
- numpy 

## Use:
The main program <lj_sim.py> must be run from terminal along with two auxiliary
files that set the parameters of the simulation, and optionally a third parameter
to specify the name of the output file (defaults to traj.xyz if not specified)

> py lj_sim.py setup.dat data.dat

Although the name of the input files does not matter, their formatting does. 
Do not change the order of any variables as compared to the provided input files.
Simply modify the uncommented lines to vary the simulation parameters. 
The calculated observables are selected within the equivalent of 'data.dat',
 in order to calculate/log the observable give it a value of 1, if not, then 0.
All other numerical simulation parameters are in their expected type, 
either int or float respectively. 

## Reduced Units:
All inputs, outputs and variables within the code are in the following 
reduced units:
- length: r* = r / σ
- energy: E* = E / ε
- mass: m* = m / m0  (for particles of mass m0)
- time: t* = t / τ where τ = σ*sqrt{m/ε}

Where σ and ε are parameters in the LJ potential.
(these depend on what is being simulated)

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

25/2
- changed how the coefficient is calculated for RDF (now done in rdf_norm)
- changed how forces are calculated by using [N,N,3] matrices
