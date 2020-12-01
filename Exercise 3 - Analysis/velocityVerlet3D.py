"""
CMod Ex3: Velocity Verlet time integration of 
two particles moving in a Morse potential.

Produces plots of the particle separations
and the total energy, both as function of time. Also
saves both to file.

The potential is U = De*((1-np.exp(-a*(r12-re)))**2-1)
where re, De and a are read in from an input file
along with the initial properties of the particles
and passed to the functions that
calculate force and potential energy.
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D


def force_morse(p1, p2, re, De, a) -> np.array:
    """
    Method to return the force on particle 1
    in a morse potential (caused by two particles).
    Force is given by
    F1(r1,r2) = 2*a*De[1-exp(-a(|r1-r2| -re))]exp[-a(|r1-r2| -re)] in the unit direction of (r2-r1)

    :param p1: Particle3D instance 
    :param p2: Particle3D instance
    :param re: parameter re from potential
    :param De: parameter De from potential
    :param a: parameter a from potential
    :return: force acting on particle as Numpy array
    """
    vec_r12: np.array = p2.pos - p1.pos
    r12 = np.linalg.norm(vec_r12)
    F = 2*a*De*(1-np.exp(-a*(r12-re)))*np.exp(-a*(r12-re))*vec_r12/r12
    return F


def pot_energy_morse(p1: Particle3D, p2: Particle3D, re: float, De: float, a: float) -> float:
    """
    Method to return potential energy 
    of particles in morse potential
    U_M(r1,r2)=De[(1-exp(-a(|r1-r2| -re)))^2-1]

    :param p1: Particle3D instance 
    :param p2: Particle3D instance
    :param re: parameter re from potential
    :param De: parameter De from potential
    :param a: parameter a from potential
    :return: potential energy of particle as float
    """
    r12 = np.linalg.norm(p2.pos-p1.pos)
    U = De*((1-np.exp(-a*(r12-re)))**2-1)
    return U


def separation(p1: Particle3D, p2: Particle3D) -> float:
    """
    Calculates the distance between two Particle3D instances
    :param p1: first particle
    :param p2: second particle
    :return: separation distance
    """
    return np.linalg.norm(p1.pos-p2.pos)


# Begin main code

def main(inputfile_name:str, dt:float):

    # Set up simulation time parameters
    total_t = 20
    numstep = int(total_t/dt)
    t = 0.0

    # Set the decimal point precission for the output file
    dp = len(str(dt)[str(dt).find('.'):])-1

    # Get Morse Potential Parameters from input file
    inputfile = open(inputfile_name, 'r')
    lines = inputfile.readlines()
    # Parameters are stored in line 4 (index 3)
    re, De, a = [float(i) for i in lines[3].split(' ')]
    inputfile.close()

    # Create particles from input file
    inputfile = open(inputfile_name, 'r')
    # Particles are in lines 7-8 (newparticle method takes care adjusting indices)
    p1, p2 = Particle3D.new_particle(inputfile, 7, 8)
    inputfile.close()

    # Open output files
    outfile_sep = open(f'separation.txt', 'w')
    outfile_energy = open(f'energy.txt', 'w')

    # Write out initial conditions
    energy = Particle3D.sys_kinetic(
        [p1, p2])+pot_energy_morse(p1, p2, re, De, a)

    outfile_sep.write(f'{t:.{dp}f} {separation(p1,p2):.15f}\n')
    outfile_energy.write(f'{t:.{dp}f} {energy:.15f}\n')

    # Initialise data lists for plotting later
    t_list = [t]
    sep_list = [separation(p1, p2)]
    energy_list = [energy]

    # Start the t integration loop
    for _ in range(numstep):
        # Update particle positions
        force = force_morse(p1,p2,re,De,a)
        p1.update_pos_2nd(dt, force)
        p2.update_pos_2nd(dt, -force)

        # Update force
        force_new = force_morse(p1,p2, re, De, a)
        # Update particle velocity by averaging
        # current and new forces
        p1.update_vel(dt, 0.5*(force+force_new))
        p2.update_vel(dt, -0.5*(force+force_new))

        # Re-define force value
        force = force_new

        # Increase t
        t += dt

        # Output particle information
        energy = Particle3D.sys_kinetic([p1,p2]) + pot_energy_morse(p1,p2,re,De,a)
        outfile_sep.write(f'{t:.{dp}f} {separation(p1,p2):.15f}\n')
        outfile_energy.write(f'{t:.{dp}f} {energy:.15f}\n')

        # Append information to data lists
        t_list.append(t)
        sep_list.append(separation(p1,p2))
        energy_list.append(energy)

    # Post-simulation:
    # Close output file
    outfile_sep.close()
    outfile_energy.close()

    '''
    # Plot particle trajectory to screen
    pyplot.title(f'Velocity Verlet: separation vs t (dt = {dt})')
    pyplot.xlabel('t')
    pyplot.ylabel('Separation')
    pyplot.plot(t_list, sep_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title(f'Velocity Verlet: total energy vs time (dt = {dt})')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(t_list, energy_list)
    pyplot.show()
    '''
    
