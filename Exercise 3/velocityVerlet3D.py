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


def main():
    # Read name of output file from command line
    if len(sys.argv) != 3:
        print("Wrong number of arguments.")
        print(f"Usage: {sys.argv[0]} <input file> <output file prefix>")
        quit()
    else:
        inputfile_name = sys.argv[1]
        outfile_name = sys.argv[2]

    # Open input file
    inputfile = open(inputfile_name, "r")

    # Set up simulation parameters
    total_time = 20
    dt = 0.01
    numstep = int(total_time/dt)
    time = 0.0

    # Set the decimal point precission for the output file
    dp = len(str(dt)[str(dt).find('.'):])-1

    # Set up simulation parameters for potential and force and create particles
    lines = inputfile.readlines()
    first_line = True
    particle_list = []
    for line in lines:
        if not line.startswith('#'):  # If the line is not a comment
            # First line has the potential and and force parameters (re,De,a)
            if first_line:
                split = line.split(' ')
                if len(split) != 3:
                    raise IndexError(
                        'First uncommented line in input file needs to have three space sparated values: <re> <De> <a>')
                global re, De, a
                try:
                    re, De, a = [float(i) for i in split]
                except:
                    TypeError(
                        'The 3 parameters in the first uncommented line of the input file must be numbers')
                first_line = False  # if the variables were initialised then stop looking for the first line
            # if its not the (re De a) line create a new particle
            else: particle_list.append(Particle3D.new_particle(line))

    p1 = particle_list[0]
    p2 = particle_list[1]

    inputfile.close()

    # Open output file
    outfile_sep = open(f'{outfile_name}_VV_sep_{dt}.txt', 'w')
    outfile_energy = open(f'{outfile_name}_VV_energy_{dt}.txt', 'w')

    # Write out initial conditions
    energy = Particle3D.sys_kinetic(
        [p1, p2])+2*pot_energy_morse(p1, p2, re, De, a)

    outfile_sep.write(f'{time:.{dp}f} {separation(p1,p2):.15f}\n')
    outfile_energy.write(f'{time:.{dp}f} {energy:.15f}\n')

    # Initialise data lists for plotting later
    time_list = [time]
    sep_list = [separation(p1, p2)]
    energy_list = [energy]

    # Start the time integration loop
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

        # Increase time
        time += dt

        # Output particle information
        energy = Particle3D.sys_kinetic([p1,p2]) + 2*pot_energy_morse(p1,p2,re,De,a)
        outfile_sep.write(f'{time:.{dp}f} {separation(p1,p2):.15f}\n')
        outfile_energy.write(f'{time:.{dp}f} {energy:.15f}\n')

        # Append information to data lists
        time_list.append(time)
        sep_list.append(separation(p1,p2))
        energy_list.append(energy)

    # Post-simulation:
    # Close output file
    outfile_sep.close()
    outfile_energy.close()

    # Plot particle trajectory to screen
    pyplot.title(f'Symplectic Euler: separation vs time, dt = {dt}')
    pyplot.xlabel('Time')
    pyplot.ylabel('Separation')
    pyplot.plot(time_list, sep_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title(f'Symplectic Euler: total energy vs time, dt = {dt}')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
