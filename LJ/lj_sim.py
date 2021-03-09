"""
CompMod Final Project: LJ_Sim, a program to simulate particles interacting
via a Lennard-Jones potential with periodic boundary conditions in the shape
of a cube. Can also calculate:
- Kinetic Energy
- Potential Energy
- Total Energy
- MSD
- RDF

This main program must be run from terminal along with two auxiliary files
that set the parameters of the simulation, and optionally a third parameter
to specify the name of the output file (defaults to traj.xyz)

> py lj_sim.py setup.dat data.dat

All inputs, outputs and variables within the code are in the following 
reduced units:
- length: r* = r / σ
- energy: E* = E / ε
- mass: m* = m / m0  (for particles of mass m0)
- time: t* = t / τ where τ = σ*sqrt{m/ε}

Where σ and ε are parameters in the LJ potential,
which depend on the simulated matter.

author: -redacted-

"""

# IMPORTS
from os import sys

import numpy as np

import mdutilities as mdu
import observables as obs
from particle3D import Particle3D as p3d
import pbc
import plotting


def force_LJ(mic_pair_separation: np.array) -> np.array:
    """
    Calculates the MIC-force on particle i due to particle j)
    given their mic separation vector r_ij

    :param mic_pair_separation: 3 vector or mic separation between two particles
    :return: LJ force
    """
    r = np.linalg.norm(mic_pair_separation)
    if r > 3.5:  # Implement cutoff radius
        return np.zeros(3)
    else:
        F = (
            -48 * (r ** (-14) - r ** (-8) / 2) * mic_pair_separation
        )  # negative due to using r_ij here
        return F


def force_matrix(mic_separations: np.ndarray) -> np.ndarray:
    """
    Calculates the [N,N,3] matrix of forces F_ij
    (force on p_i due to p_j) for i < j

    :param mic_separations: [N,N,3] matrix of mic particle separations
    :return: upper triangular [N,N,3] force matrix
    """
    N = len(mic_separations)
    F_matrix = np.zeros([N, N, 3])
    for i in range(N - 1):
        for j in range(i + 1, N):
            F_matrix[i][j] = force_LJ(mic_separations[i][j])
    return F_matrix


def net_forces(F_matrix: np.ndarray) -> np.ndarray:
    """
    Calulates the net forces on every particle given an [N,N,3] matrix of forces

    :param F_matrix: [N,N,3] matrix of forces
    :return: [N,3] matrix of net forces
    """
    N = len(F_matrix)
    F_list = []
    for i in range(N):
        F_i = np.zeros(3)
        for j in range(N):
            if i < j:
                F_i += F_matrix[i][j]
            elif i > j:
                F_i -= F_matrix[j][i]
        F_list.append(F_i)
    return np.array(F_list)


def mic_matrix(matrix: np.ndarray, l: float) -> np.ndarray:
    """
    Turns the [i][j] elements of an [N,N,3] matrix into their mic image

    :param matrix: the [N,N,3] input matrix
    :param l: the cube box length to apply mic
    :return: the mic matrix
    """
    N = len(matrix)
    mic_matrix = np.zeros([N, N, 3])
    for i in range(N):
        for j in range(N):
            mic_matrix[i][j] = pbc.minimum_image(matrix[i][j], l)
    return mic_matrix


def main():

    # SETTING VARIABLES AND PARAMETERS ========================================

    print("---\nReading from input files")
    # Read needed information from command line
    if len(sys.argv) not in [3, 4]:
        print("Wrong number of arguments")
        print(
            f"Usage: {sys.argv[0]} <setup file> <data file> <OPTIONAL: output file>")
        quit()
    setup_file = sys.argv[1]
    data_file = sys.argv[2]
    traj_file = sys.argv[3] if len(sys.argv) == 4 else 'traj.xyz'

    # Get information from input files
    with open(setup_file, "r") as setup:
        lines = setup.readlines()
        try:
            TEMP, MASS, DENSITY = map(float, lines[1::2])
        except:
            print(
                f'Error in {setup_file}\ncheck all variables are'
                f'present and their types are correct'
            )
            quit()

    with open(data_file, "r") as data:
        lines = data.readlines()
        try:
            N = int(lines[1])
            dt = float(lines[3])
            t_total = float(lines[5])

            # map '0' -> 0 -> False and '1' -> 1 -> True
            LOG_KINETIC_E, LOG_POTENTIAL_E, LOG_TOTAL_E, LOG_MSD, LOG_RDF = map(
                bool, map(int, lines[7:16:2])
            )
            t_multiple = int(lines[17])
            res = int(lines[19])
        except:
            print(
                f"Error in {data_file}\ncheck all variables"
                "are present and their types are correct"
            )
            quit()

        if LOG_TOTAL_E == True:  # Override energies
            LOG_KINETIC_E = LOG_POTENTIAL_E = True

    # check that all the required variables are positive
    positive_parameters = [N, dt, t_total, t_multiple, res]
    for param in positive_parameters:
        if param <= 0:
            print(f"Error: tried to initialize a parameter "
                  "with a non-positive value")
            quit()

    print(" all variables have been set correctly")

    # Initialization ==========================================================

    print(f"Initialising {N} particles:")

    # Create the N particles with dummy positisons and velocities
    p3d_list = p3d.new_particles(N, MASS)

    # Set initial conditions for particles
    cell_length = mdu.set_initial_positions(DENSITY, p3d_list)[0]
    mdu.set_initial_velocities(TEMP, p3d_list)

    # Initialize remaining objects for use in simulation

    t = 0.0
    total_steps = int(t_total / dt)

    separation_array = p3d.pair_separations(p3d_list)
    mic_separations = mic_matrix(separation_array, cell_length)

    F_matrix = force_matrix(mic_separations)
    F_new_list = net_forces(F_matrix)

    # SIMULATION ==============================================================

    # Open output files, print and intialize content for enabled observables:

    print("Logged/Calculated Observables:")

    print("-Particle Trajectories")
    outfile_traj = open(f"output\\{traj_file}", "w")

    if LOG_KINETIC_E:
        print("-Kinetic Energy")
        outfile_kinetic_e = open("output\\energy_kinetic.txt", "w")
    if LOG_POTENTIAL_E:
        print("-Potential Energy")
        outfile_potential_e = open("output\\energy_potential.txt", "w")
    if LOG_TOTAL_E:
        print("-Total Energy")
        outfile_total_e = open("output\\energy_total.csv", "w")
    if LOG_MSD:
        print(
            f'-Mean Squared Displacement\n   (Every '
            f'{t_multiple} step{"s" if t_multiple > 1 else ""})'
        )
        initial_positions = np.array([particle.pos for particle in p3d_list])
        outfile_msd = open("output\\msd.csv", "w")
    if LOG_RDF:
        print(
            f'-Radial Distribution Function\n   (Every '
            f'{t_multiple} step{"s" if t_multiple > 1 else ""} with {res} columns)'
        )
        histogram = np.zeros(res)
        max_distance = cell_length * (
            np.sqrt(3) / 2 + 0.001
        )  # added extra distance to avoid floating point rounding overflows
        dr = max_distance / res
        outfile_rdf = open("output\\rdf.csv", "w")

    print(f"\nTotal simulation time: {t_total}")
    print(f"Time step: {dt}")

    # Simulation Loop =========================================================
    print("\nSimulating...")

    for t_step in range(total_steps):

        # --- Logging Observables ---------------------------------------------

        # Trajectories:
        outfile_traj.write(f"{N}\n")
        outfile_traj.write(f"time = {t:.4f}\n")
        for particle in p3d_list:
            outfile_traj.write(f"{particle}\n")

        # Energies:
        sys_kinetic = p3d.sys_kinetic(p3d_list)
        sys_potential = obs.potential_energy(mic_separations)

        if LOG_KINETIC_E:
            outfile_kinetic_e.write(f"{t:.3f} {sys_kinetic:.4f}\n")

        if LOG_POTENTIAL_E:
            outfile_potential_e.write(f"{t:.3f} {sys_potential:.4f}\n")

        if LOG_TOTAL_E:
            outfile_total_e.write(
                f"{t:.3f} {(sys_kinetic+sys_potential):.4f}\n")

        # Mean Squared Displacement:
        if LOG_MSD and (t_step % t_multiple) == 0:
            outfile_msd.write(
                f"{t:.3f}, {obs.msd(initial_positions, p3d_list, cell_length):.5f}\n")

        # Radial Distribution Function:
        if LOG_RDF and (t_step % t_multiple) == 0:
            histogram += obs.rdf(mic_separations, max_distance, res)

        # --- Particle Time Integrator (Velocity Verlet)-----------------------

        # update particle positions
        F_list = F_new_list
        for i, particle in enumerate(p3d_list):
            particle.update_pos_2nd(cell_length, dt, F_list[i])

        # update forces
        separation_array = p3d.pair_separations(p3d_list)
        mic_separations = mic_matrix(separation_array, cell_length)

        F_matrix = force_matrix(mic_separations)
        F_new_list = net_forces(F_matrix)

        # update particle velocities
        for i, particle in enumerate(p3d_list):
            F_avg = 0.5 * (F_list + F_new_list)
            particle.update_vel(dt, F_avg[i])

        # update time
        t += dt
    # =========================================================================

    # Close all simulation files
    outfile_traj.close()
    if LOG_KINETIC_E:
        outfile_kinetic_e.close()
    if LOG_POTENTIAL_E:
        outfile_potential_e.close()
    if LOG_TOTAL_E:
        outfile_total_e.close()
    if LOG_MSD:
        outfile_msd.close()

    # --- POST SIMULATION CALCULATIONS ----------------------------------------

    # Normalize the Radial Distribution Function
    if LOG_RDF:
        t_count = total_steps / t_multiple
        histogram = obs.rdf_normalize(histogram, t_count, dr, N, DENSITY)

        # Write out the histogram to a file
        for i in range(res):
            outfile_rdf.write(f"{(i/res)*max_distance:.4f}, {histogram[i]}\n")
        outfile_rdf.close()

    print(" done.")

    # Plot the logged observables
    plotting.plot(LOG_TOTAL_E, LOG_MSD, LOG_RDF)


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
