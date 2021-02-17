
# IMPORTS
from os import sep, sys
import numpy as np 
import sys 

import mdutilities as mdu
import observables as obs 
from particle3D import Particle3D as p3d 
import pbc 

def force_LJ(pair_separation:np.array, l:float) -> np.array:
  '''
  Calculates the MIC-force (on particle i due to particle j)
  given their real separation vector r_ij
  '''
  mic_pair_separation = pbc.minimum_image(pair_separation,l)
  r = np.linalg.norm(mic_pair_separation)
  F = -48*(r**(-14)-0.5*r**(-8))*mic_pair_separation # -ve due to using r_ij in formula
  return F


def force_matrix(separation_matrix:np.ndarray, l:float) -> np.ndarray:
  '''
  Calculates the [N,N,3] matrix of forces F_ij 
  (force on p_i due to p_j) for i < j
  '''
  N = len(separation_matrix)
  F_matrix = np.zeros([N,N,3])
  for i in range(N-1):
    for j in range(i+1,N):
      F_matrix[i][j] = force_LJ(separation_matrix[i][j],l)
  return F_matrix


def main():

  #--- SETTING VARIABLES AND PARAMETERS ---------------------------------- 

  print('---\nReading from input files...')
  # Read needed information from command line
  '''
  if len(sys.argv) not in [3,4]:
    print("Wrong number of arguments")
    print(f"Usage: {sys.argv[0]} <setup file> <data file> <OPTIONAL: output file>")
    quit()
  setup_file = f'output/{sys.argv[1]}'
  data_file = f'output/{sys.argv[2]}'
  traj_file = f'output/{sys.argv[3]}' if len(sys.argv) == 4 else 'traj.xyz'
  '''
  setup_file = 'setup.dat'
  data_file = 'data.dat'
  traj_file = 'output\\traj.xyz'

  # Get information from input files
  with open(setup_file, 'r') as setup:
    lines = setup.readlines()
    try:
      temp, mass, density = map(float, lines[1::2]) 
    except: 
      print(f'error in {setup_file}\ncheck all variables are present and their types are correct')
      quit()

  with open(data_file, 'r') as data:
    lines = data.readlines()
    try:
      N = int(lines[1]) 
      dt, t_total = float(lines[3]),float(lines[5])

      # map '0' -> 0 -> False and '1' -> 1 -> True
      log_kinetic_e, log_potential_e, log_total_e, log_msd, log_rdf = map(
          bool, map(int, lines[7::2])) 
    except:
      print(f'error in {data_file}\ncheck all variables are present and their types are correct')
      quit()

    if log_total_e == True: # Override energies
      log_kinetic_e = True
      log_potential_e = True
  
  print(' all variables have been set correctly')

  '''
  # debugging
  print(temp, mass, density)
  print(N, dt,t_total)
  print(log_kinetic_e, log_potential_e, log_total_e, log_msd, log_rdf)
  '''

  # --- INITIALISATION -------------------------------------------------

  print(f'Initialising {N} particles:')
  # Create the N particles with dummy pos and vel
  p3d_list = p3d.new_particles(N, mass)
  
  # Set initial conditions for particles
  cell_length, full_lattice = mdu.set_initial_positions(density, p3d_list)
  mdu.set_initial_velocities(temp, p3d_list)

  # --- SIMULATION -------------------------------------------------------

  # Open output files:

  outfile_traj = open(traj_file,'w')

  if log_total_e:
    outfile_total_e = open('output\\energy_total.csv','w')
  if log_kinetic_e:
    outfile_kinetic_e = open('output\\energy_kinetic.txt', 'w')
  if log_potential_e:
    outfile_potential_e = open('output\\energy_potential.txt', 'w')
  if log_msd:
    initial_positions = np.array([particle.pos for particle in p3d_list])
    outfile_msd = open('output\\msd.csv','w')
  if log_rdf:
    outfile_rdf = open('output\\rdf.csv','w')

  # Initialize objects for use in simulation

  t = 0.0
  total_steps = int(t_total / dt)
  res = 1000
  if log_rdf:
    # Histogram of 500 columns (distances from 0 to l*sqrt(3))
    histogram = np.zeros(res)

  # Print what gets logged and calculated
  print('Logged/Calculated Observables:')
  print('-Particle Trajectories')
  if log_kinetic_e: print('-Kinetic Energy')
  if log_potential_e: print('-Potential Energy')
  if log_total_e: print('-Total Energy')
  if log_msd: print('-Mean Squared Displacement')
  if log_rdf: print('-Radial Distribution Function')

  # Simulation Loop
  print('\nSimulating...')

  for time_step in range(total_steps):

    # Calculate the Separation array
    separation_matrix = p3d.pair_separations(p3d_list)


    # --- Logging Observables ---------------------------------------------

    # Trajectories:
    outfile_traj.write(f'{N}\n')
    outfile_traj.write(f'time = {t}\n')
    for particle in p3d_list:
      outfile_traj.write(f'{particle}\n')

    # Energies:
    if log_kinetic_e:
      outfile_kinetic_e.write(f'{t:.4} {p3d.sys_kinetic(p3d_list)}\n')

    if log_potential_e:
      outfile_potential_e.write(f'{t:.4} {obs.potential_energy(separation_matrix)}\n')

    # Mean Squared Displacement:
    if log_msd and (time_step % 10) == 0:
      particle_positions = np.array([particle.pos for particle in p3d_list])
      outfile_msd.write(f'{t:.4f}, {obs.msd(initial_positions, particle_positions, cell_length):.5f}\n')

    # Radial Distribution Function:
    if log_rdf and (time_step % 10) == 0:
      histogram += obs.rdf(separation_matrix, cell_length, res)

    # --- Particle Time Integrator ----------------------------------------

    # calculate all the net forces for each particle
    F_matrix = force_matrix(separation_matrix,cell_length)
    
    '''
    print(f't = {t}\n') #%
    print(f'separation:\n{separation_matrix}')
    print(f'forces:\n{F_matrix}\n')
    '''

    F_list = []
    for i in range(N):
      F_i = np.zeros(3)
      for j in range(N):
        if i < j: 
          F_i += F_matrix[i][j]
        elif i > j: 
          F_i -= F_matrix[j][i]
      F_list.append(F_i)
    
    # print(f'f_list:\n{F_list}') #%

    # update particle positions
    for i, particle in enumerate(p3d_list):
      particle.update_pos_2nd(cell_length,dt,F_list[i])

    # update forces
    separation_matrix = p3d.pair_separations(p3d_list)
    F_matrix = force_matrix(separation_matrix, cell_length)
    F_new_list = []
    for i in range(N):
      F_i = np.zeros(3)
      for j in range(N):
        if i < j:
          F_i += F_matrix[i][j]
        elif i > j:
          F_i -= F_matrix[j][i]
      F_new_list.append(F_i)

    # update particle velocities
    for i, particle in enumerate(p3d_list):
      particle.update_vel(dt,0.5*(F_list[i]+F_new_list[i]))

    # update time
    t += dt
    

  # Close all simulation files

  outfile_traj.close()
  if log_kinetic_e: outfile_kinetic_e.close()
  if log_potential_e: outfile_potential_e.close()
  if log_msd: outfile_msd.close()

  # --- POST SIMULATION CALCULATIONS --------------------------------

  # Total Energy
  if log_total_e:
    outfile_kinetic_e = open('output\\energy_kinetic.txt','r')
    outfile_potential_e = open('output\\energy_potential.txt','r')
    kinetic_lines = outfile_kinetic_e.readlines()
    potential_lines = outfile_potential_e.readlines()
    for i in range(total_steps):
      t, kinetic_e = kinetic_lines[i].split(' ')
      potential_e = potential_lines[i].split(' ')[1]
      outfile_total_e.write(f'{t}, {float(kinetic_e)+float(potential_e)}\n')

    outfile_total_e.close()

  # Radial Distribution Function
  if log_rdf:
    histogram = obs.rdf_normalize(histogram, cell_length, res)
    max_length = cell_length*(np.sqrt(3)/2+0.01)
    for i in range(res):
      outfile_rdf.write(f'{i*max_length/res:.4f}, {histogram[i]}\n')

    outfile_rdf.close()

  print(' done.')

  


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
