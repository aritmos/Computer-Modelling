
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
  Calculates the MIC-force between two particles
  given their real separation vector
  '''
  mic_pair_separation = np.array(pbc.minimum_image(pair_separation,l))
  r = np.linalg.norm(mic_pair_separation)
  F = -48*(r**(-14)-0.5*r**(-8))*mic_pair_separation
  return F


def force_matrix(separation_array:np.ndarray, l:float) -> np.ndarray:
  '''
  Calculates the [N,N,3] matrix of forces F_ij 
  (force on p_i due to p_j) for i < j
  '''
  N = len(separation_array)
  F_matrix = np.zeros([N,N,3])
  for i in range(N-1):
    for j in range(i+1,N):
      F_matrix[i][j] = force_LJ(separation_array[i][j],l)
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
  setup_file = sys.argv[1]
  data_file = sys.argv[2]
  traj_file = sys.argv[3] if len(sys.argv) == 4 else 'traj.xyz'
  '''
  setup_file = 'setup.dat'
  data_file = 'data.dat'
  traj_file = 'traj.xyz'

  # Get information from input files
  with open(setup_file, 'r') as setup:
    lines = setup.readlines()
    try:
      temp, mass, density, epsilon, sigma = map(float, lines[1::2]) 
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
  print(temp, mass, density, epsilon, sigma)
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

  x = p3d_list[1].pos[1]
  for particle in p3d_list:
    particle.pos /= x


  # --- SIMULATION -------------------------------------------------------

  # Open output files:

  outfile_traj = open(traj_file,'w')

  if log_total_e:
    outfile_total_e = open('total_e.csv','w')
  if log_kinetic_e:
    outfile_kinetic_e = open('kinetic_e.txt', 'w')
  if log_potential_e:
    outfile_potential_e = open('potential_e.txt', 'w')

  outfile_vel = open('velocities.txt','w')

  # Initialize objects for use in simulation
  t = 0.0
  t_step = int(t_total / dt)

  # simulation loop

  for _ in range(t_step):

    separation_array = p3d.pair_separations(p3d_list)

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
      outfile_potential_e.write(f'{t:.4} {obs.potential_energy(separation_array)}\n')

    outfile_vel.write(f'time = {t:.4}\n')
    for particle in p3d_list:
      outfile_vel.write(f'{particle.vel}\n')
      

    # --- Particle Time Integrator ----------------------------------------

    # calculate all the net forces for each particle
    F_matrix = force_matrix(separation_array,cell_length)
    
    print(f't = {t}\n') #%
    print(f'separation:\n{separation_array}')
    print(f'forces:\n{F_matrix}\n')
    

    F_list = []
    for i in range(N):
      F_i = np.zeros(3)
      for j in range(N):
        if i < j: 
          F_i += F_matrix[i][j]
        elif i > j: 
          F_i -= F_matrix[j][i]
      F_list.append(F_i)
    
    print(f'f_list:\n{F_list}') #%

    # update particle positions
    for i, particle in enumerate(p3d_list):
      particle.update_pos_2nd(cell_length,dt,F_list[i])

    # update forces
    separation_array = p3d.pair_separations(p3d_list)
    F_matrix = force_matrix(separation_array, cell_length)
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

  # --- POST SIMULATION CALCULATIONS --------------------------------

  # Total Energy
  if log_total_e:
    outfile_kinetic_e = open('kinetic_e.txt','r')
    outfile_potential_e = open('potential_e.txt','r')
    kinetic_lines = outfile_kinetic_e.readlines()
    potential_lines = outfile_potential_e.readlines()
    for i in range(t_step):
      t, kinetic_e = kinetic_lines[i].split(' ')
      potential_e = potential_lines[i].split(' ')[1]
      outfile_total_e.write(f'{t}, {float(kinetic_e)+float(potential_e)}\n')

    outfile_total_e.close()


  print('DONE ! ')










# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
