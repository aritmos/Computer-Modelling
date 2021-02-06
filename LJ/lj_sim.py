
# IMPORTS
import math
from os import sys
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
  pair_separation = np.array(pbc.minimum_image(pair_separation,l))
  r = np.linalg.norm(pair_separation)
  F = -48*(r**(-12)-r**(-6))*r**(-1)*pair_separation
  return F


def force_i(N:int, i:int, l:float, separation_array:np.ndarray) -> np.array:
  '''
  Calculates the net force on particle i
  * requires l to pass to force_LJ to use MIC
  '''
  F = np.zeros(3)
  for j in range(N):
    if i < j:
      F += force_LJ(separation_array[i][j],l)
    elif i > j:
      # separation_array only carries forces i < j
      # -> for i > j we subtract the force of particle j due to i 
      F -= force_LJ(separation_array[j][i],l)
  return F


def main():

  #--- SETTING VARIABLES AND PARAMETERS --- 

  print('---\nReading from input files...')
  # Read needed information from command line
  if len(sys.argv) not in [3,4]:
    print("Wrong number of arguments")
    print(f"Usage: {sys.argv[0]} <setup file> <data file> <OPTIONAL: output file>")
    quit()
  setup_file = sys.argv[1]
  data_file = sys.argv[2]
  traj_file = sys.argv[3] if len(sys.argv) == 4 else 'traj.xyz'

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

  # --- INITIALISATION --- 

  print(f'Initialising {N} particles:')
  # Create the N particles with dummy pos and vel
  p3d_list = p3d.new_particles(N, mass)
  
  # Set initial conditions for particles
  cell_length, full_lattice = mdu.set_initial_positions(density, p3d_list)
  mdu.set_initial_velocities(temp, p3d_list)

  '''
  # debugging
  print(cell_length, full_lattice)
  for i in p3d_list:
    print(i)
  '''

  # --- SIMULATION ---

  # Open output files:

  outfile_traj = open(traj_file,'w')

  # Initialize objects for use in simulation
  t = 0

  # simulation loop

  t_step = int(t_total / dt)
  for step in range(t_step):

    # --- Logging Observables ---

    # Trajectories:
    outfile_traj.write(f'{N}\n')
    outfile_traj.write(f'Step = {step}\n')
    for particle in p3d_list:
      outfile_traj.write(f'{particle}\n')

    # --- Particle Time Integrator ---

    # calculate all the net forces for each particle
    separation_array = p3d.pair_separations(p3d_list)
    force_list = np.zeros([N,3])
    for i in range(N):
      force_list[i] = force_i(N, i, cell_length, separation_array)

    # update particle positions
    for i, particle in enumerate(p3d_list):
      particle.update_pos_2nd(cell_length,dt,force_list[i])

    # update forces
    separation_array = p3d.pair_separations(p3d_list)
    force_new_list = np.zeros([N,3])
    for i in range(N):
      force_new_list[i] = force_i(N, i, cell_length, separation_array)

    # update particle velocities
    for i, particle in enumerate(p3d_list):
      particle.update_vel(dt,0.5*(force_list[i]+force_new_list[i]))

    # update time
    t += dt
    


  outfile_traj.close()










# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
