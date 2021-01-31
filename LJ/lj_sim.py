
# IMPORTS
import math
import numpy as np 
import sys 

import mdutilities as mdu
import observables as obs 
from particle3D import Particle3D as p3d 
import pbc 



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
      dt, total_t = float(lines[3]),float(lines[5])

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
  print(N, dt,total_t)
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

  









# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
