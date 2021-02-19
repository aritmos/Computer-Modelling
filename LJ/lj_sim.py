
# IMPORTS
from os import sep, sys
import numpy as np 
import sys 

import mdutilities as mdu
import observables as obs 
from particle3D import Particle3D as p3d 
import pbc 
from tqdm import tqdm


def force_LJ(pair_separation:np.array, l:float) -> np.array:
  '''
  Calculates the MIC-force (on particle i due to particle j)
  given their real separation vector r_ij
  '''
  mic_pair_separation = pbc.minimum_image(pair_separation,l)
  r = np.linalg.norm(mic_pair_separation)
  F = -48*(r**(-14)-r**(-8)/2)*mic_pair_separation # -ve due to using r_ij in the formula
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


def net_forces(F_matrix: np.ndarray) -> np.ndarray:
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


def main():

  #--- SETTING VARIABLES AND PARAMETERS ---------------------------------- 

  print('---\nReading from input files...')
  # Read needed information from command line
  
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
  traj_file = 'output\\traj.xyz'
  '''

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
      dt = float(lines[3])
      t_total = float(lines[5])

      # map '0' -> 0 -> False and '1' -> 1 -> True
      LOG_KINETIC_E, LOG_POTENTIAL_E, LOG_TOTAL_E, LOG_MSD, LOG_RDF = map(
          bool, map(int, lines[7:16:2])) 
      t_multiple = float(lines[17])
      res = int(lines[19])
    except:
      print(f'error in {data_file}\ncheck all variables are present and their types are correct')
      quit()

    if LOG_TOTAL_E == True: # Override energies
      LOG_KINETIC_E = LOG_POTENTIAL_E = True
      
  
  # check that all the required variables are positive
  if not all([N, dt, t_total, t_multiple, res]) > 0:
    print(f'error: tried to initialize a parameter with a non-positive value')
    quit()

  print(' all variables have been set correctly')


  # --- INITIALISATION --------------------------------------------------------

  print(f'Initialising {N} particles:')
  # Create the N particles with dummy pos and vel
  p3d_list = p3d.new_particles(N, mass)
  
  # Set initial conditions for particles
  cell_length = mdu.set_initial_positions(density, p3d_list)[0]
  mdu.set_initial_velocities(temp, p3d_list)

  # --- SIMULATION ------------------------------------------------------------

  # Open output files and intialize content for observables:

  outfile_traj = open(traj_file,'w')

  if LOG_TOTAL_E:
    outfile_total_e = open('output\\energy_total.csv','w')
  if LOG_KINETIC_E:
    outfile_kinetic_e = open('output\\energy_kinetic.txt', 'w')
  if LOG_POTENTIAL_E:
    outfile_potential_e = open('output\\energy_potential.txt', 'w')
  if LOG_MSD:
    initial_positions = np.array([particle.pos for particle in p3d_list])
    outfile_msd = open('output\\msd.csv','w')
  if LOG_RDF:
    histogram = np.zeros(res)
    outfile_rdf = open('output\\rdf.csv','w')

  # Initialize remaining objects for use in simulation

  t = 0.0
  total_steps = int(t_total / dt)

  separation_array = p3d.pair_separations(p3d_list)
  F_matrix = force_matrix(separation_array, cell_length)
  F_new_list = net_forces(F_matrix)

  # Print what gets logged and calculated
  print('Logged/Calculated Observables:')
  print('-Particle Trajectories')
  if LOG_KINETIC_E: print('-Kinetic Energy')
  if LOG_POTENTIAL_E: print('-Potential Energy')
  if LOG_TOTAL_E: print('-Total Energy')
  if LOG_MSD:
    print(
        f'-Mean Squared Displacement (Every {"" if t_multiple == 1 else t_multiple} step{"s" if t_multiple > 1 else ""}')
  if LOG_RDF: 
    print(
      f'-Radial Distribution Function (Every {"" if t_multiple == 1 else t_multiple} step{"s" if t_multiple > 1 else ""}, with a resolution of {res})')

  # Simulation Loop ===========================================================
  print('\nSimulating...')

  for time_step in tqdm(range(total_steps)):

    # --- Logging Observables -------------------------------------------------

    # Trajectories:
    outfile_traj.write(f'{N}\n')
    outfile_traj.write(f'time = {t:.4f}\n')
    for particle in p3d_list:
      outfile_traj.write(f'{particle}\n')

    # Energies:
    sys_kinetic = p3d.sys_kinetic(p3d_list)
    sys_potential = obs.potential_energy(separation_array, cell_length)

    if LOG_KINETIC_E:
      outfile_kinetic_e.write(f'{t:.3} {sys_kinetic:.4f}\n')

    if LOG_POTENTIAL_E:
      outfile_potential_e.write(f'{t:.3f} {sys_potential:.4f}\n')

    if LOG_TOTAL_E:
      outfile_total_e.write(f'{t:.3f} {(sys_kinetic+sys_potential):.4f}\n')

    # Mean Squared Displacement:
    if LOG_MSD and (time_step % t_multiple) == 0:
      particle_positions = np.array([particle.pos for particle in p3d_list])
      outfile_msd.write(f'{t:.3f}, {obs.msd(initial_positions, particle_positions, cell_length):.5f}\n')

    # Radial Distribution Function:
    if LOG_RDF and (time_step % t_multiple) == 0:
      histogram += obs.rdf(separation_array, cell_length, res)

    # --- Particle Time Integrator --------------------------------------------

    # update particle positions
    F_list = F_new_list
    for i, particle in enumerate(p3d_list):
      particle.update_pos_2nd(cell_length,dt,F_list[i])

    # update forces
    separation_array = p3d.pair_separations(p3d_list)
    F_matrix = force_matrix(separation_array, cell_length)
    F_new_list = net_forces(F_matrix)

    # update particle velocities
    for i, particle in enumerate(p3d_list):
      F_avg = 0.5*(F_list+F_new_list)
      particle.update_vel(dt,F_avg[i])

    # update time
    t += dt
  # ===========================================================================

  # Close all simulation files
  outfile_traj.close()
  if LOG_KINETIC_E: outfile_kinetic_e.close()
  if LOG_POTENTIAL_E: outfile_potential_e.close()
  if LOG_TOTAL_E: outfile_total_e.close()
  if LOG_MSD: outfile_msd.close()

  # --- POST SIMULATION CALCULATIONS ------------------------------------------


  # Radial Distribution Function
  if LOG_RDF:
    histogram = obs.rdf_normalize(histogram, t_total)
    max_length = cell_length*(np.sqrt(3)/2+0.01)
    for i in range(res):
      outfile_rdf.write(f'{i*max_length/res:.4f}, {histogram[i]}\n')
    outfile_rdf.close()

  print(' done.')

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
