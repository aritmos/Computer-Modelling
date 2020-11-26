import sys
import time
import frequency_finder
from particle3D import Particle3D

def main():
  dt_list = [0.1]
  N = 10 # Number of simulations per dt

  # Symplectic Euler Simulations

  # N2 Simulations
  for dt in dt_list:
    print(f'Running {N} simulations for N2 particles with dt = {dt} ...')
    simulation_times = [] # Simulation times 
    periods = [] # Average vibration for each simulation in cm^-1
    for _ in range(N):
      script = open('symplecticEuler3D.py')
      script_read = script.read()
      sys.argv = ['symplecticEuler3D.py', 'input_N2.txt', 'N2']
      t0 = time.time()
      exec(script_read)
      t1 = time.time()
      script.close()
      simulation_time = t1 - t0
      simulation_times.append(simulation_time)
      period = frequency_finder.main()
      periods.append(period)
    avg_time = sum(simulation_times)/N
    min_time = min(simulation_times)
    avg_period = sum(periods)/N
    print(f'''
    Results for N2 with dt = {dt}:
    Average simulation compute time: {avg_time:.4f}
    Minimum simulation compute time: {min_time:.4f}
    Average period for all simulations: {avg_period} cm^-1
    ''')
      

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
