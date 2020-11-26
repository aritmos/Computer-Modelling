import time
import frequency_finder
import velocityVerlet3D
import symplecticEuler3D

def main():
  # Specify parameters to use
  # Simulation time is fixed to 20 units
  particle_type = 'O2' # (N2 or O2)
  time_integ_short = 'SE'  # (VV or SE)
  dt_list = [0.5, 0.1, 10**-2, 10**-3, 10**-4] # dt to simulate (units in ~10.18fs)
  N = 10  # Number of simulations per dt

  # Set the correct naming for simulation and output
  if particle_type not in ('N2','O2'): print(f'No data exists on {particle_type}');quit()
  if time_integ_short not in ['VV', 'SE']:
    print(f'{time_integ_short} is shorthand for no known time integrator')
    quit()
  time_integrator_name = 'Velocity Verlet' if time_integ_short == 'VV' else 'Symplectic Euler'
  time_integrator = {'VV':velocityVerlet3D.main,'SE':symplecticEuler3D.main}
  print(f'-RUNNING SIMULATIONS FOR {particle_type}-')

  # Velocity Verlet Simulations
  print(f'Running simulations using {time_integrator_name}:')

  # N2 Simulations
  for dt in dt_list:
    print(f' Running {N} with dt = {dt}:')
    simulation_times = [] # Simulation times 
    periods = [] # Average vibration for each simulation in cm^-1
    for _ in range(N):
    
      t0 = time.time()
      time_integrator[time_integ_short](f'Input_{particle_type}.txt', dt)
      t1 = time.time()
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
      Average period for all simulations: {avg_period:.1f} cm^-1
    ''')
      

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
