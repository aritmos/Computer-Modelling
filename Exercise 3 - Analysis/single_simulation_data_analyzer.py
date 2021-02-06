import time
import frequency_finder
import velocityVerlet3D
import symplecticEuler3D

def main(input_file:str, dt:float, time_integ_short:str):
  # Specify parameters to use
  # Simulation time is fixed to 20 units
  '''
  input_file = 'O2.txt' 
  time_integ_short = 'SE'  # (VV or SE)
  '''
  N = 1  # Number of simulations per dt


  # Set the correct naming for simulation and output
  if input_file not in ('N2.txt','O2.txt','N2_roto.txt','O2_roto.txt'): print(f'No data exists on {input_file}');quit()
  if time_integ_short not in ['VV', 'SE']:
    print(f'{time_integ_short} is shorthand for no known time integrator')
    quit()
  time_integrator_name = 'Velocity Verlet' if time_integ_short == 'VV' else 'Symplectic Euler'
  dictionary = {'VV':velocityVerlet3D.main,'SE':symplecticEuler3D.main}
  print(
      f'-RUNNING SIMULATIONS FOR {input_file[:-4]} using {time_integrator_name}-')
  freq_list = []
  # Simulations
  print(f' Running {N} simulations with dt = {dt}:')
  simulation_times = [] # Simulation times 
  periods = [] # Average vibration for each simulation in cm^-1
  t0 = time.time()
  dictionary[time_integ_short](input_file, dt)
  t1 = time.time()
  simulation_time = t1 - t0
  simulation_times.append(simulation_time)
  period = frequency_finder.main()
  periods.append(period)
  avg_time = sum(simulation_times)/N
  min_time = min(simulation_times)
  avg_period = sum(periods)/N
  freq_list.append(f'{avg_period:.1f}')
  print(f'''
    Average simulation compute time: {avg_time:.4f}
    Minimum simulation compute time: {min_time:.4f}
    Average vibrational period: {avg_period:.1f} cm^-1
  ''')
  return avg_period
  #print(freq_list)
