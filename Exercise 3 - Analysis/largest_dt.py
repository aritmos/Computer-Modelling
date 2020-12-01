import single_simulation_data_analyzer as ssda

dt = 0.05
input_file = 'N2.txt'
time_integrator_short = 'VV'


freq = ssda.main(input_file,dt,time_integrator_short)
v0 = 2140 # 1433.0 or 2350.4 
while 0.995*v0<freq<1.005*v0:
  dt+=0.001
  freq = ssda.main(input_file, dt, time_integrator_short)

