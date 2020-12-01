import velocityVerlet3D as vv
import symplecticEuler3D as se

input_file = 'N2.txt'

out = []
for dt in [0.5,0.4, 0.3, 0.2, 0.1]:
  se.main(input_file,dt)
  data  = open('energy.txt','r')
  lines = data.readlines()
  E_0 = float(lines[0].split(' ')[1])
  energy_error_list = []
  for line in lines:
    Et = float(line.split(' ')[1])
    energy_error_list.append(abs((E_0-Et)/E_0))
  out.append(max(energy_error_list))
print(out)
