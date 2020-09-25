import particle3D as p3d
import numpy as np


p3d_list = []

with open('particles.txt','r') as f:
  lines = f.readlines()
  for line in lines:
    p3d_list.append(p3d.Particle3D.new_particle(line))

print(p3d.Particle3D.com_velocity(p3d_list))
