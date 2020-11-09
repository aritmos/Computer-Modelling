"""
 CompMod Ex2: Auxiliary tester program for the 
 new_particle method in the Particle3D class

author: Sebastian Garcia (s1910157)

"""

from particle3D import Particle3D as p3d

with open('particles.txt','r') as input_file: #automatically closes the file
  p3d_list = p3d.new_particle(input_file)

for i in p3d_list:print(i)
print(p3d.com_velocity(p3d_list))
