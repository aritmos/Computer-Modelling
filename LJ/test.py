import sys
import time
import numpy as np
from particle3D import Particle3D as p3d
import mdutilities as mdu
import pbc 

l = 2
x = 1.78
res = 10
histogram = np.zeros(res)
histogram[int(x/l * res)] += 1


def update_progress(progress):
    print ('\r[{0}] {1}%'.format('#'*(progress/10), progress))

for i in range(101):
  update_progress(i)


