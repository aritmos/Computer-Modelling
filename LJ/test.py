from abc import abstractproperty
from lj_sim import force_matrix, net_forces
import sys
import time
import numpy as np
from particle3D import Particle3D as p3d
import mdutilities as mdu
import pbc 
from tqdm import tqdm

l = 2
x = 1.78
res = 10
histogram = np.zeros(res)
histogram[int(x/l * res)] += 1

for i,j in range(10):
  print(i,j)

mic_matrix[i][j] = pair