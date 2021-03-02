'''
from abc import abstractproperty
from lj_sim import force_matrix, net_forces
import sys
import time
import numpy as np
from particle3D import Particle3D as p3d
import mdutilities as mdu
import pbc 
from tqdm import tqdm
'''
from math import prod

def persistence(n):
    i = 0
    while len(str(n)) > 1:
      n = prod([int(i) for i in str(n)])
      i += 1
    return i

print(persistence(999))