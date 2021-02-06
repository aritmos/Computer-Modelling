import numpy as np
from particle3D import Particle3D as p3d
import mdutilities as mdu

N = 5

def func(x:int)-> int:
  global N
  return x*2*N

print(func(1))