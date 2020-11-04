"""
CMod Exercise 1: Tester of the numpy library as 
an alternative to vector.py
Creates a triplet of random vectors (in R^3) then
runs through functions implemented in numpy and
checks vector identities

Author: Sebastian Garcia (s1910157)
Version: 09/2020
"""
import numpy as np

def main():
  [u,v,w] = np.random.rand(3,3)

  print(f'''
  u = {u}
  v = {v}
  w = {w}

  ||u|| = {np.linalg.norm(u)}
  ||v|| = {np.linalg.norm(v)}
  ||w|| = {np.linalg.norm(w)}

  u + v = {u+v}
  u * v = {np.inner(u,v)}
  u x v = {np.cross(u,v)}

  Testing Vector Identities:
  {np.cross(u,v)} = {np.cross(-v,u)}
  {np.cross(u,v+w)} = {np.cross(u,v)+np.cross(u,w)}
  {np.cross(u,np.cross(v,w))} = {np.inner(u,w)*v-np.inner(u,v)*w}

  Testing equality for the last vector identity:
  {np.allclose(np.cross(u,np.cross(v,w)),np.inner(u,w)*v-np.inner(u,v)*w)}
  ''')
  # default arguments (rtol & atol) for np.allclose worked well so they arent passed

# Execute main method, but only if it is invoked directly
if __name__ == "__main__":
  main()
