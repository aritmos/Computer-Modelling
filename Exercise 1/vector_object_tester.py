"""
CMod Exercise 1: Tester of the vector.py module.
Creates a triplet of random vectors (in R^3) then
runs through functions implemented in vector.py and
checks vector identities

Author: Sebastian Garcia
Version: 09/2020
"""
import math
import random
import vector_object as vec

def main():
  u, v, w = [vec.vector([random.random() for _ in range(3)]) for _ in range(3)]
  print(f'''
  u = {u}
  v = {v}
  w = {w}

  ||u|| = {u.norm}
  ||v|| = {v.norm}
  ||w|| = {w.norm}

  u + v = {u+v}
  u * v = {u*v}
  u x v = {u%v}

  Testing Vector Identities:
  {u%v} = {-v%u}
  {u%(v+w)} = {u%v+u%w}
  {u%(v%w)} = {(u*w)*v-(u*v)*w}
  ''')


# Execute main method, but only if it is invoked directly
if __name__ == "__main__":
  main()
