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
import vector as vec

def main():
  u,v,w = [[random.random() for _ in range(3)] for _ in range(3)]
  print(f'''
  u = {u}
  v = {v}
  w = {w}

  ||u|| = {vec.norm(u)}
  ||v|| = {vec.norm(v)}
  ||w|| = {vec.norm(w)}

  u + v = {vec.add(u,v)}
  u * v = {vec.dot_product(u,v)}
  u x v = {vec.cross_product(u,v)}

  Testing Vector Identities:
  {vec.cross_product(u,v)} = {vec.scalar_mult(vec.cross_product(v,u),-1)}
  {vec.cross_product(u,vec.add(v,w))} = {vec.add(vec.cross_product(u,v),vec.cross_product(u,w))}
  {vec.cross_product(u,vec.cross_product(v,w))} = {vec.sub(vec.scalar_mult(v,vec.dot_product(u,w)),vec.scalar_mult(w,vec.dot_product(u,v)))}

  Testing equality for the last vector identity:
  {vec.equal(vec.cross_product(u,vec.cross_product(v,w)),vec.sub(vec.scalar_mult(v,vec.dot_product(u,w)),vec.scalar_mult(w,vec.dot_product(u,v))))}

  ''')
# Execute main method, but only if it is invoked directly
if __name__ == "__main__": main()