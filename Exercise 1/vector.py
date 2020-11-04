"""
CMod Exercise 1: vector.py, a module that provides a set of functions to
manipulate vectors, expressed as lists of their components.

Author: Sebastian Garcia (s1910157)
Version: 10/2020
"""
import math as m

def norm_sq(v: list) -> float:
  """
  Square modulus
  
  :param v: vector [v1,...,vn]
  :return: square modulus v1**2+...+vn**2
  """
  return sum([i**2 for i in v])


def norm(v: list) -> float:
  """
  Modulus
  
  :param v: vector [v1,...,vn]
  :return: modulus (v1**2+...+vn**2)**(1/2)
  """
  return m.sqrt(norm_sq(v))


def scalar_mult(v: list, scalar: float) -> list:
  """
  Multiplication of vector by scalar
  
  :param v: vector [v1,...,vn]
  :param scalar: scalar factor
  :return: scaled vector [v1*scalar,..., vn*scalar]
  """
  return [i*scalar for i in v]


def scalar_div(v: list,scalar: float) -> list:
  """
  Division of vector by scalar
  
  :param v: vector [v1,...,vn]
  :param scalar: scalar factor
  :return: scaled vector [v1/scalar,..., vn/scalar]
  """
  return scalar_mult(v,1/scalar)


def add(v: list, w: list) -> list:
  """
  Vector addition
  
  :param v: first vector
  :param w: second vector
  :return: vector sum v+w
  """
  if len(v)==len(w):
    return [i+j for i,j in zip(v,w)] 
  raise Exception('Dimension Mismatch')


def sub(v: list, w: list) -> list:
  """
  Vector subtraction
  
  :param v: first vector
  :param w: second vector
  :return: difference v-w
  """
  return add(v,scalar_mult(w,-1))


def dot_product(v: list, w: list) -> float:
  """
  Vector dot product
  
  :param v: first vector
  :param w: second vector
  :return: dot product v*w
  """
  if len(v) == len(w):
    return sum([i*j for i, j in zip(v, w)]) 
  raise Exception('Dimension Mismatch')


def cross_product(v: list, w:list) -> list:
  """
  Vector cross product
  
  :param v: first vector
  :param w: second vector
  :return: v x w only if both vectors have dimension 3
  """
  if len(v) == len(w) == 3:
    return [v[1]*w[2]-v[2]*w[1],-v[0]*w[2]+v[2]*w[0],v[0]*w[1]-v[1]*w[0]]
  raise Exception('Dimension Error')


def equal(v: list,w: list) -> bool:
  """
  Vector equality (with tolerance)

  :param v: first vector
  :param w: second vector
  :return: v==w given a specific tolerance
  """
  return all(abs(i-j)<min(abs(i),abs(j))*10**-10 for i,j in zip(v,w))