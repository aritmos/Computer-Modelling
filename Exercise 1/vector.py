"""

CMod Exercise 1: vector.py, a module that provides a set of functions to
manipulate vectors, expressed as lists.

Author: Sebastian Garcia
Version: 09/2020

"""
import math

def norm_sq(v):
  """
  Square modulus
  
  :param v: vector [v1,...,vn]
  :return: square modulus v1**2+...+vn**2
  """
  return sum([i**2 for i in v])


def norm(v):
  """
  Modulus
  
  :param v: vector [v1,...,vn]
  :return: modulus (v1**2+...+vn**2)**(1/2)
  """
  return math.sqrt(norm_sq(v))


def scalar_mult(v, scalar):
  """
  Multiplication of vector by scalar
  
  :param v: vector [v1,...,vn]
  :param scalar: scalar factor
  :return: scaled vector [v1*scalar,..., vn*scalar]
  """
  return [i*scalar for i in v]


def scalar_div(v,scalar):
  """
  Division of vector by scalar
  
  :param v: vector [v1,...,vn]
  :param scalar: scalar factor
  :return: scaled vector [v1/scalar,..., vn/scalar]
  """
  return scalar_mult(v,1/scalar)


def add(v, w):
  """
  Vector addition
  
  :param v: first vector
  :param w: second vector
  :return: vector sum v+w
  """
  if len(v)==len(w):
    return [i+j for i,j in zip(v,w)] 
  raise Exception('Dimension Mismatch')


def sub(v, w):
  """
  Vector subtraction
  
  :param v: first vector
  :param w: second vector
  :return: difference v-w
  """
  return add(v,scalar_mult(w,-1))


def dot_product(v, w):
  """
  Vector dot product
  
  :param v: first vector
  :param w: second vector
  :return: dot product v*w
  """
  if len(v) == len(w):
    return sum([i*j for i, j in zip(v, w)]) 
  raise Exception('Dimension Mismatch')


def cross_product(v, w):
  """
  Vector cross product
  
  :param v: first vector
  :param w: second vector
  :return: v x w only if both vectors have dimension 3
  """
  if len(v) == len(w) == 3:
    return [v[1]*w[2]-v[2]*w[1],-v[0]*w[2]+v[2]*w[0],v[0]*w[1]-v[1]*w[0]]
  raise Exception('Dimension Error')


def equal(v,w):
  """
  Vector equality
  :param v: first vector
  :param w: second vector
  :return: True if v==w given a specific tolerance
  """
  return all(abs(i-j)<min(abs(i)*10**-10,abs(j)*10**-10) for i,j in zip(v,w))