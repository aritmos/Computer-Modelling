"""
CMod Exercise 1: Vector class that mimics numpy's built-in
functionalities including operator overloading (except for *)
only using dunder methods.

Note: for two vectors 'u','v' and a scalar 'b':
  the scalar prodcut is defined as: b*u or u*b
  the inner product is defined as: u*v
  the cross product is defined as: u%v

Author: Sebastian Garcia
Version: 09/2020
"""
import math 

class vector():
  def __init__(self,values):
    """
    :param self: vector
    :param values: list of the components of the vector
    :attribute values: list of the components of the vector
    :attribute dim: dimension of the vector
    """
    self.values = values
    self.dim = len(self.values)
    self.norm = math.sqrt(sum(i**2 for i in self.values))

  def __add__(self,other):
    """
    :param self: vector
    :param other: other vector
    :return: a vector whose components are the sum of self and other
    """
    if self.dim==other.dim:
      new_values = [i+j for i,j in zip(self.values,other.values)]
      return vector(new_values)
    raise Exception('Dimension Mismatch')

  def __sub__(self,other):
    """
    :param self: vector
    :param other: other vector
    :return: a vector whose components are the subtraction of self and other
    """
    return self+(other*(-1))

  def __mul__(self,other):
    """
    :param self: vector
    :param other: either a scalar or a vector
    :return: if other is a scalar: the scaled vector
             if other is a vector: the inner product
    or the corresponding component of other
    """
    if isinstance(other, (int,float)):
      new_values = [other*i for i in self.values]
      return vector(new_values)
    elif isinstance(other, vector):
      if self.dim==other.dim:
        return sum([i*j for i, j in zip(self.values, other.values)])
      raise Exception('Dimension Mismatch')
    raise Exception('Type Mismatch')

  def __rmul__(self,other):
    """
    In the case of ()*vector where () isn't a vector
    run vector*() as is defined in __mul__
    """
    return (self*other)

  def __str__(self):
    """
    :param self: vector
    :return: self.values as str
    """
    return f'{self.values}'

  def __mod__(self,other):
    """
    :param self: vector
    :param other: other vector
    :return: cross product of self and other
    """
    if self.dim==other.dim==3:
      a = self.values  # rename values to a,b (arrays)
      b = other.values # for ease of readablity
      return vector([a[1]*b[2]-a[2]*b[1],
                     a[2]*b[0]-a[0]*b[2],
                     a[0]*b[1]-a[1]*b[0]])
    raise Exception('Dimension Error')

  def __neg__(self):
    return self*(-1)
