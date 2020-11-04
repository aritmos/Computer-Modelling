"""
CMod Exercise 1: A collection of functions for use in 
Periodic Boundary Condition problems for a cube of length 'l'

Author: Sebastian Garcia (s1910157)
Version: 10/2020
"""

def periodic_image(x:list,l:float) -> list:
  """
  Periodic Boundary
  :param x: point in R^3 
  :param l: cube side length 
  :return: image of the point inside the 'original' cube 
  """
  return [i % l for i in x]

def minimum_image(x:list,l:float) -> list:
  """
  Minimum Image Convention
  :param x: point in R^3
  :param l: cube length
  :return: image of the point closest to the origin
  """
  x = periodic_image(x,l)
  return [i if i <= l/2 else i-l for i in x]
