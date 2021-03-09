"""
CompMod Final Project: (auxiliary module) PBC, a collection of functions for use in 
Periodic Boundary Condition problems for a cube of length 'l'

Author: -redacted-

"""
import numpy as np


def periodic_image(x: np.array, l: float) -> np.array:
    """
    Periodic Boundary
    :param x: point in R^3 
    :param l: cube side length 
    :return: image of the point inside the 'original' cube 
    """
    return np.array([i % l for i in x.tolist()])


def minimum_image(x: np.array, l: float) -> np.array:
    """
    Minimum Image Convention
    :param x: point in R^3
    :param l: cube side length
    :return: image of the point closest to the origin
    """
    x = periodic_image(x, l)
    return np.array([i if i <= l/2 else i-l for i in x])
