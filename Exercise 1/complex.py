"""

CMod Ex. 2: complex.py, a module that provides a set of functions to
manipulate complex numbers, expressed as lists.

Author: A. Hermann
Version: 06/2017

"""
import math


def conj(c):
    """
    Complex conjugate

    :param c: complex number (c1, c2)
    :return: complex conjugate (c1, -c2)
    """
    return [c[0], -c[1]]


def norm_sq(c):
    """
    Square modulus
    
    :param c: complex number (c1, c2)
    :return: square modulus c1**2+c2**2
    """
    return c[0]**2 + c[1]**2


def norm(c):
    """
    Modulus
    
    :param c: complex number (c1, c2)
    :return: modulus (c1**2+c2**2)**(1/2)
    """
    return math.sqrt(norm_sq(c))


def scale(c, scalar):
    """
    Multiplication of complex with scalar
    
    :param c: complex number (c1, c2)
    :param scalar: scalar factor
    :return: scaled complex number (c1*scalar, c2*scalar)
    """
    return [c[0]*scalar, c[1]*scalar]


def add(c1, c2):
    """
    Complex addition
    
    :param c1: first complex number
    :param c2: second complex number
    :return: complex sum c1+c2
    """
    return [c1[0]+c2[0], c1[1]+c2[1]]


def sub(c1, c2):
    """
    Complex subtraction
    
    :param c1: first complex number
    :param c2: second complex number
    :return: complex difference c1-c2
    """
    return [c1[0]-c2[0], c1[1]-c2[1]]


def mul(c1, c2):
    """
    Complex multiplication
    
    :param c1: first complex factor
    :param c2: second complex factor
    :return: complex product c1*c2
    """
    real_part = c1[0]*c2[0] - c1[1]*c2[1]
    imag_part = c1[0]*c2[1] + c1[1]*c2[0]
    return [real_part, imag_part]


def div(c1, c2):
    """
    Complex division
    
    :param c1: complex dividend
    :param c2: complex divisor
    :return: c1/c2 = c1*conj(c2)/|c2|**2
    """
    z = mul(c1,conj(c2))
    return scale(z,1.0/norm_sq(c2))


