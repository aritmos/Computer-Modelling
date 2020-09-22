"""
CMod Ex2: tester of the complex.py module.
Creates a pair of random complex numbers and
runs through all functions implemented in complex.py

Author: A. Hermann
Version: 06/2017
"""

import math
import random
import complex as cplx

# Main method:
def main():
    # Create two random complex numbers

    c1 = [random.random(), random.random()]
    c2 = [random.random(), random.random()]

    # Print out both numbers
    print("c1 = ", c1)
    print("c2 = ", c2)

    # Test unary functions from complex.py:
    # conjugation, modulus (squared), scaling
    print("conj(c2) = ", cplx.conj(c2))
    print("|c2|^2 = ", cplx.norm_sq(c2))
    print("|c2|   = ", cplx.norm(c2))
    print("3*c2   = ", cplx.scale(c2,3.0))

    # Test binary functions from complex.py:
    # addition, subtraction, multiplication, division

    print("c1+c2 = ", cplx.add(c1,c2))
    print("c1-c2 = ", cplx.sub(c1,c2))
    print("c1*c2 = ", cplx.mul(c1,c2))
    print("c1/c2 = ", cplx.div(c1,c2))


# Execute main method, but only if it is invoked directly
if __name__ == "__main__":
    main()
