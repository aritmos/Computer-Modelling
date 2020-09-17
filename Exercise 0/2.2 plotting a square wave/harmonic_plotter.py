"""
Simple Python code that plots the cosine function
"""

# Import relevant python modules
import math as m
import matplotlib.pyplot as pyplot
import sys

# Define the cosine function
def my_function(N,x):
    out = 0
    for i in range(1,N+1):
        out += ((2*i-1)**-1)*m.sin((2*i-1)*x)
    return out

# Main method
def main():
    # make sure correct parameters are being passed
    # these should be an integer and a path to the output file
    n = len(sys.argv)
    if n!=3:
        raise TypeError(f'Program takes 2 positional arguments but {n-1} were given')
    if sys.argv[1].isnumeric()==False or sys.argv[2].count('.')!=1:
        raise TypeError("""
    Program takes an 2 positional arguments
    these are of type "int" and "str" respectively
    where the "str" parameter must be a path""")

    N = int(sys.argv[1])
    # number of data points
    n_loop = 100

    # open output file
    out_file = open(sys.argv[2],"w")

    # prepare data lists
    x_values = []
    y_values = []

    # obtain function values and write them to file
    for i in range(n_loop):
        x = 2*m.pi*i/n_loop - m.pi
        f = my_function(N,x)
    
        # append data to lists and output file
        x_values.append(x)
        y_values.append(f)
    
        out_file.write(str(x) + " " + str(f) + "\n")

    # close output file
    out_file.close()

    # plot result
    pyplot.plot(x_values,y_values)
    pyplot.suptitle('Plotting the cosine function')
    pyplot.xlabel('X')
    pyplot.ylabel('Cos(X)')
    pyplot.show()


# Execute main method
if __name__ == "__main__": main()
