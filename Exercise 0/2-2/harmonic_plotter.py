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
        print(f'IndexError:\nProgram takes 2 positional arguments but {n-1} were given')
        sys.exit(-1)
    if sys.argv[1].isnumeric()==False or sys.argv[2].count('.')!=1:
        print("""
    TypeError:
    Program takes 2 positional arguments
    of type "int" and "str" respectively.
    The "str" parameter must be a path""")
        sys.exit(-1)

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
