import sys
def main():
  # Handle arguments passed in the command line
  n = len(sys.argv)
  if n!=2:
    print(f'IndexError:\nProgram takes 1 positional argument but {n-1} were given')
    sys.exit(-1)
  N = sys.argv[1] # define the limit for fizzbuzz
  if N.isnumeric()==False:
    print('ValueError:\nPositional argument must be of type "int"')
    sys.exit(-1)
  # fizzbuzz program
  for i in range(1,int(N)+1):
    print('Fizz'*(i%3==0)+'Buzz'*(i%5==0) or i)
    # adds two strings together
    # the first string is "Fizz" if 3|i else ""
    # the second string is "Buzz" if 5|i else ""
    # if the combined string is not empty, print the string
    # else print i
if __name__ == '__main__':main()
