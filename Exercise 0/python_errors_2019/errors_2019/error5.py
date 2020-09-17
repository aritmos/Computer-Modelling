"""
Create a list of random integers, then delete the even entries,
then print the remaining entries.
"""

import random

# Create and print a list of 20 integers
# where all entries are 0<=n<=100
numbers = [random.randint(0,100) for n in range(20)]
print("Original list:")
print(numbers)

#index loops and deleting entries doesn't mix well
#for small arrays its fine to create a new filtered array
#for example through list comprehension 
odd_numbers = [i for i in numbers if i%2!=0]

print("Odd entries:")
print(odd_numbers)

