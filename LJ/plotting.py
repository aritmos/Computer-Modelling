from matplotlib import pyplot as plt
import numpy as np


data = np.genfromtxt('output\\energy_total.csv',delimiter=',',names=['t','E'])
plt.title('Total Energy')
plt.plot(data['t'],data['E'])
plt.show() 

data = np.genfromtxt('output\\msd.csv',delimiter=',',names=['t','E'])
plt.plot(data['t'],data['E'])
plt.title('MSD')
plt.show() 

data = np.genfromtxt('output\\rdf.csv',delimiter=',',names=['t','E'])
plt.plot(data['t'],data['E'])
plt.title('RDF')
plt.show()
