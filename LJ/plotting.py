"""
CompMod Final Project: (auxiliary module) Plotting, used for plotting
the calculated observables: Total energy, MDF and RDF.

author: -redacted-

"""

from matplotlib import pyplot as plt
import numpy as np


def plot(LOG_TOTAL_E, LOG_MSD, LOG_RDF):

    if LOG_TOTAL_E:
        dataET = np.genfromtxt('output\\energy_total.csv',
                               delimiter=' ', names=['t', 'E'])
        plt.title('Total Energy')
        plt.plot(dataET['t'], dataET['E'])
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'Energy ($\epsilon$)')
        plt.show()

    if LOG_MSD:
        data = np.genfromtxt(
            'output\\msd.csv', delimiter=',', names=['t', 'MSD'])
        plt.plot(data['t'], data['MSD'])
        plt.title('Mean Squared Displacement')
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'MSD ($\sigma^2$)')
        plt.show()

    if LOG_RDF:
        data = np.genfromtxt(
            'output\\rdf.csv', delimiter=',', names=['t', 'RDF'])
        plt.plot(data['t'], data['RDF'])
        plt.title('Radial Distribution Function')
        plt.xlabel(r'Distance ($r$)')
        plt.ylabel(r'RDF ($r$)')
        plt.show()
