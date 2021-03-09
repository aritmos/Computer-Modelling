from matplotlib import pyplot as plt
import numpy as np


def main(LOG_TOTAL_E, LOG_MSD, LOG_RDF):

    if LOG_TOTAL_E:
        dataET = np.genfromtxt('output\\energy_total.csv',
                               delimiter=' ', names=['t', 'E'])
        plt.title('Total Energy')
        plt.plot(dataET['t'], dataET['E'])
        plt.show()

    if LOG_MSD:
        data = np.genfromtxt(
            'output\\msd.csv', delimiter=',', names=['t', 'MSD'])
        plt.plot(data['t'], data['MSD'])
        plt.title('MSD')
        plt.show()

    if LOG_RDF:
        data = np.genfromtxt(
            'output\\rdf.csv', delimiter=',', names=['t', 'RDF'])
        plt.plot(data['t'], data['RDF'])
        plt.title('RDF')
        plt.show()


# main(True, False, False)
