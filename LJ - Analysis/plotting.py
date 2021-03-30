from matplotlib import pyplot as plt
import numpy as np


def main(directory, booleans):

    LOG_TOTAL_E, LOG_MSD, LOG_RDF = booleans

    if LOG_TOTAL_E:
        dataET = np.genfromtxt(
            f'{directory}\\energy_total.csv', delimiter=' ', names=['t', 'E'])
        plt.title('Total Energy')
        plt.plot(dataET['t'], dataET['E'])
        plt.savefig(f'{directory}\\ENERGY.png')
        plt.close()

    if LOG_MSD:
        data = np.genfromtxt(
            f'{directory}\\msd.csv', delimiter=',', names=['t', 'MSD'])
        plt.plot(data['t'], data['MSD'])
        plt.title('MSD')
        plt.savefig(f'{directory}\\MSD.png')
        plt.close()

    if LOG_RDF:
        data = np.genfromtxt(
            f'{directory}\\rdf.csv', delimiter=',', names=['t', 'RDF'])
        plt.plot(data['t'], data['RDF'])
        plt.title('RDF')
        plt.savefig(f'{directory}\\RDF.png')
        plt.close()
