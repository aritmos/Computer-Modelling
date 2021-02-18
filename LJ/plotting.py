from matplotlib import pyplot as plt
import numpy as np

def main(log_total_e,log_msd,log_rdf):

  if log_total_e:
    dataET = np.genfromtxt('output\\energy_total.csv',delimiter=' ',names=['t','E'])
    dataEK = np.genfromtxt('output\\energy_kinetic.txt',
                           delimiter=' ', names=['t', 'E_K'])
    dataEP = np.genfromtxt('output\\energy_potential.txt',
                           delimiter=' ', names=['t', 'E_P'])
    plt.title('Energies')
    
    #plt.plot(dataET['t'],dataET['E'])
    plt.plot(dataEK['t'], dataEK['E_K']-dataEK['E_K'][0])
    plt.plot(dataEP['t'], -(dataEP['E_P']-dataEP['E_P'][0]))
    plt.show()


  if log_msd:
    data = np.genfromtxt('output\\msd.csv',delimiter=',',names=['t','MSD'])
    plt.plot(data['t'],data['MSD'])
    plt.title('MSD')
    plt.show() 

  if log_rdf:
    data = np.genfromtxt('output\\rdf.csv',delimiter=',',names=['t','RDF'])
    plt.plot(data['t'],data['RDF'])
    plt.title('RDF')
    plt.show()

main(True,False,False)
