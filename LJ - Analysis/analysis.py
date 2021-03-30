from matplotlib import pyplot as plt
import numpy as np

# KINETIC TEMPERATURE OF MELTING STUFF

temps = [1.0, 1.1, 1.2, 1.25, 1.275, 1.28, 1.29, 1.3]

def str_float(s): return str(s).replace('.', '_')

'''
print('Temp : Kinetic Energy')

for T in temps:
    with open(f'output\\32-{str_float(T)}-0_85-{30 if T > 1.0 else 10}.0-0_0025\\energy_kinetic.txt', 'r') as data:
        lines = data.readlines()
        average_final_ke = sum([float(i.split(' ')[1])
                                for i in lines[::-1][:4000]])/4000
        print(f'{T:.3f} : {average_final_ke:.4f}')
'''

# TIME STEP CONVERGENCE

directory = [f'output\\30-1_0-0_05-']

dataET = np.genfromtxt(
    f'{directory}\\energy_total.csv', delimiter=' ', names=['t', 'E'])
plt.title('Total Energy')
plt.plot(dataET['t'], dataET['E'])
plt.savefig(f'{directory}\\ENERGY.png')
plt.close()
