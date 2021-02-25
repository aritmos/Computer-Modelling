import numpy as np
import pbc 

def potential_energy(separation_array:np.ndarray, l:float) -> float:
  total_potential_energy = 0
  N = len(separation_array)
  for i in range(N-1):
    for j in range(i+1,N):
      mic_distance = pbc.minimum_image(separation_array[i][j],l)
      r = np.linalg.norm(mic_distance)
      total_potential_energy += 4*(r**(-12)-r**(-6))
  return total_potential_energy

def msd(initial_positions:np.array, current_positions:np.array, l:float) -> float:
  N = len(initial_positions)
  separations = initial_positions - current_positions
  mic_separations = np.array([pbc.minimum_image(vector,l) for vector in separations])
  mic_separation_distances = np.linalg.norm(mic_separations, axis=1)
  return sum(mic_separation_distances**2)/N

def rdf(particle_separations:np.ndarray, rho:float, l:float, res:int) -> np.array:
  '''
  Creates a histogram with the partice distances from p0
  the histogram has range [0,l*sqrt(3)/2+dl]
  and <res> number of columns
  the values are r^-2, the remaining constant is added in the normalisation
  '''
  N = len(particle_separations)
  histogram = np.zeros(res)
  max_distance = l*(np.sqrt(3)/2+0.001)
  for i,j in zip(range(N),range(N)):
    if i > j:
      i, j = j, i
    if i != j:
      separation = particle_separations[i][j]
      mic_separation = pbc.minimum_image(separation,l)
      r = np.linalg.norm(mic_separation)
      index = int((r / max_distance) * res)
      histogram[index] += r**-2
  return histogram

def rdf_normalize(hist:np.array, t_count:float, dr:float, N:int, rho:float) -> np.array:
  C = (4*N*np.pi*rho*dr)**-1
  return C*hist/t_count

