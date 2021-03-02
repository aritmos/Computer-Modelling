import numpy as np
import pbc 

def potential_energy(mic_separations:np.ndarray) -> float:
  total_potential_energy = 0
  N = len(mic_separations)
  for i in range(N-1):
    for j in range(i+1,N):
      r = np.linalg.norm(mic_separations[i][j])
      if r > 3.5: # Implement cutoff radius
        total_potential_energy -= 0.00217
      else:
        total_potential_energy += 4*(r**(-12)-r**(-6))
  return total_potential_energy

def msd(mic_separations:np.ndarray) -> float:
  N = len(mic_separations)
  mic_distances = np.linalg.norm(mic_separations, axis=1)
  return sum(mic_distances**2)/N

def rdf(mic_separations:np.ndarray, max_distance:float, res:int) -> np.array:
  N = len(mic_separations)
  histogram = np.zeros(res)
  mic_distances = np.linalg.norm(mic_separations, axis=2)
  for i in range(N):
    for j in range(N):
      if i > j:
        i, j = j, i
      if i != j:
        r = mic_distances[i][j]
        index = int((r / max_distance) * res)
        histogram[index] += r**-2
  return histogram

def rdf_normalize(hist:np.array, t_count:float, dr:float, N:int, rho:float) -> np.array:
  C = (N*4*np.pi*rho*dr)**-1
  return C*hist/t_count

