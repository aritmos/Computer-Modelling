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

def rdf(particle_separations:np.ndarray, l:float, res:int) -> np.array:
  N = len(particle_separations)
  histogram = np.zeros(res)
  max_distance = l*(np.sqrt(3)/2+0.001)
  for i in range(1,N):
    # Reference particle is p0
    separation = particle_separations[0][i]
    mic_separation = pbc.minimum_image(separation,l)
    r = np.linalg.norm(mic_separation)
    value = (N*4*np.pi*(r**2)*(max_distance/res))**-1
    histogram[int((r / max_distance) * res)] += value

  return histogram

def rdf_normalize(histogram:np.array, total_t:int) -> np.array:
  return histogram/(total_t/5)

