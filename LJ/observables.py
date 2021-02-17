import numpy as np
import pbc 

def potential_energy(separation_array:np.ndarray) -> float:
  total_potential_energy = 0
  N = len(separation_array)
  for i in range(N-1):
    for j in range(i+1,N):
      r = np.linalg.norm(separation_array[i][j])
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
  # Maximum MIC distance is l*sqrt(3)/2,
  # due to floating point rounding in calculations we add 0.001*l
  # Max distance will be the range of the histogram
  # the histogram will have (res) number of bars
  max_distance = l*(np.sqrt(3)/2+0.001)

  for i in range(1,N):
    # Reference particle is p0
    separation = particle_separations[0][i]
    mic_separation = pbc.minimum_image(separation,l)
    mic_separation_distance = np.linalg.norm(mic_separation)
    histogram[int((mic_separation_distance/ max_distance) * res)] += 1

  return histogram

def rdf_normalize(histogram:np.array,l:float, res:int) -> np.array:
  max_length = l*(np.sqrt(3)/2+0.001)
  area = sum(histogram)*(max_length/res)
  return histogram/area

