import numpy as np

def potential_energy(separation_array:np.ndarray) -> float:
  total_potential_energy = 0
  N = len(separation_array)
  for i in range(N-1):
    for j in range(i+1,N):
      r = np.linalg.norm(separation_array[i][j])
      total_potential_energy += 4*(r**(-12)-r**(-6))
  return total_potential_energy

