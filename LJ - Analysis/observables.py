import numpy as np
import pbc


def potential_energy(mic_separations: np.ndarray) -> float:
    """
    Calculates the total potential energy of a system of N particles
    interacting via a LJ potential

    :param mic_separations: [N,N,3] matrix of particle separations
    :return: total potential energy
    """
    total_potential_energy = 0
    N = len(mic_separations)
    for i in range(N-1):
        for j in range(i+1, N):
            r = np.linalg.norm(mic_separations[i][j])
            if r > 3.5:  # Implement cutoff radius
                total_potential_energy -= 0.00217
            else:
                total_potential_energy += 4*(r**(-12)-r**(-6))
    return total_potential_energy


def msd(initial_positions: np.ndarray, p3d_list: 'list[Particle3D]', l: float) -> float:
    """
    Calculates the MSD of a system of N particles

    :param mic_separations: [N,3] matrix of the displacement of the particles
    from their initial positions
    :return: Total mean square displacement
    """
    N = len(initial_positions)
    current_positions = np.array([particle.pos for particle in p3d_list])
    displacements = current_positions - initial_positions
    mic_displacements = np.array([pbc.minimum_image(i, l)
                                  for i in displacements])
    mic_distances = np.linalg.norm(mic_displacements, axis=1)
    return sum(mic_distances**2)/N


def rdf(mic_separations: np.ndarray, max_distance: float, res: int) -> np.array:
    """
    Calculates a discretized version of a function proportional to the RDF as a histogram

    :param mic_separations: [N,N,3] matrix of particle separations
    :param max_distance: maximum distance to plot for the histogram
    :param res: number of columns for the histogram
    :return: histograms of the separations weighted as 1/r^2
    """
    N = len(mic_separations)
    # turn vector separation into scalar distances
    mic_distances = np.linalg.norm(mic_separations, axis=2)
    distances = mic_distances.reshape(N*N, 1)  # reshape into a list
    distances = distances[distances != 0]  # remove all the zeros
    weights = 2*distances**-2  # 2* since 'distances' only has i < j but we just want i!=j
    histogram = np.histogram(distances, bins=res, range=(
        0, max_distance), weights=weights)
    return histogram[0]


def rdf_normalize(hist: np.array, t_count: float, dr: float, N: int, rho: float) -> np.array:
    """
    Normalizes the cumulative histogram produced by obs.rdf()
    Implementing the remaining coefficients and the time average
    for the histogram to be the function g(r) 

    :param hist: the cumulative histogram 
    :param t_count: number of times that histograms were added to make the cumulative histogram
    :param dr: the column width of the histogram
    :param N: the number of particles
    :param rho: the particle density
    :return: the correct normalized histogram representing a discrete RDF for the simulation
    """
    C = (N*4*np.pi*rho*dr)**-1
    return C*hist/t_count
