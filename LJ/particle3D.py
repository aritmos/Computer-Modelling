"""
 CompMod Ex2: Particle3D, a class to describe point particles in 3D space
 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

author: Sebastian Garcia (s1910157)

"""

import numpy as np
import pbc 

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label - name of the particle
    mass  - mass of the particle
    pos   - position of the particle
    vel   - velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e      - computes the kinetic energy
    momentum       - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel     - updates the velocity

        Static Methods:
    new_particle - initializes P3D instances from a file handle
    sys_kinetic  - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """


    def __init__(self, label: str, mass: float, pos: list, vel: list) -> None:
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = label
        self.mass  = mass
        self.pos   = np.array(pos)
        self.vel   = np.array(vel)


    @staticmethod
    def new_particles(N:int,mass:float) -> list:
        p3d_list = []
        for i in range(0,N):
            p3d_list.append(Particle3D(f'p{i}',mass,[0,0,0],[0,0,0]))
        return p3d_list
        

    def __str__(self) -> str:
        """
        XYZ-compliant string. The format is
        <label>  <x>  <y>  <z>
        """
        return f'{self.label} {" ".join([str(i) for i in self.pos])}'


    def kinetic_e(self) -> float:
        """
        Returns the kinetic energy of a Particle3D instance

        :return: kinetic energy (float)
        """
        return 0.5*self.mass*np.sum(self.vel**2)


    def momentum(self) -> np.ndarray:
        """
        Returns the linear momentum of a Particle3D instance
        :return p: [3] float numpy array, m*v
        """
        return self.mass*self.vel


    def update_pos_1st(self, dt:float) -> None:
        """
        1st order position update

        :param dt: timestep
        """
        self.pos = self.pos + dt*self.vel


    def update_pos_2nd(self, l:float, dt:float, force:np.array) -> None:
        """
        2nd order position update using PBC

        :param l: box length
        :param dt: timestep
        :param force: [3] float numpy array, the total force acting on the particle
        """
        self.pos += dt*self.vel + (dt**2)*force/(2*self.mass)
        # Apply PBC
        self.pos = np.array(pbc.periodic_image(self.pos, l))

    def update_vel(self, dt: float, force: np.ndarray) -> None:
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float numpy array, the total force acting on the particle
        """
        self.vel = self.vel + dt*force/(self.mass)


    @staticmethod
    def sys_kinetic(p3d_list: 'list[Particle3D]') -> float:
        """
        Returns the kinetic energy of the whole given system

        :param p3d_list: list of P3D instances
        :return: total kinetic energy of the system
        """
        return sum([particle.kinetic_e() for particle in p3d_list])


    @staticmethod
    def com_velocity(p3d_list: 'list[Particle3D]') -> tuple:
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list of P3D instances
        :return total_mass: The total mass of the system 
        :return: (mass,Centre-of-mass velocity) touple 
        """
        M = sum([particle.mass for particle in p3d_list])
        return M, sum([particle.momentum() for particle in p3d_list])/M

    @staticmethod
    def pair_separations(p3d_list:'list[Particle3D]') -> np.ndarray:
        """
        Computes an [N,N,3] array of particle separations where

        separations[i][j] = r_{ij}

        :param p3d_list: list of P3D instances
        :return separations: [N,N,3] separation array 
        """
        N = len(p3d_list)
        separations = np.zeros([N,N,3])
        for i in range(N):
            for j in range(i+1,N):
                separations[i][j] = p3d_list[j].pos-p3d_list[i].pos
        return separations
