import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Velocity:
    def __init__(self, arr=(0, 0, 0)):
        assert len(arr) == 3, 'Wrong dimensionality'
        self.components = arr

    def __repr__(self):
        return f'velocity[{self[0]}, {self[1]}, {self[2]}]'

    def __getitem__(self, item):
        assert item in (0,1,2), 'Wrong index'
        return self.components[item]

    def __eq__(self, other):
        for i in range(3):
            if self[i] != other[i]:
                return False
        return True

    def __add__(self, other):
        return Velocity((self[0] + other[0], self[1] + other[1], self[2] + other[2]))

    def __mul__(self, const):
        return Velocity((self[0]*const, self[1]*const, self[2]*const))

    def __rmul__(self, const):
        return self*const

    def __sub__(self, other):
        return self + other * (-1)

    def __truediv__(self, const):
        assert const != 0, 'Division by zero'
        return self * (1/const)

    def __abs__(self):
        return (self[0]**2 + self[1]**2 + self[2]**2)**0.5


class Sphere:
    def __init__(self, center_coords, radius, mass, v0=(0, 0, 0)):
        self.coords = center_coords  # array like
        self.radius = radius
        self.mass = mass
        self.v = Velocity(v0)

    def __repr__(self):
        return f"sphere({self.coords}-{self.radius})"

    def get_props(self):
        return {
            'x': self.coords[0],
            'y': self.coords[1],
            'z': self.coords[2],
            'radius': self.radius,
            'mass': self.mass,
            'v': self.v
        }

    def distance(self, other):
        x1, y1 = self.coords[0], self.coords[1]
        x2, y2 = other.get_props()['x'], other.get_props()['y']
        return ((x1 - x2)**2 + (y1 - y2)**2)**0.5

    def are_colliding(self, other):
        return self.distance(other) <= (self.radius + other.radius)


def create_spheres(n, mass_min=0.1, mass_max=10.1, z_max=100., v0_min=-10., v0_max=10.):
    i = 0
    spheres = []

    while i < n:  # place spheres
        params = {}
        radius_ = radius_min + (radius_max - radius_min) * random.random()
        params['radius'] = radius_
        params['mass'] = mass_min + (mass_max - mass_min) * random.random()

        # z_max is None => z_max == +inf, but it's wrong place where to think about it
        params['center_coords'] = [
            x_min + radius_ + (x_max - x_min - 2 * radius_) * random.random(),
            y_min + radius_ + (y_max - y_min - 2 * radius_) * random.random(),
            z_min + radius_ + (z_max - z_min - 2 * radius_) * random.random()
        ]

        params['v0'] = [
            v0_min + (v0_max - v0_min) * random.random(),
            v0_min + (v0_max - v0_min) * random.random(),
            v0_min + (v0_max - v0_min) * random.random()
        ]

        s = Sphere(**params)

        for other in spheres:
            if s.are_colliding(other):
                continue

        spheres += [s]
        i += 1

    return spheres


def mass2color(mass):
    """
    Returns color in RGBA format (without alpha). The less the mass of sphere the closer its color to blue.
    The bigger the mass the closer sphere color to red.
    :param mass: mass of the sphere
    :return: tuple(3): rgba color
    """
    def sigmoid(x):  # returns number from 0 to 1 if x from 0 to +inf
        return 2 / (1 + np.exp(-x)) - 1

    mass = sigmoid(mass)
    return mass, 0, 1 - mass


def draw_sphere(s, ax):
    x_center, y_center, z_center = s.get_props()['x'], s.get_props()['y'], s.get_props()['z']
    r = s.get_props()['radius']
    m = s.get_props()['mass']
    # velocity = s.get_props()['v']

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = x_center + r * np.outer(np.cos(u), np.sin(v))
    y = y_center + r * np.outer(np.sin(u), np.sin(v))
    z = z_center + r * np.outer(np.ones(np.size(u)), np.cos(v))

    color = mass2color(m)
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color=color, linewidth=0, alpha=0.5)


if __name__ == '__main__':
    spheres_num = 10
    g = 9.8
    radius_min, radius_max = (1, 5)
    x_min, x_max = 0, 20
    y_min, y_max = 0, 20
    z_min, z_max = 0, 20
    log_param = 100
    t_max = 10.0
    random_state = 227

    random.seed(random_state)

    tau = t_max / 1e3
    current_t = 0
    current_step = 0

    spheres = create_spheres(spheres_num)

    # thanks to https://stackoverflow.com/questions/32424670/python-matplotlib-drawing-3d-sphere-with-circumferences
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('auto')

    for s in spheres:
        draw_sphere(s, ax)

    plt.show()



