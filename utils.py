import random
import numpy as np
import matplotlib.pyplot as plt

from classes import Sphere, Velocity


def create_spheres(n,
                   radius_min, radius_max,
                   x_min, x_max, y_min, y_max, z_min,  z_max=100.,
                   mass_min=0.1, mass_max=10.1,
                   v0_min=-10., v0_max=10.,
                   iter_limit=100,
                   g=9.8):
    i = 0
    j = 0
    spheres = []

    while i < n and j < iter_limit:  # place spheres
        params = {}
        radius_ = radius_min + (radius_max - radius_min) * random.random()
        params['radius'] = radius_
        params['mass'] = mass_min + (mass_max - mass_min) * random.random()
        params['g'] = g

        # z_max is None => z_max == +inf, but it's a wrong place to think about it
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
            if s.collides_with(other):
                continue

        spheres += [s]
        i += 1
        j += 1

    assert j < iter_limit, 'Spheres placement is failed! Location search timed out.'

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
    # thanks to https://stackoverflow.com/questions/32424670/python-matplotlib-drawing-3d-sphere-with-circumferences
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


def plot(arr, x_min=0, x_max=20, y_min=0, y_max=20, z_min=0, z_max=20):
    ax = plt.axes(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)

    for s in arr:
        draw_sphere(s, ax)

    plt.show()