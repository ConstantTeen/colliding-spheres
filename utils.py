import numpy as np
import matplotlib.pyplot as plt


def mass2color(mass, light_color=(0, 0, 1), heavy_color=(1, 0, 0)):
    """
    Returns color in RGBA format (without alpha). The less the mass of sphere the closer its color to blue.
    The bigger the mass the closer sphere color to red.
    :param mass: mass of the sphere
    :param light_color: color of object of mass = 0
    :param heavy_color: color of object of mass = +inf
    :return: tuple(3): rgba color
    """
    def sigmoid(x):  # returns number from 0 to 1 if x from 0 to +inf
        return 2 / (1 + np.exp(-x)) - 1

    mass = sigmoid(mass)
    return (1 - mass) * light_color[0] + mass * heavy_color[0], \
           (1 - mass) * light_color[1] + mass * heavy_color[1], \
           (1 - mass) * light_color[2] + mass * heavy_color[2]


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

    color = mass2color(m, light_color=(237/255, 255/255, 198/255), heavy_color=(0/255, 0/255, 0/255))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color=color, linewidth=0, alpha=0.5)


def plot(arr, x_min=0, x_max=20, y_min=0, y_max=20, z_min=0, z_max=20):
    if z_max is None:
        z_max = z_min + max((x_max-x_min), (y_max-y_min))

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
