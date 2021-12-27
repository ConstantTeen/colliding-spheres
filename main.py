import random
import numpy as np
import matplotlib.pyplot as plt

from utils import  draw_sphere, plot
from classes import SystemBuilder, Sphere, Velocity


if __name__ == '__main__':
    spheres_num = 10
    g = 9.8
    radius_min, radius_max = 1, 5
    x_min, x_max = 0, 20
    y_min, y_max = 0, 20
    z_min, z_max = 0, None
    log_param = 100
    t_max = 10.0

    random_state = 227
    random.seed(random_state)

    params = {
        'spheres_num': spheres_num,
        'g': g,
        'radius_min': radius_min,
        'radius_max': radius_max,
        'x_min': x_min,
        'x_max': x_max,
        'y_min': y_min,
        'y_max': y_max,
        'z_min': z_min,
        'z_max': z_max,
        'log_param': log_param,
        't_max': t_max,
        'random_state': random_state
    }

    system = SystemBuilder(**params)

    spheres = system.get_spheres()

    plt.rcParams['figure.figsize'] = (5, 5)
    plot(spheres, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, z_min=z_min, z_max=z_max)
    system.start_loop(1000)


# TODO:
#  - discover how to calculate velocities after collision
#  - discover how to locate cube borders and deploy their interaction with spheres
#  - deploy a mechanism of logging
