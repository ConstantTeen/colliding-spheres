import random
import numpy as np
import matplotlib.pyplot as plt

from utils import create_spheres, draw_sphere, plot
from classes import Sphere, Velocity


if __name__ == '__main__':
    spheres_num = 10
    g = 9.8
    radius_min, radius_max = 1, 5
    x_min, x_max = 0, 20
    y_min, y_max = 0, 20
    z_min, z_max = 0, 20
    log_param = 100
    t_max = 10.0

    random_state = 227
    random.seed(random_state)

    spheres = create_spheres(n=spheres_num, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, z_min=z_min, z_max=z_max, radius_min=radius_min, radius_max=radius_max)

    plt.rcParams['figure.figsize'] = (5, 5)
    plot(spheres, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, z_min=z_min, z_max=z_max)

    tau = t_max / 2e1
    current_t = 0
    current_step = 0

    while current_t < t_max:
        dv_arr = [Velocity() for _ in range(spheres_num)]

        for i in range(spheres_num):
            for j in range(i+1, spheres_num):
                if spheres[i].collides_with(spheres[j]):
                    ui, uj = spheres[i].velocities_after_collision(spheres[j])
                    dv_arr[i] = dv_arr[i] + ui
                    dv_arr[j] = dv_arr[j] + uj

        for i, dv in enumerate(dv_arr):
            spheres[i].push_to(tau=tau, dv=dv)

        plot(spheres, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, z_min=z_min, z_max=z_max)

        if current_step % log_param == 0:
            # TODO: запилить логирование
            pass

        current_t += tau
        current_step += 1


# TODO:
#  - discover how to calculate velocities after collision
#  - discover how to locate cube borders and deploy their interaction with spheres
#  - deploy a mechanism of logging
