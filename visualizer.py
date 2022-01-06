import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from os import listdir, unlink
from os.path import join
from utils import *
import sys


def read_data(folder):
    data = []

    for file_name in [fn for fn in listdir(folder) if fn.endswith('.plt')]:
        with open(join(folder, file_name)) as f:
            current_t, x_min, x_max, y_min, y_max, z_min, z_max, spheres_num = map(lambda x: x.split('=')[1],
                                                                                   f.readline()[:-2].split(';'))
            current_t, x_min, x_max, y_min, y_max, z_min, z_max = map(float,
                                                                      [current_t, x_min, x_max, y_min, y_max, z_min,
                                                                       z_max])
            spheres_num = int(spheres_num)
            f.readline()  # VARIABLES = "X" "Y" "Z" "R" "M"
            data += [np.genfromtxt(f)]

    return {
        'data': data,
        'x_min': x_min,
        'x_max': x_max,
        'y_min': y_min,
        'y_max': y_max,
        'z_min': z_min,
        'z_max': z_max,
        'spheres_num': spheres_num
    }


def update_lines(num, data, lines):
    for line, sphere in zip(lines, data):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(sphere[num, :2].T)
        line.set_3d_properties(sphere[num, 2])

    return lines


if __name__ == '__main__':
    folder = sys.argv[1]
    print(f'Retrieving data from {folder}...')
    all_data_dict = read_data(folder)

    data = all_data_dict['data']
    spheres_num = all_data_dict['spheres_num']

    data = [np.array([data[j][i, :] for j in range(len(data))]) for i in
            range(spheres_num)]  # shape = (spheres_num, observation_num, 5), type = list of np.arrays
    print(f'Animation construction...')
    # animation
    fig = plt.figure()
    # fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax = fig.add_subplot(111, projection='3d', autoscale_on=False)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(all_data_dict['x_min'], all_data_dict['x_max'])
    ax.set_ylim(all_data_dict['y_min'], all_data_dict['y_max'])
    ax.set_zlim(all_data_dict['z_min'], all_data_dict['z_max'])

    colors = [mass2color(data[i][0][-1]) for i in range(spheres_num)]
    radiuses = [data[i][0][-2] for i in range(spheres_num)]
    num_steps = len(data[0])

    lines = [ax.plot([], [], [], 'o', color=colors[i], markersize=1e1 * radiuses[i], alpha=0.5)[0] for i in
             range(spheres_num)]

    anime = animation.FuncAnimation(fig, update_lines, num_steps, fargs=(data, lines), interval=10)
    
    gif_path = 'my_balls.gif'
    print(f'Saving gif to file {gif_path}...')
    anime.save(gif_path, fps=30)
    print(f'Done!')