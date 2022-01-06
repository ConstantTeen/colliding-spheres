import matplotlib.pyplot as plt
import sys
from utils import plot
from classes import SystemBuilder


if __name__ == '__main__':
    input_file_path = sys.argv[1]

    params = dict()

    with open(input_file_path, 'r') as f:
        params['spheres_num'] = int(f.readline())
        params['g'] = float(f.readline())
        params['radius_min'] = float(f.readline())
        params['radius_max'] = float(f.readline())
        params['mass_min'] = float(f.readline())
        params['mass_max'] = float(f.readline())

        x_line = [f.readline().rstrip(), f.readline(), f.readline()]
        if x_line[0] not in ('|--', '--|', '|-|', '---'):
            raise Exception(f'Only "|--", "--|", "---" and "|-|" borders types are allowed for x-axes, not {x_line[0]}')
        params['x_border_type'], params['x_min'], params['x_max'] = x_line[0], float(x_line[1]), float(x_line[2])

        y_line = [f.readline().rstrip(), f.readline(), f.readline()]
        if y_line[0] not in ('|--', '--|', '|-|', '---'):
            raise Exception(f'Only "|--", "--|", "---" and "|-|" borders types are allowed for y-axes, not {y_line[0]}')
        params['y_border_type'], params['y_min'], params['y_max'] = y_line[0], float(y_line[1]), float(y_line[2])

        z_line = [f.readline().rstrip(), f.readline(), f.readline()]
        if z_line[0] not in ('|--', '--|', '|-|', '---'):
            raise Exception(f'Only "|--", "--|", "---" and "|-|" borders types are allowed for z-axes, not {z_line[0]}')
        params['z_border_type'], params['z_min'], params['z_max'] = z_line[0], float(z_line[1]), float(z_line[2])

        params['log_param'] = int(f.readline())
        params['t_max'] = float(f.readline())
        params['random_state'] = 228

    system = SystemBuilder(**params)
    system.start_loop()
    # plot(system.spheres)

