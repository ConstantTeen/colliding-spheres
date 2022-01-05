import random
from utils import plot
from os.path import join


class Velocity:
    def __init__(self, arr=(0., 0., 0.)):
        assert len(arr) == 3, 'Wrong dimensionality'
        self.components = list(arr)

    def scalar_product(self, other):
        other = Velocity(other)
        return sum(self.components[i] * other.components[i] for i in range(3))

    def __repr__(self):
        return f'velocity[{self[0]}, {self[1]}, {self[2]}]'

    def __iter__(self):
        yield self[0]
        yield self[1]
        yield self[2]

    def __len__(self):
        return 3

    def __getitem__(self, item):
        assert item in (0, 1, 2), f'Wrong index {item}'
        return self.components[item]

    def __eq__(self, other):
        for i in range(3):
            if self[i] != other[i]:
                return False
        return True

    def __neg__(self):
        return (-1) * self

    def __add__(self, other):
        other = Velocity(other)
        return Velocity((self[0] + other[0], self[1] + other[1], self[2] + other[2]))

    def __mul__(self, const):
        return Velocity((self[0]*const, self[1]*const, self[2]*const))

    def __rmul__(self, const):
        return self*const

    def __sub__(self, other):
        other = Velocity(other)
        return self + other * (-1)

    def __truediv__(self, const):
        assert const != 0, 'Division by zero'
        return self * (1/const)

    def __abs__(self):
        return (self[0]**2 + self[1]**2 + self[2]**2)**0.5


class System:
    def __init__(self, **kwargs):
        self.spheres = []
        self.params = kwargs
        random.seed(kwargs['random_state'])
        self.__create_spheres()

    def __create_spheres(self):
        spheres_num = self.params['spheres_num']
        iter_limit = self.params['iter_limit']
        radius_min = self.params['radius_min']
        radius_max = self.params['radius_max']
        mass_min = self.params['mass_min']
        mass_max = self.params['mass_max']
        x_min = self.params['x_min']
        y_min = self.params['y_min']
        z_min = self.params['z_min']
        x_max = self.params['x_max']
        y_max = self.params['y_max']
        z_max = self.params['z_max']
        v0_min = self.params['v0_min']
        v0_max = self.params['v0_max']

        ind = 0
        iter_ind = 0

        while ind < spheres_num and iter_ind < iter_limit:  # place spheres
            params = {}
            radius_ = radius_min + (radius_max - radius_min) * random.random()
            params['radius'] = radius_
            params['mass'] = mass_min + (mass_max - mass_min) * random.random()

            params['center_coords'] = [None, None, None]

            for i, min_, max_ in zip([0, 1, 2], [x_min, y_min, z_min], [x_max, y_max, z_max]):
                params['center_coords'][i] = min_ + radius_ + (max_ - min_ - 2 * radius_) * random.random()

            params['v'] = Velocity([
                v0_min + (v0_max - v0_min) * random.random(),
                v0_min + (v0_max - v0_min) * random.random(),
                v0_min + (v0_max - v0_min) * random.random()
            ])

            s = Sphere(**params)

            flag = False

            for other in self.spheres:
                if self.collides_with(s, other):
                    flag = True
                    break

            if flag:
                continue

            self.spheres += [s]
            ind += 1
            iter_ind += 1

        if iter_ind >= iter_limit:
            raise Exception('Spheres placement is failed! Location search timed out.')

    def get_spheres(self):
        return self.spheres

    def collides_with(self, sphere1, sphere2):
        assert isinstance(sphere1, Sphere)
        assert isinstance(sphere2, Sphere)

        return self.distance(sphere1, sphere2) <= sphere1.radius + sphere2.radius

    def get_colliding_spheres(self):
        spheres_num = self.params['spheres_num']

        for i in range(spheres_num):
            for j in range(i+1, spheres_num):
                if self.collides_with(self.spheres[i], self.spheres[j]):
                    yield i, j

    def get_out_of_bounds(self):
        for i, s in enumerate(self.spheres):
            bounds = []

            if self.params['x_border_type'][0] == '|' and s.coords[0] <= self.params['x_min'] + s.radius:
                bounds += ['yz']
            if self.params['x_border_type'][2] == '|' and s.coords[0] >= self.params['x_max'] - s.radius:
                bounds += ['YZ']
            if self.params['y_border_type'][0] == '|' and s.coords[1] <= self.params['y_min'] + s.radius:
                bounds += ['xz']
            if self.params['y_border_type'][2] == '|' and s.coords[1] >= self.params['y_max'] - s.radius:
                bounds += ['XZ']
            if self.params['z_border_type'][0] == '|' and s.coords[2] <= self.params['z_min'] + s.radius:
                bounds += ['xy']
            if self.params['z_border_type'][2] == '|' and s.coords[2] >= self.params['z_max'] - s.radius:
                bounds += ['XY']

            assert len(bounds) <= 3

            yield i, bounds

    def start_loop(self, observations_num=10000, outdir='outdir'):
        t_max = self.params['t_max']

        tau = t_max / observations_num
        current_t = 0
        current_step = 0

        while current_t < t_max:
            current_step += 1
            current_t += tau

            for s in self.spheres:
                s.coords[0] += tau * s.v[0]
                s.coords[1] += tau * s.v[1]
                s.coords[2] += tau * s.v[2]

            for i, j in self.get_colliding_spheres():
                si, sj = self.spheres[i], self.spheres[j]
                mi, mj = si.mass, sj.mass

                r_rel = [si.coords[0] - sj.coords[0], si.coords[1] - sj.coords[1], si.coords[2] - sj.coords[2]]
                v_rel = si.v - sj.v
                v_cm = (mi * si.v + mj * sj.v) / (mi + mj)

                v_rel = 2 * Velocity(r_rel) * Velocity(r_rel).scalar_product(v_rel) / abs(Velocity(r_rel))**2 - v_rel

                si.v = v_cm + v_rel * mj / (mi + mj)
                sj.v = v_cm - v_rel * mi / (mi + mj)

            for i, bounds in [(i, bounds) for i, bounds in self.get_out_of_bounds() if len(bounds) > 0]:
                s = self.spheres[i]

                if 'yz' in bounds:
                    s.coords[0] = self.params['x_min'] + s.radius
                    s.v = Velocity([-s.v[0], s.v[1], s.v[2]])
                if 'YZ' in bounds:
                    s.coords[0] = self.params['x_max'] - s.radius
                    s.v = Velocity([-s.v[0], s.v[1], s.v[2]])
                if 'xz' in bounds:
                    s.coords[1] = self.params['y_min'] + s.radius
                    s.v = Velocity([s.v[0], -s.v[1], s.v[2]])
                if 'XZ' in bounds:
                    s.coords[1] = self.params['y_max'] - s.radius
                    s.v = Velocity([s.v[0], -s.v[1], s.v[2]])
                if 'xy' in bounds:
                    s.coords[2] = self.params['z_min'] + s.radius
                    s.v = Velocity([s.v[0], s.v[1], -s.v[2]])
                if 'XY' in bounds:
                    s.coords[2] = self.params['z_max'] - s.radius
                    s.v = Velocity([s.v[0], s.v[1], -s.v[2]])

            for s in self.spheres:
                s.v = s.v - Velocity([0, 0, s.mass * self.params['g'] * tau])

            if current_step % self.params['log_param'] == 0:
                plot(self.spheres, **{
                    'x_min': self.params['x_min'],
                    'x_max': self.params['x_max'],
                    'y_min': self.params['y_min'],
                    'y_max': self.params['y_max'],
                    'z_min': self.params['z_min'],
                    'z_max': self.params['z_max'],
                })

                # print(f'{current_step=}; {current_t=:.4f}')
                #
                # with open(join(outdir, f'{current_step // self.params["log_param"]:#08}.plt'), 'w') as f:
                #     f.write(
                #         f'TITLE "Title, {current_t=:.4f}; x_min={self.params["x_min"]:.4f}; x_max={self.params["x_max"]:.4f}; y_min={self.params["y_min"]:.4f}; y_max={self.params["y_max"]:.4f}; z_min={self.params["z_min"]:.4f}; z_max={"+inf" if None is self.params["z_max"] else self.params["z_max"]:.4f}; spheres_num={len(self.spheres)}"\n'
                #     )
                #     f.write(f'VARIABLES = "X" "Y" "Z" "R" "M"\n')
                #
                #     for s in self.spheres:
                #         f.write(f'{s.coords[0]:.4f} {s.coords[1]:.4f} {s.coords[2]:.4f} {s.radius:.4f} {s.mass:.4f}\n')

    @staticmethod
    def distance(sphere1, sphere2):
        x1, y1, z1 = sphere1.coords[0], sphere1.coords[1], sphere1.coords[2]
        x2, y2, z2 = sphere2.coords[0], sphere2.coords[1], sphere2.coords[2]
        return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5


@lambda c: c()
class SystemBuilder:
    system = None

    def __call__(self, random_state=228, v0_min=-10., v0_max=10., iter_limit=100, *args, **kwargs):
        if self.system is not None:
            return self.system

        kwargs['random_state'] = random_state
        kwargs['v0_min'] = v0_min
        kwargs['v0_max'] = v0_max
        kwargs['iter_limit'] = iter_limit

        self.system = System(**kwargs)
        return self.system


class Sphere:
    def __init__(self, center_coords, radius, mass, v=(0, 0, 0)):
        self.coords = center_coords  # array like
        self.radius = radius
        self.mass = mass
        self.v = Velocity(v)
        self.time_after_last_hit = 0

    def __repr__(self):
        return f"sphere<({self.coords[0]:.2}, {self.coords[1]:.2}, {self.coords[2]:.2}); r={self.radius:.2}>\n"

    def get_props(self):
        return {
            'x': self.coords[0],
            'y': self.coords[1],
            'z': self.coords[2],
            'radius': self.radius,
            'mass': self.mass,
            'v0': self.v,
            't': self.time_after_last_hit
        }

