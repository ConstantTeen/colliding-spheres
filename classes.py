import numpy as np
import random
from utils import plot


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


class System:
    params = dict()
    spheres = []

    def __init__(self, **kwargs):
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

        i = 0
        j = 0

        while i < spheres_num and j < iter_limit:  # place spheres
            params = {}
            radius_ = radius_min + (radius_max - radius_min) * random.random()
            params['radius'] = radius_
            params['mass'] = mass_min + (mass_max - mass_min) * random.random()

            if z_max is None:
                params['center_coords'] = [
                    x_min + radius_ + (x_max - x_min - 2 * radius_) * random.random(),
                    y_min + radius_ + (y_max - y_min - 2 * radius_) * random.random(),
                    z_min + radius_ + max((x_max - x_min), (y_max - y_min)) * random.random()
                ]
            else:
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

            for other in self.spheres:
                if self.collides_with(s, other):
                    continue

            self.spheres += [s]
            i += 1
            j += 1

        assert j < iter_limit, 'Spheres placement is failed! Location search timed out.'

    def get_spheres(self):
        return self.spheres

    def __get_colliding_spheres(self):
        spheres_num = self.params['spheres_num']

        for i in range(spheres_num):
            for j in range(i+1, spheres_num):
                if self.collides_with(self.spheres[i], self.spheres[j]):
                    yield i, j

    def start_loop(self, observations_num=500):
        t_max = self.params['t_max']
        spheres_num = self.params['spheres_num']

        tau = t_max / observations_num
        current_t = 0
        current_step = 0

        while current_t < t_max:
            # Сферы делятся на несколько типов
            # 1. ничего и никого не касаются. Такие сферы сохраняют траекторию движения, а значит можно просто += v0z*tau + g/2*tau*(2t+tau), v0xy*tau; t+=tau
            # 2. сфера касается другой сферы или нескольких. Нужно посчитать результирующий вектор скорости после столкновения. Тогда эти tau времени сфера полетит по новой траектории
            #     при этом этот результирующий вектор будет вектором начальной скорости. v0=v; r += v0*tau; t = tau
            # 3. Сфера касается стенки. Считаем вектор скорости после столкновения и v0=v; r += v0*tau; t = tau
            # 4. Сфера касается и стенки и других сфер. тоже самое только результ вектор складывается по стенкам и сферам
            dv_arr = [Velocity() for _ in range(spheres_num)]
            borders = ('xy', 'XY', 'xz', 'XZ', 'yz', 'YZ')  # TODO: вынести в константное поле класса

            for i, s in enumerate(self.spheres):
                for border in borders:
                    if self.collides_with(s, border):
                        dv_arr[i], _ = self.velocities_after_collision(s, border)

            colliding_pairs = self.__get_colliding_spheres()

            for i, j in colliding_pairs:
                ui, uj = self.velocities_after_collision(self.spheres[i], self.spheres[j])
                dv_arr[i] = dv_arr[i] + ui
                dv_arr[j] = dv_arr[j] + uj

            for i, s in enumerate(self.spheres):
                self.make_a_move(s, tau, dv_arr[i])

            plot(self.spheres, **{
                'x_min': self.params['x_min'],
                'x_max': self.params['x_max'],
                'y_min': self.params['y_min'],
                'y_max': self.params['y_max'],
                'z_min': self.params['z_min'],
                'z_max': self.params['z_max'],
            })

            current_step += 1
            current_t += tau

    def make_a_move(self, sphere, tau, dv=(0, 0, 0)):
        if dv == (0, 0, 0):
            sphere.coords[0] += sphere.v0[0] * tau
            sphere.coords[1] += sphere.v0[1] * tau
            sphere.coords[2] += sphere.v0[2] * tau + self.params['g'] / 2 * tau * (2 * sphere.time_after_last_hit + tau)
            sphere.time_after_last_hit += tau
            return
        sphere.v0 = dv
        sphere.coords[0] += sphere.v0[0] * tau
        sphere.coords[1] += sphere.v0[1] * tau
        sphere.coords[2] += sphere.v0[2] * tau
        sphere.time_after_last_hit = tau

    def get_current_v(self, sphere):
        return sphere.v0 + Velocity((0, 0, self.params['g'] * sphere.time_after_last_hit))

    @staticmethod
    def distance(sphere1, sphere2):
        x1, y1, z1 = sphere1.coords[0], sphere1.coords[1], sphere1.coords[2]
        x2, y2, z2 = sphere2.coords[0], sphere2.coords[1], sphere2.coords[2]
        return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5

    def collides_with(self, sphere1, sphere2):
        assert \
            (isinstance(sphere1, Sphere) or isinstance(sphere1, str)) and \
            (isinstance(sphere2, Sphere) or isinstance(sphere2, str)), \
            f'Wrong type of input data! {type(sphere1)}, {type(sphere2)}'

        if isinstance(sphere1, Sphere) and isinstance(sphere2, Sphere):
            return self.distance(sphere1, sphere2) <= (sphere1.radius + sphere2.radius)
        if isinstance(sphere1, Sphere) and isinstance(sphere2, str):
            assert sphere2 in ('xy', 'XY', 'xz', 'XZ', 'yz', 'YZ'), f'{sphere2} is not in available axes names'
            if sphere2 == 'xy':
                return sphere1.coords[2] <= sphere1.radius + self.params['z_min']
            if sphere2 == 'XY':
                if self.params['z_max'] is None:
                    return False
                return sphere1.coords[2] >= self.params['z_max'] - sphere1.radius
            if sphere2 == 'xz':
                return sphere1.coords[1] <= sphere1.radius + self.params['y_min']
            if sphere2 == 'XZ':
                return sphere1.coords[1] <= self.params['y_max'] - sphere1.radius
            if sphere2 == 'yz':
                return sphere1.coords[0] <= sphere1.radius + self.params['x_min']
            if sphere2 == 'YZ':
                return sphere1.coords[0] <= self.params['x_max'] - sphere1.radius
        assert 1 != 1, f'{type(sphere1)} is not a sphere!'

    def velocities_after_collision(self, sphere1, sphere2):
        assert self.collides_with(sphere1, sphere2), 'Objects are not colliding!'

        if isinstance(sphere1, Sphere) and isinstance(sphere2, Sphere):
            # TODO: допилить физику collision physics https://ru.wikipedia.org/wiki/%D0%A3%D0%B4%D0%B0%D1%80
            return -self.get_current_v(sphere1), -self.get_current_v(sphere2)
        if isinstance(sphere1, Sphere) and isinstance(sphere2, str):
            assert sphere2 in ('xy', 'XY', 'xz', 'XZ', 'yz', 'YZ'), f'{sphere2} is not in available axes names'

            dict_n = {
                'xy': Velocity((0, 0, -1)),
                'XY': Velocity((0, 0, 1)),
                'xz': Velocity((0, -1, 0)),
                'XZ': Velocity((0, 1, 0)),
                'yz': Velocity((-1, 0, 0)),
                'YZ': Velocity((1, 0, 0))
            }
            v = self.get_current_v(sphere1)
            u = v - 2*v.scalar_product(dict_n[sphere2])*dict_n[sphere2]

            return u, None

        assert 1 != 1, 'Wrong object name. Spheres or strings allowed only'


@lambda c: c()
class SystemBuilder:
    system = None

    def __call__(
                self, random_state=228,
                z_max=None,
                mass_min=0.1, mass_max=10.1,
                v0_min=-10., v0_max=10.,
                iter_limit=100,
                g=9.8,
                *args, **kwargs):
        if self.system is not None:
            return self.system

        if 'random_state' not in kwargs:
            kwargs['random_state'] = random_state
        if 'z_max' not in kwargs:
            kwargs['z_max'] = z_max
        if 'mass_min' not in kwargs:
            kwargs['mass_min'] = mass_min
        if 'mass_max' not in kwargs:
            kwargs['mass_max'] = mass_max
        if 'v0_min' not in kwargs:
            kwargs['v0_min'] = v0_min
        if 'v0_max' not in kwargs:
            kwargs['v0_max'] = v0_max
        if 'iter_limit' not in kwargs:
            kwargs['iter_limit'] = iter_limit
        if 'g' not in kwargs:
            kwargs['g'] = g

        self.system = System(**kwargs)
        return self.system


class Sphere:
    def __init__(self, center_coords, radius, mass, v0=(0, 0, 0)):
        self.coords = center_coords  # array like
        self.radius = radius
        self.mass = mass
        self.v0 = Velocity(v0)
        self.time_after_last_hit = 0






    def __repr__(self):
        return f"sphere({self.coords}-{self.radius})"

    def get_props(self):
        return {
            'x': self.coords[0],
            'y': self.coords[1],
            'z': self.coords[2],
            'radius': self.radius,
            'mass': self.mass,
            'v': self.v0
        }





# TODO: перенести все методы в систему


