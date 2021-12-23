import numpy as np


class Velocity:
    def __init__(self, arr=(0., 0., 0.)):
        assert len(arr) == 3, 'Wrong dimensionality'
        self.components = list(arr)

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


class Sphere:
    def __init__(self, center_coords, radius, mass, g, v0=(0, 0, 0)):
        self.coords = center_coords  # array like
        self.radius = radius
        self.mass = mass
        self.v = Velocity(v0)
        self.time_after_last_hit = 0
        self.g = g

    def push_to(self, tau, dv=(0, 0, 0)):
        if dv != (0, 0, 0):
            self.v = self.v + Velocity(dv)
            self.time_after_last_hit = 0

        self.coords[0] += self.v[0] * tau
        self.coords[1] += self.v[1] * tau
        self.coords[2] += self.g / 2 * (2 * self.time_after_last_hit * tau + tau ** 2)
        self.time_after_last_hit += tau

    def velocities_after_collision(self, other):
        assert isinstance(other, Sphere)
        # TODO: допилить физику collision physics https://ru.wikipedia.org/wiki/%D0%A3%D0%B4%D0%B0%D1%80
        return -self.v, -other.v

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
        x1, y1, z1 = self.coords[0], self.coords[1], self.coords[2]
        x2, y2, z2 = other.get_props()['x'], other.get_props()['y'], other.get_props()['z']
        return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5

    def collides_with(self, other):
        return self.distance(other) <= (self.radius + other.radius)
