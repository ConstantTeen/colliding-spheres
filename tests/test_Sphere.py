from classes import Velocity, Sphere
import random


def test_sphere():
    for _ in range(100):
        x, y, z = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)
        radius = random.randint(0, 100)
        mass = random.randint(0, 100)
        vx, vy, vz = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)

        s = Sphere(center_coords=(x, y, z), radius=radius, mass=mass, v0=(vx, vy, vz))
        props = s.get_props()

        assert props['x'] == x
        assert props['y'] == y
        assert props['z'] == z
        assert props['radius'] == radius
        assert props['mass'] == mass
        assert props['v0'] == (vx, vy, vz)
        assert props['t'] == 0
