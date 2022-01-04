from classes import Velocity
import random


def test_velocity():
    for _ in range(100):
        x, y, z = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)
        v = Velocity((x, y, z))

        assert len(v) == 3
        assert (v[0], v[1], v[2]) == (x, y, z)
        assert v == (x, y, z)
        assert v != (x+1, y, z)
        assert 2*v == (2*x, 2*y, 2*z)
        assert v*2 == (2*x, 2*y, 2*z)
        assert v*0 == (0, 0, 0)
        assert abs(v/123 - (x/123, y/123, z/123)) < 1e-5
        assert abs(v / (-123) - (x / (-123), y / (-123), z / (-123))) < 1e-5

        for i, el in enumerate(v):
            assert el == [x, y, z][i]

        assert abs(v) == (x**2 + y**2 + z**2)**0.5

    for _ in range(100):
        x1, y1, z1 = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)
        x2, y2, z2 = random.randint(-100, 100), random.randint(-100, 100), random.randint(-100, 100)
        v1 = Velocity((x1, y1, z1))
        v2 = Velocity((x2, y2, z2))

        assert v1.scalar_product(v2) == x1*x2 + y1*y2 + z1*z2
        assert v2.scalar_product(v1) == x1*x2 + y1*y2 + z1*z2
        assert v1 + v2 == (x1+x2, y1+y2, z1+z2)
        assert v1 - v2 == (x1-x2, y1-y2, z1-z2)
