from utils import mass2color


def test_mass2color():
    assert mass2color(0) == (0, 0, 1)
    assert mass2color(0, light_color=(0.7, 0.1, 0.4)) == (0.7, 0.1, 0.4)
