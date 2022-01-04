from classes import System


def test_sphere():
    params = {
        'spheres_num': 1,
        'g': 9.8,
        'radius_min': 1,
        'radius_max': 5,
        'x_min': 0,
        'x_max': 20,
        'y_min': 0,
        'y_max': 20,
        'z_min': 0,
        'z_max': 20,
        'log_param': 100,
        't_max': 10,
        'random_state': 228,
        'mass_min': 0.1,
        'mass_max': 10.1,
        'v0_min': -10.,
        'v0_max': 10.,
        'iter_limit': 100
    }

    for k in range(1, 50):
        params['spheres_num'] = k

        system = System(**params)

        assert system.params['g'] == params['g']
        assert system.params['x_min'] == params['x_min']
        assert system.params['x_max'] == params['x_max']
        assert system.params['y_min'] == params['y_min']
        assert system.params['y_max'] == params['y_max']
        assert system.params['z_min'] == params['z_min']
        assert system.params['z_max'] == params['z_max']
        assert system.params['log_param'] == params['log_param']
        assert system.params['t_max'] == params['t_max']
        assert system.params['random_state'] == params['random_state']
        assert system.params['mass_min'] == params['mass_min']
        assert system.params['mass_max'] == params['mass_max']
        assert system.params['v0_min'] == params['v0_min']
        assert system.params['v0_max'] == params['v0_max']

        spheres = system.spheres

        assert len(spheres) == params['spheres_num']

        for s in spheres:
            assert s.radius >= params['radius_min']
            assert s.radius <= params['radius_max']
            assert s.coords[0] >= params['x_min'] + s.radius
            assert s.coords[0] <= params['x_max'] - s.radius
            assert s.coords[1] >= params['y_min'] + s.radius
            assert s.coords[1] <= params['y_max'] - s.radius
            assert s.coords[2] >= params['z_min'] + s.radius
            assert s.coords[2] <= params['z_max'] - s.radius
            assert s.v0[0] >= params['v0_min']
            assert s.v0[0] <= params['v0_max']
            assert s.v0[1] >= params['v0_min']
            assert s.v0[1] <= params['v0_max']
            assert s.v0[2] >= params['v0_min']
            assert s.v0[2] <= params['v0_max']
            assert s.mass >= params['mass_min']
            assert s.mass <= params['mass_max']
            assert s.time_after_last_hit == 0

    params['spheres_num'] = 1
    params['v0_min'] = 0
    params['v0_max'] = 0

    system = System(**params)
    s = system.spheres[0]

    assert system.get_current_v(s) == [0, 0, 0]
    assert s.time_after_last_hit == 0

    tau = params['t_max'] / 100

    x0, y0, z0 = s.coords
    system.make_a_move(s, tau)

    assert s.time_after_last_hit == tau
    assert s.coords[0] == x0
    assert s.coords[1] == y0
    assert abs(s.coords[2] - z0 + params['g'] * tau * tau / 2) < 1e-7

