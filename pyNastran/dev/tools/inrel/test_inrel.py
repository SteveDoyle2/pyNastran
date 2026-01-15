import unittest
import numpy as np
from pyNastran.dev.tools.inrel.inertia_relief2 import (
    inertia_relief)


class TestInrel(unittest.TestCase):
    def test_inrel_bar_force_constant(self):
        nnode = 11
        mass_total = 11.
        xyz, mass = get_bar(
            mass_total=mass_total, nnode=nnode)
        gx = 0.1
        gy = 0.0
        gz = 0.5
        fx_expected = -gx * mass[0]
        fy_expected = -gy * mass[0]
        fz_expected = -gz * mass[0]

        fx = gx * mass
        fy = gy * mass
        fz = gz * mass
        force = np.column_stack([fx, fy, fz])
        moment = force * 0.0
        inertia = np.zeros((nnode, 6), dtype='float64')
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment
        # print('df:\n', dforce)
        # print('dm:\n', dmoment)
        assert np.allclose(dforce[0, 0], fx_expected), (dforce[0, 0], fx_expected)
        assert np.allclose(dforce[0, 1], fy_expected), (dforce[0, 1], fy_expected)
        assert np.allclose(dforce[0, 2], fz_expected), (dforce[0, 2], fz_expected)

        assert np.allclose(dforce[:, 0].max(), fx_expected), (dforce[:, 0].max(), fx_expected)
        assert np.allclose(dforce[:, 1].max(), fy_expected), (dforce[:, 1].max(), fy_expected)
        assert np.allclose(dforce[:, 2].max(), fz_expected), (dforce[:, 2].max(), fz_expected)

        assert np.allclose(dforce[:, 0].min(), fx_expected), (dforce[:, 0].min(), fx_expected)
        assert np.allclose(dforce[:, 1].min(), fy_expected), (dforce[:, 1].min(), fy_expected)
        assert np.allclose(dforce[:, 2].min(), fz_expected), (dforce[:, 2].min(), fz_expected)

        assert np.allclose(force_out.sum(), 0.)
        assert np.allclose(moment_out.sum(), 0.)
        # print('force_out:\n', force_out)

    def test_inrel_bar_moment_constant(self):
        nnode = 5
        mass_total = 0.02
        xyz, mass = get_bar(
            mass_total=mass_total, nnode=nnode)

        moment = np.zeros((nnode, 3))
        moment[:, 1] = 10.
        force = np.zeros((nnode, 3))
        inertia = np.zeros((nnode, 6), dtype='float64')
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment
        # assert np.allclose(dforce[0, 0], gx/2)
        # assert np.allclose(dforce[0, 1], gy/2)
        # assert np.allclose(dforce[0, 2], gz/2)
        #
        # assert np.allclose(dforce[0, :].max(), gx/2)
        # assert np.allclose(dforce[1, :].max(), gy/2)
        # assert np.allclose(dforce[2, :].max(), gz/2)

        # assert np.allclose(dforce[0, :].min(), gx/2)
        # assert np.allclose(dforce[1, :].min(), gy/2)
        # assert np.allclose(dforce[2, :].min(), gz/2)

        assert np.allclose(force_out.sum(), 0.)
        assert np.allclose(moment_out.sum(), 0.)
        # print('df:\n', dforce)
        # print('dm:\n', dmoment)
        # print('force_out:\n', force_out)
        # print('moment_out:\n', moment_out)

    def test_inrel_bar_force_linear(self):
        """inertia1 = [0, 6250, 6250, 0, 0, 0]"""
        nnode = 5
        mass_total = 5.
        xyz, mass = get_bar(
            mass_total=mass_total, nnode=nnode)
        # mass[0] = 5.

        linspace = np.linspace(-.5, .5, num=nnode)
        fx = 0.0 * linspace
        fy = 0.0 * mass
        fz = linspace
        # print(f'mass = {mass}; sum={mass.sum():g}')
        # print(f'fz = {fz}; sum={fz.sum():g}')
        force = np.column_stack([fx, fy, fz])
        nz_expected = fz.sum() / mass_total
        # print(f'nz_expected = {nz_expected:g}')
        moment = force * 0.0
        inertia = np.zeros((nnode, 6), dtype='float64')
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment
        # print('df:\n', dforce)
        # print('dm:\n', dmoment)
        assert np.allclose(force_out.sum(), 0.)
        assert np.allclose(moment_out.sum(), 0.)
        # print('force_out:\n', force_out)

    def test_inrel_bar_force_linear_rotated(self):  # pragma: no cover
        """
        inertia1 = [370.646, 4576.892, 7552.462, 0, 0, 0]
        inertia_rotated = [0, 6250, 6250, 0, 0, 0]
        """
        nnode = 5
        mass_total = 5.
        xyz, mass = get_bar(
            mass_total=mass_total, nnode=nnode)
        # mass[0] = 5.

        r = xyz[:, 0].copy()
        v = ([1., 2., 3.])
        v /= np.linalg.norm(v)
        p0 = np.zeros(3)
        xyz = p0[np.newaxis, :] + (v[np.newaxis, :] * r[:, np.newaxis])
        # r2 = np.linalg.norm(xyz, axis=1)
        # print('r2',r2)

        # theta = np.radians(45.)
        # xyz[:, 0] = r * np.cos(theta)
        # xyz[:, 1] = r * np.sin(theta)

        linspace = np.linspace(-.5, .5, num=nnode)
        fx = 0.0 * linspace
        fy = 0.0 * mass
        fz = linspace
        # print(f'mass = {mass}; sum={mass.sum():g}')
        # print(f'fz = {fz}; sum={fz.sum():g}')
        force = np.column_stack([fx, fy, fz])
        nz_expected = fz.sum() / mass_total
        # print(f'nz_expected = {nz_expected:g}')
        moment = force * 0.0
        inertia = np.zeros((nnode, 6), dtype='float64')
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment, debug=True)
        force_out = force + dforce
        moment_out = moment + dmoment
        # print('df:\n', dforce)
        print('dm:\n', dmoment)

        # [16.8622951, 4.21557378, 0, 4.21557378, 16.8622951]
        print(f'dm_norm = {np.linalg.norm(dmoment, axis=1)}')
        assert np.allclose(force_out.sum(), 0.), force_out.sum()
        assert np.allclose(moment_out.sum(), 0.), moment_out.sum()
        # print('force_out:\n', force_out)


def get_bar(mass_total: float=11.0,
            nnode: int=11):
    x = np.linspace(0., 100., num=nnode)
    xyz = np.zeros((nnode, 3), dtype='float64')
    xyz[:, 0] = x

    mass = np.ones(nnode, dtype='float64') * mass_total / nnode
    return xyz, mass


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
