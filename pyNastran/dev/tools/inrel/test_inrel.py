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

    def test_inrel_bar_force_linear_rotated(self):
        """
        Same bar as test_inrel_bar_force_linear, but rotated onto v=[1,2,3]
        so the principal axes no longer coincide with the global axes.

        This is the case that exercises the principal-frame transform: the
        products of inertia are non-zero, so S is not the identity.

        inertia1 = [370.646, 4576.892, 7552.462, 0, 0, 0]
        inertia_rotated = [0, 6250, 6250, 0, 0, 0]

        Note the bar is collinear, so the principal inertia about the bar
        axis is 0 (eigh returns ~1e-12); the degenerate axis must be skipped
        by a relative tolerance rather than an exact != 0 test.
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
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment

        # the net force and moment must both vanish component-by-component,
        # not just in total (a plain .sum() can cancel across components)
        assert np.allclose(force_out.sum(axis=0), 0.), force_out.sum(axis=0)
        assert np.allclose(moment_out.sum(axis=0), 0.), moment_out.sum(axis=0)

        # the applied load is self-equilibrating in force, so the linear
        # relief is zero and the whole balance is carried by alpha
        assert np.allclose(dforce, 0.), dforce

        # the relief must be a genuine rotation about the bar, not all zeros
        assert np.abs(dmoment).max() > 0., dmoment


    def test_inrel_rotated_with_self_inertia(self):
        """
        Rotated 3-D mass cloud carrying non-zero per-node self-inertia,
        including products of inertia.

        The self-inertia is supplied in the global frame while alpha is solved
        in the principal frame, so each node's tensor must be rotated before
        being contracted with alpha. Applying the raw diagonal instead mixes
        frames and leaves a visible moment residual.
        """
        rng = np.random.default_rng(11)
        nnode = 6
        xyz = rng.normal(scale=5.0, size=(nnode, 3))
        mass = np.ones(nnode)

        inertia = np.zeros((nnode, 6), dtype='float64')
        inertia[:, :3] = rng.uniform(1.0, 4.0, (nnode, 3))   # Ixx, Iyy, Izz
        inertia[:, 3:] = rng.uniform(-1.0, 1.0, (nnode, 3))  # Ixy, Ixz, Iyz

        force = rng.normal(size=(nnode, 3))
        moment = rng.normal(size=(nnode, 3))
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment

        assert np.allclose(force_out.sum(axis=0), 0.), force_out.sum(axis=0)
        assert np.allclose(moment_out.sum(axis=0), 0.), moment_out.sum(axis=0)

    def test_inrel_near_collinear_degenerate_axis(self):
        """
        A *nearly* collinear bar loaded by a moment about the bar axis.

        The principal inertia about the bar axis is zero in exact arithmetic,
        but eigh returns a tiny value of either sign (here ~-7e-13) rather
        than exactly 0. An exact ``inertia != 0`` guard admits that value and
        computes alpha = M/~0, producing an enormous bogus relief; the guard
        must instead be a tolerance relative to the largest principal inertia.

        This load case is not physically solvable: there is no inertia about
        the bar axis, so no angular acceleration can react the applied axial
        moment. The correct behavior is to warn, apply no relief about that
        axis, and leave the axial moment visibly unbalanced, rather than to
        manufacture a huge alpha that corrupts the other two axes as well.
        """
        nnode = 5
        v = np.array([1., 2., 3.])
        v /= np.linalg.norm(v)
        xyz = v[np.newaxis, :] * np.linspace(0., 100., nnode)[:, np.newaxis]
        # nudge off perfectly collinear so eigh cannot return an exact 0
        xyz = xyz + np.array([[0., 0., 1e-7]]) * np.arange(nnode)[:, np.newaxis]
        mass = np.ones(nnode)

        force = np.zeros((nnode, 3))
        # a pure moment along the bar axis -- about the zero-inertia axis
        moment_per_node = 10.
        moment = np.tile((v * moment_per_node)[np.newaxis, :], (nnode, 1))
        inertia = np.zeros((nnode, 6), dtype='float64')
        with self.assertWarns(RuntimeWarning):
            dforce, dmoment = inertia_relief(
                mass, xyz, inertia,
                force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment

        # the giveaway for a divide-by-~0 is a huge relief term
        assert np.abs(dmoment).max() < 1e-3, np.abs(dmoment).max()

        # no force was applied, so there is nothing to relieve linearly
        assert np.allclose(force_out.sum(axis=0), 0.), force_out.sum(axis=0)

        # the unreacted moment must be exactly the applied axial moment: the
        # degenerate axis is skipped, and the other two axes are untouched
        moment_expected = v * moment_per_node * nnode
        assert np.allclose(moment_out.sum(axis=0), moment_expected), (
            moment_out.sum(axis=0), moment_expected)

    def test_inrel_large_forces(self):
        """
        Guards the removal of a hardcoded ``assert abs(f).max() < 1.``, which
        tripped on any realistic load magnitude.
        """
        rng = np.random.default_rng(5)
        nnode = 5
        xyz = rng.normal(scale=5.0, size=(nnode, 3))
        mass = np.ones(nnode)
        force = rng.normal(scale=1e5, size=(nnode, 3))
        moment = np.zeros((nnode, 3))
        inertia = np.zeros((nnode, 6), dtype='float64')
        dforce, dmoment = inertia_relief(
            mass, xyz, inertia,
            force, moment)
        force_out = force + dforce
        moment_out = moment + dmoment
        # scale the tolerance to the load magnitude
        assert np.abs(force_out.sum(axis=0)).max() < 1e-6, force_out.sum(axis=0)
        assert np.abs(moment_out.sum(axis=0)).max() < 1e-6, moment_out.sum(axis=0)


def get_bar(mass_total: float=11.0,
            nnode: int=11):
    x = np.linspace(0., 100., num=nnode)
    xyz = np.zeros((nnode, 3), dtype='float64')
    xyz[:, 0] = x

    mass = np.ones(nnode, dtype='float64') * mass_total / nnode
    return xyz, mass


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
