"""Tests for inertia_loads."""
import unittest

import numpy as np
from pyNastran.dev.tools.inrel.inertia_loads import (
    inertia_loads, resultant, expected_resultant)


def get_cloud(nnode: int = 7, seed: int = 12, self_inertia: bool = True):
    """a random, deliberately off-axis mass cloud"""
    rng = np.random.default_rng(seed)
    mass = rng.uniform(0.5, 4.0, nnode)
    cg = rng.uniform(-10.0, 10.0, (nnode, 3))
    if self_inertia:
        inertia = np.zeros((nnode, 6), dtype='float64')
        # diagonal terms must dominate for a physical tensor
        inertia[:, :3] = rng.uniform(1.0, 5.0, (nnode, 3))
        inertia[:, 3:] = rng.uniform(-0.3, 0.3, (nnode, 3))
    else:
        inertia = np.zeros((nnode, 6), dtype='float64')
    return mass, cg, inertia


class TestInertiaLoads(unittest.TestCase):
    def test_pure_translation(self):
        """F = -m*a at every mass; no moment, no dependence on ref_xyz"""
        mass, cg, inertia = get_cloud()
        accel = np.array([1.0, -2.0, 3.0])
        alpha = np.zeros(3)
        ref_xyz = np.array([100.0, -50.0, 7.0])  # arbitrary; must not matter

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)

        expected = -mass[:, np.newaxis] * accel[np.newaxis, :]
        assert np.allclose(force, expected), force - expected
        assert np.allclose(moment, 0.0), moment

    def test_translation_ref_independent(self):
        """with alpha=0 the answer cannot depend on ref_xyz"""
        mass, cg, inertia = get_cloud()
        accel = np.array([0.0, 0.0, -2.0])
        alpha = np.zeros(3)

        force1, moment1 = inertia_loads(
            mass, cg, inertia, accel, alpha, np.zeros(3))
        force2, moment2 = inertia_loads(
            mass, cg, inertia, accel, alpha, np.array([13.0, -4.0, 9.0]))
        assert np.allclose(force1, force2)
        assert np.allclose(moment1, moment2)

    def test_pure_rotation_single_mass(self):
        """a point mass on the +x axis, alpha about z -> force in -y"""
        mass = np.array([2.0])
        cg = np.array([[3.0, 0.0, 0.0]])
        inertia = np.zeros((1, 6))
        accel = np.zeros(3)
        alpha = np.array([0.0, 0.0, 1.0])
        ref_xyz = np.zeros(3)

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)
        # a = alpha x r = [0,0,1] x [3,0,0] = [0,3,0]; F = -m*a
        assert np.allclose(force, [[0.0, -6.0, 0.0]]), force
        assert np.allclose(moment, 0.0), moment

    def test_self_inertia_moment(self):
        """M = -[I] alpha, with the packed off-diagonals negated"""
        mass = np.zeros(1)
        cg = np.zeros((1, 3))
        # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
        inertia = np.array([[10.0, 20.0, 30.0, 1.0, 2.0, 3.0]])
        alpha = np.array([1.0, 0.0, 0.0])

        unused_force, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), alpha, np.zeros(3))
        # first column of [[10,-1,-2],[-1,20,-3],[-2,-3,30]], negated
        assert np.allclose(moment, [[-10.0, 1.0, 2.0]]), moment

    def test_resultant_matches_closed_form(self):
        """the distributed load must sum to the rigid-body resultant"""
        mass, cg, inertia = get_cloud()
        accel = np.array([1.5, -0.5, -9.0])
        alpha = np.array([0.3, -1.2, 0.7])
        ref_xyz = np.array([2.0, -3.0, 1.0])

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)
        force_total, moment_total = resultant(force, moment, cg, ref_xyz)
        force_expected, moment_expected = expected_resultant(
            mass, cg, inertia, accel, alpha, ref_xyz)

        assert np.allclose(force_total, force_expected), (
            force_total, force_expected)
        assert np.allclose(moment_total, moment_expected), (
            moment_total, moment_expected)

    def test_resultant_about_cg(self):
        """same check with ref_xyz at the cg, where sum(F) = -m*accel exactly"""
        mass, cg, inertia = get_cloud()
        ref_xyz = (mass[:, np.newaxis] * cg).sum(axis=0) / mass.sum()
        accel = np.array([0.0, 0.0, -3.0])
        alpha = np.array([0.0, 2.0, 0.0])

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)
        force_total, moment_total = resultant(force, moment, cg, ref_xyz)
        force_expected, moment_expected = expected_resultant(
            mass, cg, inertia, accel, alpha, ref_xyz)

        assert np.allclose(force_total, -mass.sum() * accel), force_total
        assert np.allclose(force_total, force_expected)
        assert np.allclose(moment_total, moment_expected), (
            moment_total, moment_expected)

    def test_offset_couple(self):
        """reporting at the grid adds (cg - grid) x F"""
        mass = np.array([5.0])
        xyz = np.array([[0.0, 0.0, 0.0]])       # the grid
        cg = np.array([[0.0, 2.0, 0.0]])        # CONM2 offset in +y
        inertia = np.zeros((1, 6))
        accel = np.array([0.0, 0.0, -1.0])

        force_cg, moment_cg = inertia_loads(
            mass, cg, inertia, accel, np.zeros(3), np.zeros(3))
        force_grid, moment_grid = inertia_loads(
            mass, cg, inertia, accel, np.zeros(3), np.zeros(3), xyz=xyz)

        assert np.allclose(force_cg, force_grid), 'force must not move'
        assert np.allclose(moment_cg, 0.0), moment_cg
        # F = [0,0,5]; r = cg-grid = [0,2,0]; r x F = [10,0,0]
        assert np.allclose(moment_grid, [[10.0, 0.0, 0.0]]), moment_grid

    def test_offset_resultant_invariant(self):
        """the resultant is the same whether loads are reported at cg or grid"""
        mass, cg, inertia = get_cloud()
        rng = np.random.default_rng(3)
        xyz = cg + rng.uniform(-1.0, 1.0, cg.shape)  # grids offset from the cgs
        accel = np.array([2.0, 1.0, -4.0])
        alpha = np.array([-0.4, 0.9, 0.2])
        ref_xyz = np.array([1.0, 1.0, 1.0])

        force_cg, moment_cg = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)
        force_grid, moment_grid = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz, xyz=xyz)

        total_cg = resultant(force_cg, moment_cg, cg, ref_xyz)
        total_grid = resultant(force_grid, moment_grid, xyz, ref_xyz)

        assert np.allclose(total_cg[0], total_grid[0]), (total_cg[0], total_grid[0])
        assert np.allclose(total_cg[1], total_grid[1]), (total_cg[1], total_grid[1])

    def test_wtmass(self):
        """wtmass scales the whole 6x6, so it scales force and moment alike"""
        mass, cg, inertia = get_cloud()
        accel = np.array([0.0, 0.0, -386.1])
        alpha = np.array([0.1, 0.2, 0.3])
        ref_xyz = np.zeros(3)
        wtmass = 1.0 / 386.1

        force1, moment1 = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz, wtmass=wtmass)
        force2, moment2 = inertia_loads(
            mass * wtmass, cg, inertia * wtmass, accel, alpha, ref_xyz)
        assert np.allclose(force1, force2)
        assert np.allclose(moment1, moment2)

    def test_subset_of_nodes(self):
        """passing a subset of rows gives those rows' loads unchanged"""
        mass, cg, inertia = get_cloud()
        accel = np.array([1.0, 2.0, -3.0])
        alpha = np.array([0.5, 0.0, -0.25])
        ref_xyz = np.array([1.0, 0.0, 0.0])

        force_all, moment_all = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)

        i = np.array([1, 4, 5])
        force_sub, moment_sub = inertia_loads(
            mass[i], cg[i], inertia[i], accel, alpha, ref_xyz)

        assert np.allclose(force_sub, force_all[i]), force_sub
        assert np.allclose(moment_sub, moment_all[i]), moment_sub

    def test_zero_acceleration(self):
        """no acceleration, no load"""
        mass, cg, inertia = get_cloud()
        force, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), np.zeros(3))
        assert np.allclose(force, 0.0)
        assert np.allclose(moment, 0.0)

    def test_omega_none_is_zero(self):
        """omega=None must be identical to omega=0"""
        mass, cg, inertia = get_cloud()
        accel = np.array([1.0, -1.0, 2.0])
        alpha = np.array([0.2, 0.3, -0.1])
        ref_xyz = np.array([1.0, 2.0, 3.0])

        force1, moment1 = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz)
        force2, moment2 = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz, omega=np.zeros(3))
        assert np.allclose(force1, force2)
        assert np.allclose(moment1, moment2)

    def test_steady_roll_centrifugal(self):
        """
        steady roll about x: a mass at y=R sees a_y = -p^2 R (inward), so the
        d'Alembert force is +p^2*m*R (outward), and x is untouched
        """
        mass = np.array([3.0])
        cg = np.array([[7.0, 2.0, 0.0]])  # x offset must not matter
        inertia = np.zeros((1, 6))
        rate = 4.0

        force, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), np.zeros(3),
            omega=np.array([rate, 0.0, 0.0]))

        expected = 3.0 * rate ** 2 * 2.0  # m p^2 R, outward = +y
        assert np.allclose(force, [[0.0, expected, 0.0]]), force
        assert np.allclose(moment, 0.0), moment

    def test_centrifugal_is_quadratic_in_rate(self):
        """doubling omega quadruples the centrifugal force"""
        mass, cg, inertia = get_cloud()
        ref_xyz = np.zeros(3)
        omega = np.array([0.0, 0.0, 1.5])

        force1, unused = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), ref_xyz, omega=omega)
        force2, unused = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), ref_xyz,
            omega=2.0 * omega)
        assert np.allclose(force2, 4.0 * force1), (force2, force1)

    def test_centrifugal_is_radial(self):
        """
        the centrifugal force has no component along the rotation axis and
        points away from the axis
        """
        mass, cg, inertia = get_cloud()
        axis = np.array([0.0, 0.0, 1.0])
        force, unused = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), np.zeros(3),
            omega=3.0 * axis)

        assert np.allclose(force @ axis, 0.0), 'axial component must vanish'
        # radial vector from the axis, and the outward projection
        radial = cg.copy()
        radial[:, 2] = 0.0
        outward = np.einsum('ni,ni->n', force, radial)
        assert (outward > 0).all(), outward

    def test_rate_coupling_zero_on_principal_axis(self):
        """omega x ([I] omega) vanishes when omega is along a principal axis"""
        mass = np.zeros(1)
        cg = np.zeros((1, 3))
        inertia = np.array([[10.0, 20.0, 30.0, 0.0, 0.0, 0.0]])  # diagonal

        for axis in np.eye(3):
            unused, moment = inertia_loads(
                mass, cg, inertia, np.zeros(3), np.zeros(3), np.zeros(3),
                omega=2.0 * axis)
            assert np.allclose(moment, 0.0), (axis, moment)

    def test_rate_coupling_nonzero_off_principal(self):
        """
        a rate off the principal axes produces a couple; check it against the
        hand-computed omega x ([I] omega)
        """
        mass = np.zeros(1)
        cg = np.zeros((1, 3))
        inertia = np.array([[10.0, 20.0, 30.0, 0.0, 0.0, 0.0]])
        omega = np.array([1.0, 1.0, 0.0])

        unused, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), np.zeros(3),
            omega=omega)
        # [I]w = [10, 20, 0]; w x [I]w = [1,1,0] x [10,20,0] = [0,0,10]
        assert np.allclose(moment, [[0.0, 0.0, -10.0]]), moment

    def test_rate_coupling_identity_brute_force(self):
        """
        the docstring claims sum(m * r x (omega x (omega x r))) equals
        omega x ([I_parallel] omega).  Verify by brute force rather than
        trusting the identity -- this is what lets expected_resultant stand as
        an independent check.
        """
        mass, cg, inertia = get_cloud(self_inertia=False)
        ref_xyz = np.array([1.0, -2.0, 0.5])
        omega = np.array([0.7, -1.3, 0.4])

        radius = cg - ref_xyz
        # direct per-mass sum of r x (m * a_centrifugal)
        accel = np.cross(omega, np.cross(omega, radius))
        brute = np.cross(radius, mass[:, np.newaxis] * accel).sum(axis=0)

        # the closed form
        eye = np.eye(3)
        rdotr = np.einsum('ni,ni->n', radius, radius)
        imat = (np.einsum('n,n,ij->ij', mass, rdotr, eye)
                - np.einsum('n,ni,nj->ij', mass, radius, radius))
        closed = np.cross(omega, imat @ omega)

        assert np.allclose(brute, closed), (brute, closed)

    def test_steady_roll_resultant(self):
        """alpha = 0, omega != 0: the pure steady-roll case"""
        mass, cg, inertia = get_cloud()
        ref_xyz = np.array([1.0, 0.5, -2.0])
        omega = np.array([5.0, 0.0, 0.0])

        force, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), ref_xyz, omega=omega)
        total = resultant(force, moment, cg, ref_xyz)
        expected = expected_resultant(
            mass, cg, inertia, np.zeros(3), np.zeros(3), ref_xyz, omega=omega)

        assert np.allclose(total[0], expected[0]), (total[0], expected[0])
        assert np.allclose(total[1], expected[1]), (total[1], expected[1])

    def test_steady_turn_resultant(self):
        """a coupled rate plus a lateral acceleration, all terms active"""
        mass, cg, inertia = get_cloud()
        ref_xyz = np.array([2.0, -3.0, 1.0])
        accel = np.array([0.0, 6.0, -9.8])
        alpha = np.array([0.4, -0.2, 0.9])
        omega = np.array([0.8, 0.3, -1.1])

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz, omega=omega)
        total = resultant(force, moment, cg, ref_xyz)
        expected = expected_resultant(
            mass, cg, inertia, accel, alpha, ref_xyz, omega=omega)

        assert np.allclose(total[0], expected[0]), (total[0], expected[0])
        assert np.allclose(total[1], expected[1]), (total[1], expected[1])

    def test_omega_resultant_with_offset_grids(self):
        """all rate terms active and moments reported at offset grids"""
        mass, cg, inertia = get_cloud()
        rng = np.random.default_rng(9)
        xyz = cg + rng.uniform(-1.0, 1.0, cg.shape)
        ref_xyz = np.array([0.5, 0.5, 0.5])
        accel = np.array([1.0, -2.0, 3.0])
        alpha = np.array([0.1, 0.5, -0.3])
        omega = np.array([-0.6, 1.4, 0.2])

        force, moment = inertia_loads(
            mass, cg, inertia, accel, alpha, ref_xyz, omega=omega, xyz=xyz)
        total = resultant(force, moment, xyz, ref_xyz)
        expected = expected_resultant(
            mass, cg, inertia, accel, alpha, ref_xyz, omega=omega)

        assert np.allclose(total[0], expected[0]), (total[0], expected[0])
        assert np.allclose(total[1], expected[1]), (total[1], expected[1])

    def test_omega_wtmass(self):
        """wtmass scales the rate terms too"""
        mass, cg, inertia = get_cloud()
        ref_xyz = np.zeros(3)
        omega = np.array([1.0, -2.0, 0.5])
        alpha = np.array([0.3, 0.1, 0.2])
        wtmass = 1.0 / 386.1

        force1, moment1 = inertia_loads(
            mass, cg, inertia, np.zeros(3), alpha, ref_xyz,
            omega=omega, wtmass=wtmass)
        force2, moment2 = inertia_loads(
            mass * wtmass, cg, inertia * wtmass, np.zeros(3), alpha, ref_xyz,
            omega=omega)
        assert np.allclose(force1, force2)
        assert np.allclose(moment1, moment2)

    def test_shape_errors(self):
        mass, cg, inertia = get_cloud(nnode=4)
        accel = np.zeros(3)
        alpha = np.zeros(3)
        ref = np.zeros(3)
        with self.assertRaises(ValueError):
            inertia_loads(mass, cg[:3], inertia, accel, alpha, ref)
        with self.assertRaises(ValueError):
            inertia_loads(mass, cg, inertia[:, :3], accel, alpha, ref)
        with self.assertRaises(ValueError):
            inertia_loads(mass, cg, inertia, np.zeros(2), alpha, ref)
        with self.assertRaises(ValueError):
            inertia_loads(mass, cg, inertia, accel, alpha, ref, xyz=cg[:2])
        with self.assertRaises(ValueError):
            inertia_loads(mass, cg, inertia, accel, alpha, ref,
                          omega=np.zeros(4))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
