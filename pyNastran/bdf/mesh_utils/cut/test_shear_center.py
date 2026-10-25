"""Closed-form checks for ``shear_center``."""
import unittest

import numpy as np

from pyNastran.bdf.mesh_utils.cut.shear_center import shear_center

E = 1.0e7
G = 4.0e6


def build_walls(segments, ex=E, gxy=G, nsub: int=8):
    """
    Turns ``[((u1, v1), (u2, v2), thickness), ...]`` into the wall soup the
    cutter would hand over: ``(nwall, 3)`` end coordinates with the in-plane
    axes in columns 0 and 2, plus per-wall length/thickness/E/G.

    Each segment is chopped into ``nsub`` pieces so the topology has to be
    recovered by welding rather than being handed over intact.
    """
    p1, p2, thickness = [], [], []
    for start, end, ti in segments:
        start = np.asarray(start, dtype='float64')
        end = np.asarray(end, dtype='float64')
        for i in range(nsub):
            p1.append(start + (end - start) * i / nsub)
            p2.append(start + (end - start) * (i + 1) / nsub)
            thickness.append(ti)
    p1 = np.array(p1, dtype='float64')
    p2 = np.array(p2, dtype='float64')
    nwall = len(p1)

    xyz1 = np.zeros((nwall, 3), dtype='float64')
    xyz2 = np.zeros((nwall, 3), dtype='float64')
    xyz1[:, 0], xyz1[:, 2] = p1[:, 0], p1[:, 1]
    xyz2[:, 0], xyz2[:, 2] = p2[:, 0], p2[:, 1]
    length = np.linalg.norm(xyz2 - xyz1, axis=1)
    return (xyz1, xyz2, length, np.array(thickness, dtype='float64'),
            np.full(nwall, ex, dtype='float64'),
            np.full(nwall, gxy, dtype='float64'))


def channel(h: float=10., b: float=4., t: float=0.1):
    """A C-section: web on ``u=0`` from -h/2 to h/2, flanges opening toward +u."""
    return [((0., -h / 2), (0., h / 2), t),
            ((0., h / 2), (b, h / 2), t),
            ((0., -h / 2), (b, -h / 2), t)]


def box(a: float=8., c: float=5., t: float=0.1):
    """A rectangular tube, 2a wide by 2c tall, centered on the origin."""
    return [((-a, -c), (a, -c), t),
            ((a, -c), (a, c), t),
            ((a, c), (-a, c), t),
            ((-a, c), (-a, -c), t)]


class TestShearCenter(unittest.TestCase):
    def test_channel(self):
        """
        The textbook thin-walled channel.  Its shear center sits outside the
        section entirely, on the far side of the web from the flanges, at::

            e = b^2 h^2 t / (4 I)        I = t h^3/12 + b t h^2/2

        measured from the web centerline.  This is the case that actually
        distinguishes a correct solve from a plausible one: a centroid-based
        answer would put it at u = +0.889, inside the section, and a sign error
        would put it at u = +1.41.
        """
        h, b, t = 10., 4., 0.1
        result = shear_center(*build_walls(channel(h, b, t)))
        self.assertEqual(result.method, 'open')
        self.assertEqual(result.ncells, 0)

        inertia = t * h ** 3 / 12. + b * t * h ** 2 / 2.
        e_exact = b ** 2 * h ** 2 * t / (4. * inertia)
        assert np.allclose(result.xy_shear_center[0], -e_exact, rtol=1e-10), \
            (result.xy_shear_center[0], -e_exact)
        # symmetric about v, so the shear center is on the axis
        assert np.allclose(result.xy_shear_center[1], 0., atol=1e-10), \
            result.xy_shear_center
        # ...and it is nowhere near the centroid, which is the whole point
        assert np.allclose(result.xy_neutral_axis[0], 0.888888888889, rtol=1e-9), \
            result.xy_neutral_axis
        assert result.force_error < 1e-10, result.force_error

    def test_channel_is_mesh_independent(self):
        """
        The strip integrals are exact for a straight wall, so chopping a wall
        into more pieces must not change the answer at all -- only the welding
        has more work to do.  A drift here means the flow integration went
        lumped somewhere.
        """
        reference = None
        for nsub in (1, 2, 5, 40):
            result = shear_center(*build_walls(channel(), nsub=nsub))
            if reference is None:
                reference = result.xy_shear_center
            assert np.allclose(result.xy_shear_center, reference, rtol=1e-12), \
                (nsub, result.xy_shear_center, reference)

    def test_angle(self):
        """Two legs meeting at a point: the shear center is that point."""
        leg, t = 6., 0.1
        segments = [((0., 0.), (leg, 0.), t), ((0., 0.), (0., leg), t)]
        result = shear_center(*build_walls(segments))
        assert np.allclose(result.xy_shear_center, 0., atol=1e-10), \
            result.xy_shear_center
        assert result.force_error < 1e-10, result.force_error

    def test_doubly_symmetric_box(self):
        """A closed doubly symmetric section: shear center == neutral axis."""
        result = shear_center(*build_walls(box()))
        self.assertEqual(result.method, 'closed')
        self.assertEqual(result.ncells, 1)
        assert np.allclose(result.xy_shear_center, 0., atol=1e-9), \
            result.xy_shear_center
        assert np.allclose(result.xy_neutral_axis, 0., atol=1e-9), \
            result.xy_neutral_axis
        assert result.force_error < 1e-10, result.force_error

    def test_two_cell_box(self):
        """A center web makes two cells; symmetry still pins the answer."""
        a, c, t = 8., 5., 0.1
        segments = box(a, c, t) + [((0., -c), (0., c), t)]
        result = shear_center(*build_walls(segments))
        self.assertEqual(result.ncells, 2)
        assert np.allclose(result.xy_shear_center, 0., atol=1e-9), \
            result.xy_shear_center
        assert result.force_error < 1e-10, result.force_error

    def test_branched_section(self):
        """
        A closed box with two dangling flanges.  The flanges carry bending
        stress but no closed-cell flow, so the solve has to drive the flow to
        zero at each free end on its own -- there is no special case for it.
        """
        a, c, t = 8., 5., 0.1
        segments = box(a, c, t) + [((a, c), (a + 3., c), t),
                                   ((-a, c), (-a - 3., c), t)]
        result = shear_center(*build_walls(segments))
        self.assertEqual(result.ncells, 1)
        # still symmetric about u=0
        assert np.allclose(result.xy_shear_center[0], 0., atol=1e-9), \
            result.xy_shear_center
        # the ears raise the centroid but drag the shear center the other way
        assert result.xy_neutral_axis[1] > 0.5, result.xy_neutral_axis
        assert result.xy_shear_center[1] < 0., result.xy_shear_center
        assert result.force_error < 1e-10, result.force_error

    def test_shear_modulus_only_matters_when_closed(self):
        """
        An open section's shear center is a pure bending quantity: scaling E or
        G uniformly cannot move it.  A closed section redistributes flow by
        wall compliance, so stiffening one wall in shear does move it.
        """
        base = shear_center(*build_walls(channel()))
        for ex, gxy in ((10. * E, G), (E, 100. * G)):
            result = shear_center(*build_walls(channel(), ex=ex, gxy=gxy))
            assert np.allclose(result.xy_shear_center, base.xy_shear_center,
                               rtol=1e-12), (ex, gxy, result.xy_shear_center)

        a, c, t = 8., 5., 0.1
        walls = list(build_walls(box(a, c, t)))
        nbottom = len(walls[0]) // 4
        walls[5] = walls[5].copy()
        walls[5][:nbottom] *= 10.          # stiffen the bottom wall in shear
        result = shear_center(*walls)
        assert np.allclose(result.xy_shear_center[0], 0., atol=1e-9), \
            result.xy_shear_center
        assert result.xy_shear_center[1] < -1., result.xy_shear_center
        assert result.force_error < 1e-10, result.force_error

    def test_heterogeneous_modulus(self):
        """
        The neutral axis is the *modulus*-weighted centroid, not the area
        centroid.  Stiffening the bottom cap in E must pull it down even though
        the geometry is unchanged.
        """
        a, c, t = 8., 5., 0.1
        walls = list(build_walls(box(a, c, t)))
        nbottom = len(walls[0]) // 4
        walls[4] = walls[4].copy()
        walls[4][:nbottom] *= 5.
        result = shear_center(*walls)
        assert result.xy_neutral_axis[1] < -2., result.xy_neutral_axis
        assert np.allclose(result.xy_neutral_axis[0], 0., atol=1e-9), \
            result.xy_neutral_axis
        assert result.force_error < 1e-10, result.force_error

    def test_degenerate_inputs(self):
        """No walls, and walls all on one line."""
        empty = (np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0),
                 np.zeros(0), np.zeros(0), np.zeros(0))
        result = shear_center(*empty)
        self.assertEqual(result.method, 'none')
        assert np.isnan(result.xy_shear_center).all()

        flat = [((0., 0.), (10., 0.), 0.1)]
        result = shear_center(*build_walls(flat))
        self.assertEqual(result.method, 'none')
        assert np.isnan(result.xy_shear_center).all()
        # the centroid is still perfectly well defined
        assert np.allclose(result.xy_neutral_axis, [5., 0.]), \
            result.xy_neutral_axis


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
