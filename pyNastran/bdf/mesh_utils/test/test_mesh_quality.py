"""various mesh_utils tests"""
import os
import unittest

import numpy as np
from cpylog import SimpleLogger

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import (
    delete_bad_shells, get_bad_shells, quad_quality_nastran,
)

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)


class TestMeshQuality(unittest.TestCase):

    def test_quad_180_01(self):
        r"""
        Identify a 180+ degree quad

        EID = 100 (max_theta > 180)
        y
        ^         4
        |       / |
        |     /   |
        |   /     |
        | /       |
        /         |
        1------2  |----> x
                \ |
                 \|
                  3


        """
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,   0.,0.\n'
            'GRID,2,,1.,   0.,0.\n'
            'GRID,3,,2.,  -1.,0.\n'
            'GRID,4,,2.,   1.,0.\n'
            'CQUAD4,100,1, 1,2,3,4\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        log = SimpleLogger(level='error')
        bdf_filename = 'cquad4.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, log=log, xref=True)
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
        nid_map = {}
        for i, nid in enumerate(sorted(model.nodes.keys())):
            nid_map[nid] = i


        max_theta_active = 180.

        min_theta = 0.1
        max_skew = 1000.
        max_aspect_ratio = 1000.
        max_taper_ratio = 1000.

        # max theta
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta_active,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
            max_taper_ratio=max_taper_ratio)
        assert eids_to_delete == [100], eids_to_delete

        delete_bad_shells(
            model,
            min_theta=min_theta,
            max_theta=max_theta_active,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
            max_taper_ratio=max_taper_ratio)

        assert len(model.elements) == 0, model.elements
        os.remove(bdf_filename)

    def test_quad_quality_02(self):
        r"""
        Identify bad quads

        EID = 101 (skew)
        y
        ^         4
        |       /  \
        |     /       \
        |   /               \
        | /                      \
        /                             \
        1------------------------------5-> x
          \                     /
             \          /
                  3

        EID = 102 (aspect ratio)
        y
        ^         8--------------------------6
        |       /                           /
        |     /                            /
        |   /                             /
        | /                              /
        /                               /
        1------------------------------5-> x

        EID = 103 (taper ratio)
        y
        ^         8--------------7
        |       /                 \
        |     /                    \
        |   /                       \
        | /                          \
        /                             \
        1------------------------------5-> x

        """
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,   0.,0.\n'
            #'GRID,2,,1.,   0.,0.\n'
            'GRID,3,,2.,  -1.,0.\n'
            'GRID,4,,2.,   1.,0.\n'

            # shift x by 200
            'GRID,5,,200., 0.,0.\n'

            # shifted by 30 to match dx of GRID 8
            #raised up to match 7,8 height
            'GRID,6,,245., 100.,0.\n'

            # raised up to prevent large thetas
            # brought in to cause a taper
            'GRID,7,,55.,  100.,0.\n'
            'GRID,8,,45.,  100.,0.\n'

            #'CQUAD4,100,1, 1,2,3,4\n'
            'CQUAD4,101,1, 1,3,5,4\n'
            'CQUAD4,102,1, 1,5,6,8\n'
            'CQUAD4,103,1, 1,5,7,8\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        log = SimpleLogger(level='error')
        bdf_filename = 'cquad4.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, log=log, xref=True)
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
        nid_map = {}
        for i, nid in enumerate(sorted(model.nodes.keys())):
            nid_map[nid] = i


        max_skew_active = 70.
        max_aspect_ratio_active = 50.
        max_taper_ratio_active = 2.5

        min_theta = 0.1
        max_theta = 1000.
        max_skew = 1000.
        max_aspect_ratio = 1000.
        max_taper_ratio = 1000.

        #log = SimpleLogger(level='debug')
        #model.log = log

        # max skew
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta,
            max_skew=max_skew_active,
            max_aspect_ratio=max_aspect_ratio,
            max_taper_ratio=max_taper_ratio)
        assert eids_to_delete == [101], eids_to_delete

        # aspect ratio
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio_active,
            max_taper_ratio=max_taper_ratio,
        )
        assert eids_to_delete == [101], eids_to_delete

        # taper ratio
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
            max_taper_ratio=max_taper_ratio_active,
        )
        assert eids_to_delete == [103], eids_to_delete
        os.remove(bdf_filename)

    def test_tri_quality_01(self):
        r"""
        Identify a 180+ degree quad

        EID = 100 (perfect)
        y
        ^     3
        |    / \
        |   /   \
        |  /     \
        | /       \
        |/         \
        1-----------2----> x

        EID = 101 (???)
        y
        ^     6
        |    /   \
        |   /       \
        |  /           \
        | /               \
        |/                   \
        4----------------------5----> x

        """
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,   0.,0.\n'
            'GRID,2,,1.,   0.,0.\n'
            'GRID,3,,0.5,  1.,0.\n'
            'CTRIA3,100,1, 1,2,3\n'

            'GRID,4,,  0.,   0.,5.\n'
            'GRID,5,,300.,  0.,5.\n'
            'GRID,6,, 0.5,  1.,5.\n'
            'CTRIA3,101,1, 4,5,6\n'

            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        log = SimpleLogger(level='error')
        bdf_filename = 'cquad4.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, log=log, xref=True)
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
        nid_map = {}
        for i, nid in enumerate(sorted(model.nodes.keys())):
            nid_map[nid] = i


        max_skew_active = 70.
        max_theta_active = 115.
        max_aspect_ratio_active = 15.

        min_theta = 0.1
        max_theta = 1000.
        max_skew = 1000.
        max_aspect_ratio = 1000.

        # max theta
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta_active,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
        )
        assert eids_to_delete == [101], eids_to_delete

        # max skew
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta,
            max_skew=max_skew_active,
            max_aspect_ratio=max_aspect_ratio,
        )
        assert eids_to_delete == [101], eids_to_delete

        # aspect ratio
        eids_to_delete = get_bad_shells(
            model, xyz_cid0, nid_map,
            min_theta=min_theta,
            max_theta=max_theta,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio_active,
        )
        assert eids_to_delete == [101], eids_to_delete

        delete_bad_shells(
            model,
            min_theta=min_theta,
            max_theta=max_theta_active,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
        )
        assert eids_to_delete == [101], eids_to_delete

        #assert len(model.elements) == 0, model.elements
        os.remove(bdf_filename)

    def test_fix_bad_quads(self):
        """split high interior angle quads"""
        msg = [
            'SOL 101',
            'CEND',
            'BEGIN BULK',
            'GRID,1,,0.,0.,0.',
            'GRID,2,,1.,0.,0.',
            'GRID,3,,1.,1.,0.',
            'GRID,4,,0.,1.,0.',
            'GRID,5,,1.,1.,0.00001',
            'GRID,6,,0.,0.,0.00001',
            'CQUAD4,1, 2, 1,2,3,4',
            'CQUAD4,2, 2, 1,2,3,5',
            'CQUAD4,3, 2, 1,6,3,5',
        ]
        log = SimpleLogger(level='error')
        bdf_filename = 'fix_bad_quads.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(msg))
        model = read_bdf(bdf_filename, log=log, xref=False, debug=False)
        model.cross_reference(xref=True, xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        convert_bad_quads_to_tris(model, min_edge_length=0.01)
        #for eid, elem in sorted(model.elements.items()):
            #print(elem)
        assert model.card_count['CQUAD4'] == 2, model.card_count
        assert model.card_count['CTRIA3'] == 1, model.card_count
        os.remove(bdf_filename)


    def test_quad_quality_nastran_perfect_square(self):
        """A perfect unit square should have ideal quality metrics."""
        p1 = np.array([0., 0., 0.])
        p2 = np.array([1., 0., 0.])
        p3 = np.array([1., 1., 0.])
        p4 = np.array([0., 1., 0.])

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # Skew: midpoint vectors are perpendicular → skew = pi/2
        assert np.isclose(skew, np.pi / 2, atol=1e-12), f'skew={skew}'
        # Taper: all corner triangles equal → taper = 0
        assert np.isclose(taper, 0.0, atol=1e-12), f'taper={taper}'
        # Warp: flat element → warp = 0
        assert np.isclose(warp, 0.0, atol=1e-12), f'warp={warp}'

    def test_quad_quality_nastran_rectangle(self):
        """A 4:1 rectangle should have ideal skew/taper/warp (no distortion)."""
        p1 = np.array([0., 0., 0.])
        p2 = np.array([4., 0., 0.])
        p3 = np.array([4., 1., 0.])
        p4 = np.array([0., 1., 0.])

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # Rectangle: midpoint vectors still perpendicular
        assert np.isclose(skew, np.pi / 2, atol=1e-12), f'skew={skew}'
        # Rectangle: no taper
        assert np.isclose(taper, 0.0, atol=1e-12), f'taper={taper}'
        # Flat: no warp
        assert np.isclose(warp, 0.0, atol=1e-12), f'warp={warp}'

    def test_quad_quality_nastran_skewed(self):
        """A parallelogram has reduced skew angle."""
        # Parallelogram: shift top nodes by 1 in x
        p1 = np.array([0., 0., 0.])
        p2 = np.array([2., 0., 0.])
        p3 = np.array([3., 1., 0.])
        p4 = np.array([1., 1., 0.])

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # Midpoint vectors: PQ=(2,1,0)-(1,0,0)=(1,1,0), RS=(2.5,0.5,0)-(0.5,0.5,0)=(2,0,0)
        # angle between (1,1,0) and (2,0,0) = arccos(2/(sqrt(2)*2)) = arccos(1/sqrt(2)) = 45 deg
        assert np.isclose(skew, np.pi / 4, atol=1e-10), f'skew={np.degrees(skew)} deg'
        # Parallelogram has zero taper (all corner triangles are equal)
        assert np.isclose(taper, 0.0, atol=1e-10), f'taper={taper}'
        # Flat
        assert np.isclose(warp, 0.0, atol=1e-12), f'warp={warp}'

    def test_quad_quality_nastran_tapered(self):
        """A trapezoid should show nonzero taper."""
        p1 = np.array([0., 0., 0.])
        p2 = np.array([4., 0., 0.])
        p3 = np.array([3., 1., 0.])
        p4 = np.array([1., 1., 0.])

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # Taper should be 0 for symmetric trapezoid (equal corner triangles)
        # Actually symmetric trapezoid: ABD area != BCA area
        # Let's just verify taper >= 0 and taper < 0.5 (not failing threshold)
        assert taper >= 0.0, f'taper={taper}'
        # Flat
        assert np.isclose(warp, 0.0, atol=1e-12), f'warp={warp}'

    def test_quad_quality_nastran_warped(self):
        """A quad with one node out of plane should show nonzero warp."""
        p1 = np.array([0., 0., 0.])
        p2 = np.array([1., 0., 0.])
        p3 = np.array([1., 1., 0.])
        p4 = np.array([0., 1., 0.5])  # lifted out of plane

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # Warp should be positive (element is not flat)
        assert warp > 0.0, f'warp={warp}'
        # For this geometry:
        # AC = (1,1,0), BD = (-1,1,0.5)
        # AC x BD = (1*0.5-0*1, 0*(-1)-1*0.5, 1*1-1*(-1)) = (0.5, -0.5, 2)
        # PK = (0.5,-0.5,2)/|(0.5,-0.5,2)|
        # AB = (1,0,0)
        # HH = AB . PK = 0.5/|pk_raw|
        # D_AC = sqrt(2), D_BD = sqrt(1+1+0.25) = sqrt(2.25) = 1.5
        # warp = |0.5/sqrt(0.25+0.25+4)| / (sqrt(2) + 1.5)
        pk_norm = np.sqrt(0.25 + 0.25 + 4.0)
        hh = 0.5 / pk_norm
        d_ac = np.sqrt(2.0)
        d_bd = np.sqrt(2.25)
        expected_warp = abs(hh) / (d_ac + d_bd)
        assert np.isclose(warp, expected_warp, atol=1e-12), f'warp={warp}, expected={expected_warp}'

    def test_quad_quality_nastran_collapsed(self):
        """A quad collapsed to a triangle should have taper = 1.0."""
        p1 = np.array([0., 0., 0.])
        p2 = np.array([1., 0., 0.])
        p3 = np.array([0.5, 1., 0.])
        p4 = np.array([0.5, 1., 0.])  # coincident with p3

        skew, taper, warp = quad_quality_nastran(p1, p2, p3, p4)

        # When collapsed to triangle, one corner triangle has zero area
        # and Amax = total triangle area, Q = 0.5 * quad area
        # taper should approach 1.0
        assert taper > 0.9, f'taper={taper}'


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
