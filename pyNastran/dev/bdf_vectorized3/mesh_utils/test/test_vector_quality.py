"""tests mesh quality utilities: get_bad_shells, delete_bad_shells, convert_bad_quads_to_tris"""
import unittest
import numpy as np

from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.mesh_utils.quality import (
    get_bad_shells, delete_bad_shells, convert_bad_quads_to_tris)

from cpylog import SimpleLogger
log = SimpleLogger(level='error')


def _make_quad_model() -> BDF:
    """Creates a model with one good quad and one bad (high-skew) quad.

    Good quad: square 1x1
    Bad quad: severely skewed parallelogram
    """
    model = BDF(debug=False, log=log)
    # good square quad
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    # bad skewed quad
    model.add_grid(5, [2., 0., 0.])
    model.add_grid(6, [3., 0., 0.])
    model.add_grid(7, [5., 0.01, 0.])
    model.add_grid(8, [4., 0.01, 0.])

    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)

    model.add_cquad4(100, 1, [1, 2, 3, 4])
    model.add_cquad4(101, 1, [5, 6, 7, 8])
    model.setup()
    return model


def _make_tria_model() -> BDF:
    """Creates a model with one good tri and one bad (near-degenerate) tri."""
    model = BDF(debug=False, log=log)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [0.5, 1., 0.])
    # bad tri: nearly collinear
    model.add_grid(4, [2., 0., 0.])
    model.add_grid(5, [3., 0., 0.])
    model.add_grid(6, [2.5, 0.001, 0.])

    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)

    model.add_ctria3(200, 1, [1, 2, 3])
    model.add_ctria3(201, 1, [4, 5, 6])
    model.setup()
    return model


def _make_degenerate_quad_model() -> BDF:
    """Creates a model with a collapsed quad (two coincident nodes)."""
    model = BDF(debug=False, log=log)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    # node 4 == node 1 (collapsed edge)

    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)

    # quad with n4 == n1 (degenerate)
    model.add_cquad4(300, 1, [1, 2, 3, 1])
    # good quad
    model.add_grid(4, [0., 1., 0.])
    model.add_cquad4(301, 1, [1, 2, 3, 4])
    model.setup()
    return model


def _make_short_edge_quad_model() -> BDF:
    """Creates a model with a quad having a very short edge."""
    model = BDF(debug=False, log=log)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [1., 1., 0.])
    model.add_grid(4, [0., 1., 0.])
    # quad with short edge (node 5 very close to node 6)
    model.add_grid(5, [2., 0., 0.])
    model.add_grid(6, [2.0001, 0., 0.])
    model.add_grid(7, [3., 1., 0.])
    model.add_grid(8, [2., 1., 0.])

    model.add_pshell(1, mid1=1, t=0.1)
    model.add_mat1(1, 3.0e7, None, 0.3)

    model.add_cquad4(400, 1, [1, 2, 3, 4])
    model.add_cquad4(401, 1, [5, 6, 7, 8])
    model.setup()
    return model


class TestGetBadShells(unittest.TestCase):
    """Tests for get_bad_shells"""

    def test_good_quad_not_flagged(self):
        """A well-shaped quad should not be returned as bad."""
        model = _make_quad_model()
        bad_eids = get_bad_shells(model, max_skew=70.)
        assert 100 not in bad_eids, f'good quad 100 should not be bad, bad_eids={bad_eids}'

    def test_bad_skew_quad_flagged(self):
        """A severely skewed quad should be flagged."""
        model = _make_quad_model()
        bad_eids = get_bad_shells(model, max_skew=70.)
        assert 101 in bad_eids, f'bad quad 101 should be flagged, bad_eids={bad_eids}'

    def test_bad_tria_flagged(self):
        """A nearly degenerate tri (very small min_theta) should be flagged."""
        model = _make_tria_model()
        # min_theta for near-collinear tri will be very small
        bad_eids = get_bad_shells(model, min_theta=1.0)
        assert 201 in bad_eids, f'bad tri 201 should be flagged, bad_eids={bad_eids}'
        assert 200 not in bad_eids, f'good tri 200 should not be bad, bad_eids={bad_eids}'

    def test_no_bad_elements_with_loose_thresholds(self):
        """With very loose thresholds, nothing should be flagged."""
        model = _make_quad_model()
        bad_eids = get_bad_shells(
            model, min_theta=0.001, max_theta=179.99,
            max_skew=179., max_aspect_ratio=1e6,
            max_taper_ratio=1e6, max_warp=179.)
        assert len(bad_eids) == 0, f'no elements should be bad, got {bad_eids}'

    def test_non_shell_elements_not_flagged(self):
        """Non-shell elements (bars, solids) should never be in the result."""
        model = BDF(debug=False, log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [0., 0., 1.])

        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        model.add_cquad4(100, 1, [1, 2, 3, 4])

        # add a bar element
        model.add_pbar(2, 1, area=1.0)
        model.add_cbar(200, 2, [1, 2], x=[0., 0., 1.], g0=None)
        model.setup()

        bad_eids = get_bad_shells(model, max_skew=0.01)
        # bar element 200 should never appear
        assert 200 not in bad_eids


class TestDeleteBadShells(unittest.TestCase):
    """Tests for delete_bad_shells"""

    def test_delete_removes_bad_element(self):
        """delete_bad_shells should remove the bad quad and keep the good one."""
        model = _make_quad_model()
        bad_eids = delete_bad_shells(model, max_skew=70.)
        assert 101 in bad_eids
        # verify element was removed
        remaining_eids = model.cquad4.element_id
        assert 101 not in remaining_eids
        assert 100 in remaining_eids

    def test_delete_no_bad_elements(self):
        """With loose thresholds, nothing should be deleted."""
        model = _make_quad_model()
        n_before = model.cquad4.n
        bad_eids = delete_bad_shells(
            model, min_theta=0.001, max_theta=179.99,
            max_skew=179., max_aspect_ratio=1e6,
            max_taper_ratio=1e6, max_warp=179.)
        assert len(bad_eids) == 0
        assert model.cquad4.n == n_before


class TestConvertBadQuadsToTris(unittest.TestCase):
    """Tests for convert_bad_quads_to_tris"""

    def test_collapsed_quad_becomes_tri(self):
        """A quad with two identical nodes should convert to a CTRIA3.

        EID 300: nodes [1,2,3,1] -> collapsed (n4==n1)
        EID 301: nodes [1,2,3,4] -> good quad, kept
        Tolerance: node coincidence (same node ID)
        """
        model = _make_degenerate_quad_model()
        n_quads_before = model.cquad4.n
        n_tris_before = model.ctria3.n

        nconverted = convert_bad_quads_to_tris(model)

        assert nconverted == 1, f'expected 1 conversion, got {nconverted}'
        assert model.cquad4.n == n_quads_before - 1
        assert model.ctria3.n == n_tris_before + 1
        # verify EID 300 is now a tri
        assert 300 in model.ctria3.element_id
        # verify EID 301 still a quad
        assert 301 in model.cquad4.element_id

    def test_short_edge_quad_converted(self):
        """A quad with a very short edge should convert to tri when min_edge_length is set.

        EID 401: edge 5-6 length = 0.0001 (< threshold 0.01)
        EID 400: all edges ~1.0 (> threshold)
        """
        model = _make_short_edge_quad_model()
        nconverted = convert_bad_quads_to_tris(model, min_edge_length=0.01)
        assert nconverted == 1
        assert 401 not in model.cquad4.element_id
        assert 400 in model.cquad4.element_id

    def test_no_conversion_when_all_good(self):
        """No quads should be converted when all edges are long enough."""
        model = _make_quad_model()
        nconverted = convert_bad_quads_to_tris(model, min_edge_length=0.0)
        assert nconverted == 0

    def test_empty_model(self):
        """Should return 0 when model has no CQUAD4."""
        model = BDF(debug=False, log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0.5, 1., 0.])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        model.add_ctria3(100, 1, [1, 2, 3])
        model.setup()

        nconverted = convert_bad_quads_to_tris(model)
        assert nconverted == 0


if __name__ == '__main__':
    unittest.main()
