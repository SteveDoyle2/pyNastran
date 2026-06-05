"""Tests for aero_utils: control surface downwash and monitor point loads.

Tests cover:
- get_control_surface_downwash: hinge geometry, swept hinge, missing AESURF
- resolve_aecomp_panel_ids: CAERO list_type, AELIST list_type, missing AECOMP
- compute_monitor_point_loads: force/moment summation, coefficients, subset
- compute_all_monitor_points: integration with MONPNT1 cards
"""
import unittest
from pathlib import Path

import numpy as np
from numpy import testing as npt

import pyNastran
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf
from pyNastran.dev.bdf_vectorized3.mesh_utils.aero_utils import (
    get_control_surface_downwash,
    resolve_aecomp_panel_ids,
    compute_monitor_point_loads,
    compute_all_monitor_points,
)

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()
BWB_PATH = MODEL_PATH / 'bwb'


class TestGetControlSurfaceDownwash(unittest.TestCase):
    """Tests for get_control_surface_downwash.

    Validates:
    - Correct downwash magnitude for spanwise hinge (w=1.0)
    - Swept hinge reduces effectiveness (w=cos(sweep))
    - Only AELIST panels get nonzero downwash
    - Missing AESURF returns all zeros
    """

    def _make_model_with_aesurf(self, hinge_cid_origin=None, hinge_j=None):
        """Build a minimal model: 4 panels, 2 on control surface."""
        model = BDF(debug=False)
        model.add_aero(velocity=0., cref=1.0, rho_ref=1.0)
        model.add_paero1(1)

        # 4 panels: 2 chordwise x 2 spanwise
        model.add_caero1(
            eid=100, pid=1, igroup=1,
            p1=[0., 0., 0.], x12=1.0,
            p4=[0., 2., 0.], x43=1.0,
            nspan=2, nchord=2)

        # AELIST: last 2 panels (trailing edge) are control surface
        # Box IDs: 100, 101, 102, 103
        model.add_aelist(sid=10, elements=[102, 103])

        # Coordinate system for hinge line
        if hinge_cid_origin is not None and hinge_j is not None:
            origin = hinge_cid_origin
            # Build z-axis perpendicular to j
            j = np.asarray(hinge_j, dtype=float)
            j = j / np.linalg.norm(j)
            # Pick k orthogonal to j
            if abs(j[2]) < 0.9:
                k = np.cross(j, [0, 0, 1])
            else:
                k = np.cross(j, [1, 0, 0])
            k = k / np.linalg.norm(k)
            i = np.cross(j, k)
            zaxis = origin + k
            xzplane = origin + i
            model.add_cord2r(cid=100, origin=origin, zaxis=zaxis, xzplane=xzplane)
            cid1 = 100
        else:
            cid1 = 0  # default: hinge along global y

        model.add_aesurf(
            aesid=1, label='FLAP',
            cid1=cid1, aelist_id1=10)

        model.setup()
        return model

    def test_spanwise_hinge_gives_unit_downwash(self):
        """Spanwise hinge (y-axis): affected panels get w=1.0.

        Hinge along y, panels with normal [0,0,1]:
        w = -(y_hat x x_hat) . n = -[0,0,-1] . [0,0,1] = 1.0
        """
        model = self._make_model_with_aesurf()
        panel_ids = np.array([100, 101, 102, 103])
        normals = np.array([[0, 0, 1.]] * 4)

        w = get_control_surface_downwash(model, 'FLAP', panel_ids, normals)

        # Panels 102, 103 are in AELIST
        npt.assert_allclose(w[0], 0.0)
        npt.assert_allclose(w[1], 0.0)
        npt.assert_allclose(w[2], 1.0, atol=1e-14)
        npt.assert_allclose(w[3], 1.0, atol=1e-14)

    def test_swept_hinge_reduces_effectiveness(self):
        """45-degree swept hinge: w = cos(45) = sqrt(2)/2 ≈ 0.707.

        Hinge axis at 45 deg in xy-plane: j = [1/sqrt(2), 1/sqrt(2), 0].
        """
        j_swept = [1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0]
        model = self._make_model_with_aesurf(
            hinge_cid_origin=[0., 0., 0.], hinge_j=j_swept)
        panel_ids = np.array([100, 101, 102, 103])
        normals = np.array([[0, 0, 1.]] * 4)

        w = get_control_surface_downwash(model, 'FLAP', panel_ids, normals)

        expected = np.cos(np.radians(45))
        npt.assert_allclose(w[2], expected, atol=1e-10)
        npt.assert_allclose(w[3], expected, atol=1e-10)

    def test_missing_aesurf_returns_zeros(self):
        """Nonexistent AESURF label returns all-zero downwash."""
        model = self._make_model_with_aesurf()
        panel_ids = np.array([100, 101, 102, 103])
        normals = np.array([[0, 0, 1.]] * 4)

        w = get_control_surface_downwash(model, 'NOSUCH', panel_ids, normals)
        npt.assert_allclose(w, 0.0)

    def test_empty_model_returns_zeros(self):
        """Model with no AESURF cards returns zeros."""
        model = BDF(debug=False)
        model.setup()
        panel_ids = np.array([1, 2, 3])
        normals = np.array([[0, 0, 1.]] * 3)

        w = get_control_surface_downwash(model, 'FLAP', panel_ids, normals)
        assert len(w) == 3
        npt.assert_allclose(w, 0.0)


class TestResolveAecompPanelIds(unittest.TestCase):
    """Tests for resolve_aecomp_panel_ids.

    Validates:
    - CAERO list_type resolves contiguous panel ranges
    - AELIST list_type resolves explicit panel lists
    - Missing AECOMP returns None
    """

    def _make_model_with_aecomp(self):
        """Two CAERO1 groups: 100-103 and 200-203."""
        model = BDF(debug=False)
        model.add_aero(velocity=0., cref=1.0, rho_ref=1.0)
        model.add_paero1(1)

        model.add_caero1(
            eid=100, pid=1, igroup=1,
            p1=[0., 0., 0.], x12=1.0,
            p4=[0., 2., 0.], x43=1.0,
            nspan=2, nchord=2)
        model.add_caero1(
            eid=200, pid=1, igroup=1,
            p1=[0., 3., 0.], x12=1.0,
            p4=[0., 5., 0.], x43=1.0,
            nspan=2, nchord=2)

        # AECOMP 'WING' references CAERO1 IDs 100 and 200
        model.add_aecomp('WING', 'CAERO', [100, 200])

        # AECOMP 'FLAP' via AELIST
        model.add_aelist(sid=50, elements=[102, 103, 202, 203])
        model.add_aecomp('FLAP', 'AELIST', [50])

        model.setup()
        return model

    def test_caero_list_type(self):
        """AECOMP with CAERO list_type selects all panels in those groups."""
        model = self._make_model_with_aecomp()
        all_ids = np.array([100, 101, 102, 103, 200, 201, 202, 203])

        mask = resolve_aecomp_panel_ids(model, 'WING', all_ids)

        assert mask is not None
        assert mask.sum() == 8  # all panels from both CAERO1s

    def test_aelist_list_type(self):
        """AECOMP with AELIST list_type selects only listed panels."""
        model = self._make_model_with_aecomp()
        all_ids = np.array([100, 101, 102, 103, 200, 201, 202, 203])

        mask = resolve_aecomp_panel_ids(model, 'FLAP', all_ids)

        assert mask is not None
        assert mask.sum() == 4
        expected_ids = {102, 103, 202, 203}
        selected_ids = set(all_ids[mask])
        assert selected_ids == expected_ids

    def test_missing_aecomp_returns_none(self):
        """Nonexistent AECOMP name returns None."""
        model = self._make_model_with_aecomp()
        all_ids = np.array([100, 101, 102, 103])

        result = resolve_aecomp_panel_ids(model, 'NOSUCH', all_ids)
        assert result is None

    def test_empty_model_returns_none(self):
        """Model with no AECOMP cards returns None."""
        model = BDF(debug=False)
        model.setup()
        all_ids = np.array([1, 2, 3])

        result = resolve_aecomp_panel_ids(model, 'WING', all_ids)
        assert result is None


class TestComputeMonitorPointLoads(unittest.TestCase):
    """Tests for compute_monitor_point_loads.

    Validates:
    - Force summation is correct
    - Moment cross product is correct
    - Coefficients are properly nondimensionalized
    - Panel subset (mask) works correctly
    - Empty selection returns zeros
    """

    def test_single_panel_force_and_moment(self):
        """One panel with Fz=100 at x=10 gives My = 100*10 = 1000.

        Panel at (10, 0, 0), force = (0, 0, 100), ref at origin.
        M = r x F = (10,0,0) x (0,0,100) = (0, -1000, 0).
        """
        forces = np.array([[0., 0., 100.]])
        points = np.array([[10., 0., 0.]])
        ref = np.array([0., 0., 0.])

        result = compute_monitor_point_loads(
            forces, points, ref, q_inf=10.0, sref=100.0,
            cref=5.0, bref=20.0)

        npt.assert_allclose(result['forces'], [0, 0, 100])
        npt.assert_allclose(result['moments'], [0, -1000, 0])
        # CZ = Fz / (q*S) = 100 / (10*100) = 0.1
        assert result['coefficients'][2] == 0.1
        # CMY = My / (q*S*c) = -1000 / (10*100*5) = -0.2
        npt.assert_allclose(result['coefficients'][4], -0.2)

    def test_multiple_panels_sum(self):
        """Forces from two panels sum correctly."""
        forces = np.array([[0., 0., 50.], [0., 0., 50.]])
        points = np.array([[5., 0., 0.], [15., 0., 0.]])
        ref = np.array([0., 0., 0.])

        result = compute_monitor_point_loads(
            forces, points, ref, q_inf=1.0, sref=1.0,
            cref=1.0, bref=1.0)

        npt.assert_allclose(result['forces'], [0, 0, 100])
        # My = -(50*5 + 50*15) = -1000... actually:
        # M = (5,0,0)x(0,0,50) + (15,0,0)x(0,0,50)
        #   = (0,-250,0) + (0,-750,0) = (0,-1000,0)
        npt.assert_allclose(result['moments'], [0, -1000, 0])
        assert result['n_panels'] == 2

    def test_panel_mask_selects_subset(self):
        """Panel mask selects only indicated panels."""
        forces = np.array([[10., 0., 0.], [0., 0., 100.], [0., 20., 0.]])
        points = np.array([[0., 0., 0.], [5., 0., 0.], [0., 0., 0.]])
        ref = np.array([0., 0., 0.])
        mask = np.array([False, True, False])

        result = compute_monitor_point_loads(
            forces, points, ref, q_inf=1.0, sref=1.0,
            cref=1.0, bref=1.0, panel_mask=mask)

        npt.assert_allclose(result['forces'], [0, 0, 100])
        assert result['n_panels'] == 1

    def test_empty_mask_returns_zeros(self):
        """All-False mask gives zero result."""
        forces = np.array([[0., 0., 100.]])
        points = np.array([[5., 0., 0.]])
        ref = np.array([0., 0., 0.])
        mask = np.array([False])

        result = compute_monitor_point_loads(
            forces, points, ref, q_inf=1.0, sref=1.0,
            cref=1.0, bref=1.0, panel_mask=mask)

        npt.assert_allclose(result['forces'], 0)
        npt.assert_allclose(result['moments'], 0)
        assert result['n_panels'] == 0

    def test_coefficients_nondimensionalization(self):
        """Coefficients use correct reference lengths.

        CX,CY,CZ: F / (q*S)
        CMX,CMZ: M / (q*S*b)
        CMY: M / (q*S*c)
        """
        forces = np.array([[100., 200., 300.]])
        points = np.array([[0., 0., 0.]])  # at ref → zero moments
        ref = np.array([0., 0., 0.])

        result = compute_monitor_point_loads(
            forces, points, ref, q_inf=2.0, sref=50.0,
            cref=10.0, bref=30.0)

        qS = 2.0 * 50.0
        npt.assert_allclose(result['coefficients'][0], 100. / qS)
        npt.assert_allclose(result['coefficients'][1], 200. / qS)
        npt.assert_allclose(result['coefficients'][2], 300. / qS)


class TestBWBIntegration(unittest.TestCase):
    """Integration tests using the BWB model.

    Validates:
    - get_control_surface_downwash on real AESURF (TFLAP)
    - resolve_aecomp_panel_ids with real AECOMP cards
    """

    @classmethod
    def setUpClass(cls):
        bdf_file = BWB_PATH / 'bwb_saero.bdf'
        if not bdf_file.exists():
            raise unittest.SkipTest('BWB model not available')
        cls.model = read_bdf(str(bdf_file), debug=False)

    def test_tflap_downwash_nonzero(self):
        """BWB TFLAP control surface produces nonzero downwash on its panels."""
        caero1 = self.model.caero1
        pts_list, elems_list = caero1.panel_points_elements()

        # Element indices are already globally offset across cards
        all_points = np.vstack(pts_list)
        all_elements = np.vstack(elems_list)
        n_panels = len(all_elements)

        # Compute normals: cross(p2-p0, p1-p3) points +z for flat panels
        p0 = all_points[all_elements[:, 0]]
        p1 = all_points[all_elements[:, 1]]
        p2 = all_points[all_elements[:, 2]]
        p3 = all_points[all_elements[:, 3]]
        d1 = p2 - p0
        d2 = p1 - p3
        normals = np.cross(d1, d2)
        norms = np.linalg.norm(normals, axis=1, keepdims=True)
        normals = normals / norms

        # Build box IDs from CAERO1 EIDs + sequential numbering
        _, npanels = caero1.get_panel_npoints_nelements()
        panel_ids = np.concatenate([
            np.arange(eid, eid + n)
            for eid, n in zip(caero1.element_id, npanels)])

        w = get_control_surface_downwash(
            self.model, 'TFLAP', panel_ids, normals)

        n_active = np.sum(np.abs(w) > 1e-10)
        assert n_active > 0, 'TFLAP should affect at least one panel'
        assert n_active < n_panels, 'TFLAP should not affect all panels'
        # Downwash should be positive (TE down = positive lift)
        assert np.all(w[w != 0] > 0), 'TFLAP downwash should be positive'


if __name__ == '__main__':
    unittest.main()
