"""tests for structural_lcp.py (structural loads -> load collection points)"""
from pathlib import Path
import tempfile
import unittest
import numpy as np

from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.dev.tools.pressure_map.structural_lcp import (
    structural_lcp_map, structural_loads_to_lcp_model, get_nodal_loads,
    consistent_pressure_nodal_forces, _sorted_xyz)


def _capture_log(level: str='info') -> tuple[SimpleLogger, list[str]]:
    """info-level logger whose messages are appended to the returned list"""
    msgs = []
    log = SimpleLogger(level=level,
                       log_func=lambda typ, filename, lineno, msg: msgs.append(msg))
    return log, msgs


def _load_model(log: SimpleLogger) -> BDF:
    """
    fuse: CTRIA3 eid=1 (pid=10); nodes 1, 2, 3  (flat, z=0)
    wing: CQUAD4 eid=2 (pid=20); nodes 2, 4, 5, 6 trapezoid (z=0)
    node 2 (1,0,0) is shared (junction)

    sid=1: PLOAD2 p=2 on both elements
    sid=2: FORCE (0,0,10) at node 5 (wing only), FORCE (0,0,-4) at node 2 (junction)
           MOMENT (1,0,0) at node 6
    """
    model = BDF(log=log)
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(2, [1., 0., 0.])
    model.add_grid(3, [0., 1., 0.])
    # trapezoid: base 4 at y=0 (x=1..5), top 2 at y=2 (x=2..4)
    model.add_grid(4, [5., 0., 0.])
    model.add_grid(5, [4., 2., 0.])
    model.add_grid(6, [2., 2., 0.])
    model.add_mat1(1, 3.0e7, None, 0.3)
    model.add_pshell(10, mid1=1, t=0.1)
    model.add_pshell(20, mid1=1, t=0.1)
    model.add_ctria3(1, 10, [1, 2, 3])
    model.add_cquad4(2, 20, [2, 4, 5, 6])
    model.add_pload2(1, 2.0, [1, 2])
    model.add_force(2, 5, 10., [0., 0., 1.])
    model.add_force(2, 2, 4., [0., 0., -1.])
    model.add_moment(2, 6, 1., [1., 0., 0.])
    model.add_load(3, 1.0, [1.0, 1.0], [1, 2])
    model.cross_reference()
    return model


def _lcp_model(log: SimpleLogger, xyzs: dict[int, list[float]]) -> BDF:
    model = BDF(log=log)
    model.bdf_filename = 'lcp_nodes.bdf'
    for nid, xyz in xyzs.items():
        model.add_grid(nid, xyz)
    model.cross_reference()
    return model


def _cards(model_out: BDF) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    forces = {}
    moments = {}
    for load in model_out.loads.get(2, []):
        forces[load.node] = forces.get(load.node, np.zeros(3)) + load.scaled_vector
    for load in model_out.loads.get(3, []):
        moments[load.node] = moments.get(load.node, np.zeros(3)) + load.scaled_vector
    return forces, moments


def _total_about(lcp_model: BDF, forces, moments, point) -> tuple[np.ndarray, np.ndarray]:
    total_f = np.zeros(3)
    total_m = np.zeros(3)
    for nid, force in forces.items():
        xyz = lcp_model.nodes[nid].get_position()
        total_f += force
        total_m += moments[nid] + np.cross(xyz - point, force)
    return total_f, total_m


class TestConsistentPressure(unittest.TestCase):
    def test_tri_linear_pressure(self):
        """F_a = (A n)/12 (2 p_a + p_b + p_c)"""
        xyz = np.array([[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]]])
        force = consistent_pressure_nodal_forces(
            xyz, np.array([[1., 2., 3.]]), np.full((1, 3), np.nan))
        np.testing.assert_allclose(force[0, :, 2], [7/6, 8/6, 9/6], atol=1e-14)
        np.testing.assert_allclose(force[0, :, :2], 0., atol=1e-14)

    def test_trapezoid_uniform_pressure_moment(self):
        """
        non-parallelogram quad: nodal forces aren't p*A/4, but they
        reproduce F = p A n and M = c x F about the area centroid exactly
        """
        xyz = np.array([[[1., 0., 0.], [5., 0., 0.], [4., 2., 0.], [2., 2., 0.]]])
        p = 2.
        force = consistent_pressure_nodal_forces(
            xyz, np.full((1, 4), p), np.full((1, 3), np.nan))[0]
        area = 6.
        centroid = np.array([3., 8/9, 0.])
        total_f = force.sum(axis=0)
        total_m = np.cross(xyz[0], force).sum(axis=0)
        np.testing.assert_allclose(total_f, [0., 0., p * area], atol=1e-13)
        np.testing.assert_allclose(total_m, np.cross(centroid, [0., 0., p * area]), atol=1e-13)
        self.assertFalse(np.allclose(force[:, 2], p * area / 4))

    def test_direction(self):
        """PLOAD4 N: load is p*A along N (not the normal)"""
        xyz = np.array([[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]]])
        force = consistent_pressure_nodal_forces(
            xyz, np.full((1, 3), 3.), np.array([[1., 0., 0.]]))
        np.testing.assert_allclose(force[0], [[2., 0., 0.]] * 3, atol=1e-14)


class TestNodalLoads(unittest.TestCase):
    def test_force_rotated_cid_and_load_combination(self):
        """FORCE in a rotated CORD2R; LOAD scale factors multiply"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        # local x = global y
        model.add_cord2r(1, [0., 0., 0.], [0., 0., 1.], [0., 1., 0.])
        model.add_force(1, 1, 5., [1., 0., 0.], cid=1)
        model.add_load(100, 2.0, [1.5], [1])
        model.cross_reference()
        nids, xyz = _sorted_xyz(model, '')
        nodal = get_nodal_loads(model, 100, nids, xyz)
        np.testing.assert_allclose(nodal['point_force'], [[0., 15., 0.]], atol=1e-12)
        np.testing.assert_allclose(nodal['point_moment'], [[0., 0., 0.]], atol=1e-12)

    def test_pload4_corner_pressures(self):
        """PLOAD4 P1-P3 on a CTRIA3 follow the element node order"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [0., 2., 0.])
        model.add_mat1(1, 3.0e7, None, 0.3)
        model.add_pshell(10, mid1=1, t=0.1)
        model.add_ctria3(1, 10, [3, 1, 2])
        model.add_pload4(1, [1], [1., 2., 3., 1.])
        model.cross_reference()
        nids, xyz = _sorted_xyz(model, '')
        nodal = get_nodal_loads(model, 1, nids, xyz)
        self.assertEqual(nodal['pressure_nid'].tolist(), [3, 1, 2])
        np.testing.assert_allclose(nodal['pressure_force'][:, 2], [7/6, 8/6, 9/6], atol=1e-14)


class TestStructuralLCP(unittest.TestCase):
    def test_groups_by_pid(self):
        """
        fuse (pid 10) -> nid 1001; wing (pid 20) -> 2001/2002.
        Junction FORCE at node 2 goes to the first group (fuse) with a
        warning.  Every balance passes and the cards reproduce the load case.
        """
        log, msgs = _capture_log()
        load_model = _load_model(log)
        lcp_model = _lcp_model(log, {1001: [0., 0., 1.], 2001: [2., 0., 1.], 2002: [4., 2., 1.]})
        groups = {
            'fuse': {'property_ids': [10], 'nodes': [1001]},
            'wing': {'property_ids': [20], 'nodes': [2001, 2002]},
        }
        model_out = structural_loads_to_lcp_model(
            load_model, lcp_model, 3, groups, reference_point=np.array([1., 2., 3.]))
        text = '\n'.join(msgs)
        self.assertEqual(text.count('PASS:'), 3, text)
        self.assertNotIn('FAIL', text)
        self.assertIn('FORCE/MOMENT grids=[2] are used by more than one group', text)

        forces, moments = _cards(model_out)
        self.assertEqual(sorted(forces), [1001, 2001, 2002])
        # fuse: tri p=2, A=1/2 -> Fz=1; junction FORCE -4
        np.testing.assert_allclose(forces[1001], [0., 0., -3.], atol=1e-12)
        # wing: quad Fz = 2*6 = 12, FORCE +10 at node 5 -> 2002
        np.testing.assert_allclose(forces[2001] + forces[2002], [0., 0., 22.], atol=1e-12)

        # load case about the origin
        tri_f = np.array([0., 0., 1.])
        quad_f = np.array([0., 0., 12.])
        case_f = tri_f + quad_f + [0., 0., 10.] + [0., 0., -4.]
        case_m = (np.cross([1/3, 1/3, 0.], tri_f) + np.cross([3., 8/9, 0.], quad_f) +
                  np.cross([4., 2., 0.], [0., 0., 10.]) + np.cross([1., 0., 0.], [0., 0., -4.]) +
                  np.array([1., 0., 0.]))
        total_f, total_m = _total_about(lcp_model, forces, moments, np.zeros(3))
        np.testing.assert_allclose(total_f, case_f, atol=1e-12)
        np.testing.assert_allclose(total_m, case_m, atol=1e-12)

    def test_element_ids_unassigned(self):
        """an element in no group is dropped: group PASS, total FAIL + warning"""
        log, msgs = _capture_log()
        load_model = _load_model(log)
        lcp_model = _lcp_model(log, {2001: [2., 0., 1.]})
        groups = {'wing': {'element_ids': [2], 'nodes': [2001]}}
        model_out = structural_loads_to_lcp_model(load_model, lcp_model, 1, groups)
        text = '\n'.join(msgs)
        self.assertIn('PASS:', text)
        self.assertIn('FAIL:', text)
        self.assertIn('loaded element ids=[1] are not in any group', text)
        forces, unused_moments = _cards(model_out)
        np.testing.assert_allclose(forces[2001], [0., 0., 12.], atol=1e-12)

    def test_bad_group_key(self):
        log = SimpleLogger(level='error')
        load_model = _load_model(log)
        lcp_model = _lcp_model(log, {2001: [2., 0., 1.]})
        with self.assertRaises(ValueError):
            structural_loads_to_lcp_model(
                load_model, lcp_model, 1, {'wing': {'property_ids': [20], 'nids': [2001]}})

    def test_total_and_files(self):
        """no groups -> 'total'; writes the loads, main deck, and csv"""
        log = SimpleLogger(level='warning')
        load_model = _load_model(log)
        with tempfile.TemporaryDirectory() as dirname:
            dirname = Path(dirname)
            lcp_filename = dirname / 'lcp_nodes.bdf'
            _lcp_model(log, {1: [0., 0., 0.], 2: [5., 2., 0.]}).write_bdf(lcp_filename)
            lcp_load_filename = dirname / 'lcp_loads.bdf'
            model_out = structural_lcp_map(
                load_model, str(lcp_filename), 3, None,
                lcp_load_filename=str(lcp_load_filename), log=log)
            self.assertTrue(lcp_load_filename.exists())
            self.assertTrue((dirname / 'lcp_loads.force_moment.csv').exists())
            main = read_bdf(dirname / 'main_lcp_loads.bdf', log=log)
            self.assertEqual(sorted(main.loads), [2, 3])
            forces, _ = _cards(model_out)
            total = sum(forces.values())
            np.testing.assert_allclose(total, [0., 0., 19.], atol=1e-12)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
