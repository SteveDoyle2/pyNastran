"""Tests for GUI utility functions and performance fixes.

Tests the correctness of:
1. map_centroidal_result (vectorized scatter-add via np.add.at)
2. _convert_ids_to_vtk_idtypearray (numpy bulk VTK array construction)
3. eid_to_nid_map (verifies map content is correct)
4. _vectorized_shell_normals (vectorized np.cross vs per-element Normal())
5. build_offset_normals_dims (vectorized vs fallback loop, mixed element types)
"""
import unittest
import numpy as np
from pyNastran.converters.nastran.gui.utils import (
    _build_map_centroidal_result, _vectorized_shell_normals,
    build_offset_normals_dims)
from pyNastran.bdf.bdf import BDF


class TestMapCentroidalResult(unittest.TestCase):
    """Tests map_centroidal_result produces correct nodal averages from centroidal data."""

    def _make_model_with_map(self):
        """Create a 2-quad model and build the centroidal mapper.

        Grid:  1---2---3
               |   |   |
               4---5---6

        Elements: eid=1 (1,2,5,4), eid=2 (2,3,6,5)
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 1., 0.])
        model.add_grid(2, [1., 1., 0.])
        model.add_grid(3, [2., 1., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        model.add_grid(6, [2., 0., 0.])
        model.add_cquad4(1, 1, [1, 2, 5, 4])
        model.add_cquad4(2, 1, [2, 3, 6, 5])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()

        nid_map = {}
        for i, nid in enumerate(sorted(model.nodes.keys())):
            nid_map[nid] = i

        _build_map_centroidal_result(model, nid_map)
        return model, nid_map

    def test_uniform_centroidal_data(self):
        """If all elements have the same centroidal value, all nodes get that value."""
        model, nid_map = self._make_model_with_map()
        centroidal_data = np.array([10.0, 10.0], dtype='float32')
        result = model.map_centroidal_result(centroidal_data)
        np.testing.assert_allclose(result[list(nid_map.values())], 10.0, rtol=1e-6)

    def test_varying_centroidal_data(self):
        """Shared nodes get the average of adjacent elements.

        Element 1 value=2.0, element 2 value=4.0.
        Nodes 2,5 shared -> average=3.0; nodes 1,4 -> 2.0; nodes 3,6 -> 4.0
        """
        model, nid_map = self._make_model_with_map()
        centroidal_data = np.array([2.0, 4.0], dtype='float32')
        result = model.map_centroidal_result(centroidal_data)

        assert abs(result[nid_map[1]] - 2.0) < 1e-6
        assert abs(result[nid_map[4]] - 2.0) < 1e-6
        assert abs(result[nid_map[3]] - 4.0) < 1e-6
        assert abs(result[nid_map[6]] - 4.0) < 1e-6
        assert abs(result[nid_map[2]] - 3.0) < 1e-6
        assert abs(result[nid_map[5]] - 3.0) < 1e-6

    def test_single_element(self):
        """Single CTRIA3: all nodes get the element centroidal value."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0.5, 1., 0.])
        model.add_ctria3(10, 1, [1, 2, 3])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()
        nid_map = {1: 0, 2: 1, 3: 2}
        _build_map_centroidal_result(model, nid_map)

        centroidal_data = np.array([7.5], dtype='float32')
        result = model.map_centroidal_result(centroidal_data)
        np.testing.assert_allclose(result[:3], 7.5, rtol=1e-6)


class TestConvertIdsToVtkIdTypeArray(unittest.TestCase):
    """Tests that _convert_ids_to_vtk_idtypearray produces correct VTK arrays."""

    def test_numpy_array_input(self):
        """Standard numpy int array input."""
        from pyNastran.gui.utils.vtk.vtk_utils import _convert_ids_to_vtk_idtypearray
        ids = np.array([0, 5, 10, 99], dtype='int64')
        result = _convert_ids_to_vtk_idtypearray(ids)
        assert result.GetNumberOfTuples() == 4
        assert result.GetValue(0) == 0
        assert result.GetValue(1) == 5
        assert result.GetValue(2) == 10
        assert result.GetValue(3) == 99

    def test_single_int_input(self):
        """Single integer input (not array)."""
        from pyNastran.gui.utils.vtk.vtk_utils import _convert_ids_to_vtk_idtypearray
        result = _convert_ids_to_vtk_idtypearray(42)
        assert result.GetNumberOfTuples() == 1
        assert result.GetValue(0) == 42

    def test_empty_array(self):
        """Empty array produces empty VTK array."""
        from pyNastran.gui.utils.vtk.vtk_utils import _convert_ids_to_vtk_idtypearray
        ids = np.array([], dtype='int64')
        result = _convert_ids_to_vtk_idtypearray(ids)
        assert result.GetNumberOfTuples() == 0

    def test_large_ids(self):
        """Large IDs (>2^31) work correctly for 64-bit VTK builds."""
        from pyNastran.gui.utils.vtk.vtk_utils import _convert_ids_to_vtk_idtypearray
        ids = np.array([0, 1, 2**31 - 1], dtype='int64')
        result = _convert_ids_to_vtk_idtypearray(ids)
        assert result.GetNumberOfTuples() == 3
        assert result.GetValue(2) == 2**31 - 1


class TestEidToNidMap(unittest.TestCase):
    """Tests that the eid_to_nid_map is correctly built from model elements."""

    def _make_model(self):
        """Simple model: 2 CONRODs + 1 CTRIA3."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0.5, 1., 0.])
        model.add_grid(4, [2., 0., 0.])
        model.add_conrod(10, 1, [1, 2], A=1.0)
        model.add_conrod(20, 1, [2, 4], A=1.0)
        model.add_ctria3(30, 1, [1, 2, 3])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()
        return model

    def test_eid_to_nid_map_content(self):
        """Verify the eid->nodes mapping is correctly extracted."""
        model = self._make_model()
        eid_to_nid_map = {eid: elem.nodes for eid, elem in model.elements.items()}

        assert len(eid_to_nid_map) == 3
        nids_10 = [n.nid if hasattr(n, 'nid') else n for n in eid_to_nid_map[10]]
        assert nids_10 == [1, 2]
        nids_20 = [n.nid if hasattr(n, 'nid') else n for n in eid_to_nid_map[20]]
        assert nids_20 == [2, 4]
        nids_30 = [n.nid if hasattr(n, 'nid') else n for n in eid_to_nid_map[30]]
        assert nids_30 == [1, 2, 3]

    def test_eid_to_nid_map_is_invariant(self):
        """Building the map multiple times yields the same result (no mutation)."""
        model = self._make_model()
        map1 = {eid: elem.nodes for eid, elem in model.elements.items()}
        map2 = {eid: elem.nodes for eid, elem in model.elements.items()}
        assert map1.keys() == map2.keys()
        for eid in map1:
            nids1 = [n.nid if hasattr(n, 'nid') else n for n in map1[eid]]
            nids2 = [n.nid if hasattr(n, 'nid') else n for n in map2[eid]]
            assert nids1 == nids2


class TestVectorizedShellNormals(unittest.TestCase):
    """Tests _vectorized_shell_normals matches per-element Normal() results.

    Tolerances: atol=1e-5 (float32 storage precision).
    Interactions tested: tri vs quad dispatch, tilted geometry, mixed element types.
    """

    def _make_quad_model(self):
        """Two CQUAD4s in XY plane.

        Grid:  1---2---3
               |   |   |
               4---5---6
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 1., 0.])
        model.add_grid(2, [1., 1., 0.])
        model.add_grid(3, [2., 1., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        model.add_grid(6, [2., 0., 0.])
        model.add_cquad4(1, 1, [1, 2, 5, 4])
        model.add_cquad4(2, 1, [2, 3, 6, 5])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()
        return model

    def _make_tri_model(self):
        """Two CTRIA3s in XY plane."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])
        model.add_ctria3(2, 1, [1, 3, 4])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()
        return model

    def test_quad_normals_xy_plane(self):
        """CQUAD4s in XY plane -> normals along Z axis."""
        model = self._make_quad_model()
        nelements = 2
        eid_map = {1: 0, 2: 1}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}
        xyz_cid0 = np.array([
            [0., 1., 0.], [1., 1., 0.], [2., 1., 0.],
            [0., 0., 0.], [1., 0., 0.], [2., 0., 0.],
        ], dtype='float64')
        normals = np.full((nelements, 3), np.nan, dtype='float32')

        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

        for i in range(nelements):
            assert not np.isnan(normals[i, 0]), f'element {i} normal is NaN'
            np.testing.assert_allclose(abs(normals[i, 2]), 1.0, atol=1e-6)
            np.testing.assert_allclose(normals[i, :2], [0., 0.], atol=1e-6)

    def test_tri_normals_xy_plane(self):
        """CTRIA3s in XY plane -> normals along Z axis."""
        model = self._make_tri_model()
        nelements = 2
        eid_map = {1: 0, 2: 1}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3}
        xyz_cid0 = np.array([
            [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
        ], dtype='float64')
        normals = np.full((nelements, 3), np.nan, dtype='float32')

        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

        for i in range(nelements):
            assert not np.isnan(normals[i, 0])
            np.testing.assert_allclose(abs(normals[i, 2]), 1.0, atol=1e-6)
            np.testing.assert_allclose(normals[i, :2], [0., 0.], atol=1e-6)

    def test_vectorized_matches_per_element(self):
        """Vectorized normals match element.Normal() for each CQUAD4."""
        model = self._make_quad_model()
        nelements = 2
        eid_map = {1: 0, 2: 1}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}
        xyz_cid0 = np.array([
            [0., 1., 0.], [1., 1., 0.], [2., 1., 0.],
            [0., 0., 0.], [1., 0., 0.], [2., 0., 0.],
        ], dtype='float64')
        normals = np.full((nelements, 3), np.nan, dtype='float32')

        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

        for eid, elem in model.elements.items():
            expected = elem.Normal()
            ieid = eid_map[eid]
            np.testing.assert_allclose(
                normals[ieid], expected, atol=1e-5,
                err_msg=f'eid={eid} normal mismatch')

    def test_tilted_quad(self):
        """A quad tilted 45deg has correct normal direction."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 0.707, 0.707])
        model.add_grid(4, [0., 0.707, 0.707])
        model.add_cquad4(1, 1, [1, 2, 3, 4])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()

        eid_map = {1: 0}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3}
        xyz_cid0 = np.array([
            [0., 0., 0.], [1., 0., 0.],
            [1., 0.707, 0.707], [0., 0.707, 0.707],
        ], dtype='float64')
        normals = np.full((1, 3), np.nan, dtype='float32')

        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

        expected = model.elements[1].Normal()
        np.testing.assert_allclose(normals[0], expected, atol=1e-5)

    def test_mixed_tri_quad(self):
        """Model with both CTRIA3 and CQUAD4."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [0.5, 0.5, 0.])
        model.add_ctria3(1, 1, [1, 2, 5])
        model.add_cquad4(2, 1, [2, 3, 4, 1])
        model.add_pshell(1, mid1=1, t=0.1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()

        nelements = 2
        eid_map = {1: 0, 2: 1}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4}
        xyz_cid0 = np.array([
            [0., 0., 0.], [1., 0., 0.], [1., 1., 0.],
            [0., 1., 0.], [0.5, 0.5, 0.],
        ], dtype='float64')
        normals = np.full((nelements, 3), np.nan, dtype='float32')

        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

        for eid, elem in model.elements.items():
            expected = elem.Normal()
            ieid = eid_map[eid]
            np.testing.assert_allclose(
                normals[ieid], expected, atol=1e-5,
                err_msg=f'eid={eid} normal mismatch')


class TestBuildOffsetNormalsDims(unittest.TestCase):
    """Tests build_offset_normals_dims vectorized path matches fallback loop.

    Tolerances: atol=1e-5 (float32 precision).
    Interactions tested: mixed element types (shell, solid, rod, spring),
    vectorized vs non-vectorized paths produce identical results.
    """

    def _make_mixed_model(self):
        """Model with CQUAD4, CTRIA3, CONROD, CHEXA, CELAS2."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])
        model.add_grid(9, [2., 0., 0.])

        model.add_cquad4(1, 1, [1, 2, 3, 4])
        model.add_ctria3(2, 1, [1, 2, 5])
        model.add_conrod(3, 1, [1, 9], A=1.0)
        model.add_chexa(4, 3, [1, 2, 3, 4, 5, 6, 7, 8])
        model.add_celas2(5, 100.0, [1, 2])

        model.add_pshell(1, mid1=1, t=0.1)
        model.add_psolid(3, 1)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()
        return model

    def _get_maps(self, model):
        eid_map = {eid: i for i, eid in enumerate(sorted(model.elements.keys()))}
        nid_map = {nid: i for i, nid in enumerate(sorted(model.nodes.keys()))}
        nnodes = len(model.nodes)
        xyz_cid0 = np.zeros((nnodes, 3), dtype='float64')
        for nid_val, node in model.nodes.items():
            xyz_cid0[nid_map[nid_val]] = node.get_position()
        return eid_map, nid_map, xyz_cid0

    def test_vectorized_matches_fallback(self):
        """Vectorized path produces identical results to fallback loop."""
        model = self._make_mixed_model()
        nelements = len(model.elements)
        eid_map, nid_map, xyz_cid0 = self._get_maps(model)

        out_old = build_offset_normals_dims(model, eid_map, nelements)
        out_new = build_offset_normals_dims(model, eid_map, nelements,
                                            xyz_cid0=xyz_cid0, nid_map=nid_map)

        names = ['normals', 'offset', 'xoffset', 'yoffset', 'zoffset',
                 'element_dim', 'nnodes_array']
        for name, a, b in zip(names, out_old, out_new):
            nan_a = np.isnan(a) if a.dtype.kind == 'f' else np.zeros_like(a, dtype=bool)
            nan_b = np.isnan(b) if b.dtype.kind == 'f' else np.zeros_like(b, dtype=bool)
            np.testing.assert_array_equal(nan_a, nan_b, err_msg=f'{name}: NaN mismatch')
            mask = ~nan_a
            if mask.any():
                np.testing.assert_allclose(
                    a[mask], b[mask], atol=1e-5, err_msg=f'{name}: value mismatch')

    def test_element_dim_values(self):
        """Element dimensions are correct for each element type."""
        model = self._make_mixed_model()
        nelements = len(model.elements)
        eid_map, nid_map, xyz_cid0 = self._get_maps(model)

        out = build_offset_normals_dims(model, eid_map, nelements,
                                        xyz_cid0=xyz_cid0, nid_map=nid_map)
        element_dim = out[5]
        nnodes_arr = out[6]

        assert element_dim[eid_map[1]] == 2  # CQUAD4 -> 2D
        assert element_dim[eid_map[2]] == 2  # CTRIA3 -> 2D
        assert element_dim[eid_map[3]] == 1  # CONROD -> 1D
        assert element_dim[eid_map[4]] == 3  # CHEXA -> 3D
        assert element_dim[eid_map[5]] == 0  # CELAS2 -> 0D

        assert nnodes_arr[eid_map[1]] == 4  # CQUAD4
        assert nnodes_arr[eid_map[2]] == 3  # CTRIA3
        assert nnodes_arr[eid_map[3]] == 2  # CONROD
        assert nnodes_arr[eid_map[4]] == 8  # CHEXA
        assert nnodes_arr[eid_map[5]] == 2  # CELAS2

    def test_shell_offset_with_pshell(self):
        """Shell elements with PSHELL get correct z0 offset from prop.z1."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])
        model.add_pshell(1, mid1=1, t=0.5)
        model.add_mat1(1, 1e7, None, 0.3)
        model.cross_reference()

        eid_map = {1: 0}
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3}
        xyz_cid0 = np.array([
            [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]
        ], dtype='float64')

        out = build_offset_normals_dims(model, eid_map, 1,
                                        xyz_cid0=xyz_cid0, nid_map=nid_map)
        offset_val = out[1][0]
        # PSHELL z1 defaults to -t/2 = -0.25
        expected_z0 = model.properties[1].z1
        np.testing.assert_allclose(offset_val, expected_z0, atol=1e-6)


if __name__ == '__main__':
    unittest.main()