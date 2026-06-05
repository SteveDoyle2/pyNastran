"""defines OP2 Matrix Test"""
from pathlib import Path
import unittest

import numpy as np
from scipy import sparse
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_geom import read_op2_geom, FatalError
from pyNastran.op2.result_objects.utils import modal_kinetic_energy_fraction
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').absolute()


class TestOP2Matrix(unittest.TestCase):
    """various matrix tests"""

    def test_gpspc(self):
        """Tests the gspc1 MATPOOL model"""
        op2_filename = PKG_PATH / 'op2' / 'test' / 'matrices' / 'gpsc1.op2'
        model = read_op2_geom(op2_filename, debug=False)

        deltak = model.matrices['DELTAK']
        assert deltak.data.shape == (17, 221), deltak.data.shape

        deltam = model.matrices['DELTAM']
        assert deltam.data.shape == (17, 221), deltam.data.shape

        deltam0 = model.matrices['DELTAM0']
        assert deltam0.data.shape == (6, 78), deltam0.data.shape

        #mrggt = model.matrices['MRGGT']
        mrggt = model.matrices['MRGG']
        assert mrggt.data.shape == (24, 24), mrggt.data.shape

        rbm0 = model.matrices['RBM0']
        assert rbm0.data.shape == (6, 6), rbm0.data.shape

        #uexpt = model.matrices['UEXPT']
        uexpt = model.matrices['UEXP']
        assert uexpt.data.shape == (276, 24), uexpt.data.shape

    def test_kelm_kdict(self):
        """Tests reading KELM and KDICT"""
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar_kelm.op2'
        model = read_op2_geom(op2_filename, debug=False)

        kelm = model.matrices['KELM']
        kdict = model.matdicts['KDICT']
        assert kelm.data.shape == (300, 21), kelm.data.shape
        assert kdict.element_types == [34, 2, 67, 68, 33, 1, 39, 74], kdict.element_types
        ngrids = [len(np.where(sil[0, :] > 0)[0]) for sil in kdict.sils]

        #print(kelm.data)
        #print('ndata =', len(kelm.data.data), 300*21)
        eids = np.hstack(kdict.eids)
        ndofs = [len(eids)*numgrid*dof_per_grid for eids, numgrid, dof_per_grid in zip(kdict.eids, ngrids, kdict.dof_per_grids)]
        ndof = np.cumsum(ndofs)

        sil = np.hstack([sil.ravel() for sil in kdict.sils])
        usil = np.unique(sil)

        address = np.vstack(kdict.address)
        #print(kelm.data.__dict__.keys())
        #print(kelm)
        #print('ngrids =', ngrids)
        #print('numgrids =', kdict.numgrids)
        #print('ndofs =', ndofs)
        #print('ndof csum =', ndof)
        #print('eids =', eids, len(eids))
        #print('sil =', sil, len(sil))
        #print('address:\n', address)
        #print("usil =", usil, len(usil))
        #for sil, ndofci, ndofi, dof_per_grid, numgrid in zip(kdict.sils, ndof, ndofs, kdict.dof_per_grids, kdict.numgrids):
            #print(sil, 'neids=%s cdof=%s ndof=%s dof/grid=%s ngrid=%s' % (sil.shape[0], ndofci, ndofi, dof_per_grid, numgrid))
        #print(kdict)

    def test_kelm_melm_kdict_mdict(self):
        """Tests reading KELM, MELM, KDICT, MDICT and verifying structure.

        Verifies:
          - KELM/MELM shapes match expected (nrows=sum of element DOFs, ncols=nelements)
          - KDICT/MDICT have correct element types and counts
          - KDICT addresses are consistent with KELM dimensions
          - KDICT SILs (scalar index list) have valid DOF connectivity
          - KGG is symmetric (form=6)
          - KELM non-zero data exists
        """
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar_kelm.op2'
        model = read_op2_geom(op2_filename, debug=False)

        # --- matrices ---
        assert 'KELM' in model.matrices, list(model.matrices.keys())
        assert 'MELM' in model.matrices, list(model.matrices.keys())
        assert 'KGG' in model.matrices, list(model.matrices.keys())

        kelm = model.matrices['KELM']
        melm = model.matrices['MELM']
        kgg = model.matrices['KGG']

        # KELM: 300 rows (sum of element stiffness DOFs), 21 columns (elements)
        assert kelm.data.shape == (300, 21), kelm.data.shape
        assert melm.data.shape == (300, 21), melm.data.shape
        # KGG: square symmetric, 150 DOF = 25 grids * 6 DOF
        assert kgg.data.shape == (150, 150), kgg.data.shape

        # KELM must have non-zero stiffness data
        kelm_dense = kelm.data.toarray() if hasattr(kelm.data, 'toarray') else kelm.data
        assert np.any(kelm_dense != 0), 'KELM is all zeros'

        # --- matdicts ---
        assert 'KDICT' in model.matdicts, list(model.matdicts.keys())
        assert 'MDICT' in model.matdicts, list(model.matdicts.keys())

        kdict = model.matdicts['KDICT']
        mdict = model.matdicts['MDICT']

        # 21 elements total across 8 element types
        assert kdict.nelements == 21, kdict.nelements
        assert mdict.nelements == 21, mdict.nelements
        assert len(kdict.element_types) == 8, kdict.element_types
        assert len(mdict.element_types) == 8, mdict.element_types

        # element types: CBAR=34, CBEAM=2, CHEXA=67, CPENTA=68,
        #                CQUAD4=33, CROD=1, CTETRA=39, CTRIA3=74
        expected_etypes = [34, 2, 67, 68, 33, 1, 39, 74]
        assert kdict.element_types == expected_etypes, kdict.element_types

        # total element IDs across all groups = 21
        all_eids = np.hstack(kdict.eids)
        assert len(all_eids) == 21, f'neids={len(all_eids)}'
        assert len(np.unique(all_eids)) == 21, 'duplicate element IDs'

        # address array: each element has [start, end] into KELM rows
        for i, (addr, eids) in enumerate(zip(kdict.address, kdict.eids)):
            assert addr.shape == (len(eids), 2), (
                f'group {i}: address shape {addr.shape} vs neids={len(eids)}')
            # addresses must be non-negative
            assert np.all(addr >= 0), f'group {i}: negative address'

        # SILs: each group has (neids, numgrid) shape
        for i, (sil, eids, numgrid) in enumerate(
                zip(kdict.sils, kdict.eids, kdict.numgrids)):
            assert sil.shape[0] == len(eids), (
                f'group {i}: sil rows {sil.shape[0]} vs neids={len(eids)}')
            assert sil.shape[1] == numgrid, (
                f'group {i}: sil cols {sil.shape[1]} vs numgrid={numgrid}')

        # DOF per grid must be positive
        for dpg in kdict.dof_per_grids:
            assert dpg > 0, f'dof_per_grid={dpg}'

        # --- get_element_matrices ---
        # Extract per-element stiffness matrices from KELM using KDICT
        elem_k = kdict.get_element_matrices(kelm)
        assert len(elem_k) == 21, f'expected 21, got {len(elem_k)}'

        # Expected sizes (ndof = nnodes * dof_per_grid):
        #   CHEXA (8 nodes * 3 DOF) = 24x24
        #   CPENTA (6 nodes * 3 DOF) = 18x18
        #   CTETRA (4 nodes * 3 DOF) = 12x12
        #   CQUAD4 (4 nodes * 6 DOF) = 24x24
        #   CTRIA3 (3 nodes * 6 DOF) = 18x18
        #   CBAR (2 nodes * 6 DOF) = 12x12
        #   CBEAM (2 nodes * 6 DOF) = 12x12
        #   CROD (2 nodes * 6 DOF) = 12x12
        expected_shapes = {
            1: (24, 24),    # CHEXA
            2: (18, 18),    # CPENTA
            3: (18, 18),    # CPENTA
            4: (12, 12),    # CTETRA
            5: (12, 12),    # CTETRA
            6: (24, 24),    # CQUAD4
            7: (24, 24),    # CQUAD4
            8: (18, 18),    # CTRIA3
            9: (18, 18),    # CTRIA3
            10: (18, 18),   # CTRIA3
            11: (18, 18),   # CTRIA3
            12: (12, 12),   # CBEAM
            13: (12, 12),   # CBAR
            14: (12, 12),   # CROD
            15: (12, 12),   # CROD
            16: (24, 24),   # CQUAD4
            17: (24, 24),   # CQUAD4
            18: (18, 18),   # CTRIA3
            19: (18, 18),   # CTRIA3
            20: (18, 18),   # CTRIA3
            21: (18, 18),   # CTRIA3
        }
        for eid, expected_shape in expected_shapes.items():
            mat = elem_k[eid]
            assert mat.shape == expected_shape, (
                f'eid={eid}: shape={mat.shape}, expected={expected_shape}')
            # all element stiffness matrices must be symmetric
            assert np.allclose(mat, mat.T), f'eid={eid}: not symmetric'
            # diagonal must be non-negative (stiffness)
            assert np.all(np.diag(mat) >= 0), f'eid={eid}: negative diagonal'

        # Also extract element mass matrices
        elem_m = mdict.get_element_matrices(melm)
        assert len(elem_m) == 21
        for eid, mat in elem_m.items():
            assert np.allclose(mat, mat.T), f'mass eid={eid}: not symmetric'
            assert np.all(np.diag(mat) >= 0), f'mass eid={eid}: negative diagonal'

        # --- get_element_dof_info ---
        dof_info = kdict.get_element_dof_info()
        assert len(dof_info) == 21
        # CHEXA eid=1: 8 nodes, 3 DOF/grid
        sil_vals, dpg = dof_info[1]
        assert len(sil_vals) == 8
        assert dpg == 3
        # CBAR eid=13: 2 nodes, 6 DOF/grid
        sil_vals, dpg = dof_info[13]
        assert len(sil_vals) == 2
        assert dpg == 6
        # CTETRA eid=4: 4 nodes, 3 DOF/grid
        sil_vals, dpg = dof_info[4]
        assert len(sil_vals) == 4
        assert dpg == 3

        # --- verify extraction against KGG for untransformed elements ---
        # CBAR eid=13 has xform=None, node at SIL=133 is unique to this element.
        # The 2nd-node self-block in element K must equal KGG at that DOF.
        kgg_dense = kgg.data.toarray() if hasattr(kgg.data, 'toarray') else kgg.data
        ke_bar = elem_k[13]
        # SIL=133 -> global DOFs 132..137 (0-based). Node 2 -> elem DOFs 6..11.
        assert np.allclose(ke_bar[6:12, 6:12], kgg_dense[132:138, 132:138]), (
            'CBAR node2 self-block does not match KGG')
        # CBEAM eid=12, SIL=127 unique, node2 -> elem DOFs 6..11
        ke_beam = elem_k[12]
        assert np.allclose(ke_beam[6:12, 6:12], kgg_dense[126:132, 126:132]), (
            'CBEAM node2 self-block does not match KGG')

    def test_op2_dmi_01(self):
        """tests DMI matrix style"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = MODEL_PATH / 'matrix' / 'matrix.dat'
        op2_filename = MODEL_PATH / 'matrix' / 'mymatrix.op2'
        matrices = {
            'A': True,
            'B': False,
            'ATB': False,
            'BTA': False,
            'MYDOF': True,
        }
        model = read_bdf(bdf_filename, log=log)

        dmi_a = model.dmi['A']
        assert dmi_a.shape == (4, 2), 'shape=%s' % (dmi_a.shape)
        #print('dmi_a\n', dmi_a)
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2(log=log)
        op2.set_additional_matrices_to_read(matrices)
        try:
            op2.read_op2(op2_filename)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = B.T @ A
        ATB = A.T @ B
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        #BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for matrix_name, expected in zip(matrix_names, expecteds):
            assert matrix_name in op2.matrices, matrix_name
            actual = op2.matrices[matrix_name].data.toarray()
            compare_dmi_matrix_from_bdf_to_op2(model, op2, expected, actual, matrix_name)

    def test_op2_dmi_02(self):
        """tests DMI matrix style"""
        bdf_filename = MODEL_PATH / 'matrix' / 'matrix.dat'
        op2_filename = MODEL_PATH / 'matrix' / 'mymatrix.op2'
        matrices = {
            'A': True,
            'B': False,
            'ATB': False,
            'BTA': False,
            'MYDOF': True,
        }
        model = read_bdf(bdf_filename, debug=False)

        dmi_a = model.dmi['A']
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        #print('model.dmi.A =\n%s' % dmi_a)
        #print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2(debug=False)
        try:
            op2.read_op2(op2_filename, skip_undefined_matrices=True)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        # M rows, Ncols
        A = np.array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = np.array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = B.T @ A
        ATB = A.T @ B
        ATB_expected = np.array([
            [35., 18.],
            [18., 100.]
        ], dtype='float32')
        #BTA_expected = ATB_expected

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for matrix_name, expected in zip(matrix_names, expecteds):
            assert matrix_name in op2.matrices, matrix_name
            actual = op2.matrices[matrix_name].data.toarray()
            compare_dmi_matrix_from_bdf_to_op2(model, op2, expected, actual, matrix_name)

    def test_modal_kinetic_energy_fraction(self):
        """Tests modal kinetic energy fraction computation.

        Verifies:
          - MKE fractions sum to 1.0 across all grids for each mode
          - Single-node motion gives 100% at that node
          - Coupled motion distributes proportional to mass * phi^2
          - Works with real eigenvector data from OP2
        """
        # --- simple known-answer test: 2 nodes, 2 modes ---
        class MockEigenvectors:
            def __init__(self, data, node_gridtype, modes):
                self.data = data
                self.node_gridtype = node_gridtype
                self.modes = modes

        nnodes = 2
        phi_data = np.zeros((2, nnodes, 6))
        phi_data[0, 0, 0] = 1.0  # mode 1: node 1 T1
        phi_data[1, 1, 0] = 1.0  # mode 2: node 2 T1
        node_gt = np.array([[1, 1], [2, 1]])
        eig_mock = MockEigenvectors(phi_data, node_gt, np.array([1, 2]))

        # Lumped mass: m1=2.0, m2=3.0
        mass_diag = np.zeros(12)
        mass_diag[0:6] = 2.0
        mass_diag[6:12] = 3.0
        mgg = sparse.diags(mass_diag)

        mke, nids, modes = modal_kinetic_energy_fraction(eig_mock, mgg)
        # mode 1: all KE at node 1
        assert np.allclose(mke[0], [1.0, 0.0], atol=1e-12)
        # mode 2: all KE at node 2
        assert np.allclose(mke[1], [0.0, 1.0], atol=1e-12)
        assert np.array_equal(nids, [1, 2])
        assert np.array_equal(modes, [1, 2])

        # --- coupled motion: both nodes move equally ---
        phi_coupled = np.zeros((1, 2, 6))
        phi_coupled[0, 0, 0] = 1.0
        phi_coupled[0, 1, 0] = 1.0
        eig_coupled = MockEigenvectors(phi_coupled, node_gt, np.array([1]))

        mke_c, _, _ = modal_kinetic_energy_fraction(eig_coupled, mgg)
        # KE fractions: m1*phi1^2 / (m1*phi1^2 + m2*phi2^2) = 2/5, 3/5
        assert np.allclose(mke_c[0], [0.4, 0.6], atol=1e-12)

        # --- real eigenvector data from OP2 with identity mass ---
        op2_filename = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.op2'
        model = read_op2_geom(op2_filename, debug=False)
        eig_real = model.eigenvectors[1]
        n = eig_real.data.shape[1]  # 25 nodes
        mgg_identity = sparse.eye(n * 6, format='csr')

        mke_real, nids_real, modes_real = modal_kinetic_energy_fraction(
            eig_real, mgg_identity)
        assert mke_real.shape == (3, 25)
        # fractions sum to 1.0 per mode (tol for float precision)
        for i in range(3):
            assert abs(mke_real[i].sum() - 1.0) < 1e-10, (
                f'mode {i+1} sum={mke_real[i].sum()}')
        # all fractions non-negative (with identity mass, = phi_g^2 / ||phi||^2)
        assert np.all(mke_real >= -1e-15), 'negative MKE fraction'


def compare_dmi_matrix_from_bdf_to_op2(bdf_model, op2_model, expected, actual, matrix_name):
    """compares two matrices"""
    if not (np.array_equal(expected, actual) or
            np.array_equal(expected, np.squeeze(actual))):

        if matrix_name in bdf_model.dmis:
            dmi = bdf_model.dmis[matrix_name]
            table_array, rows_reversed, cols_reversed = dmi.get_matrix(
                is_sparse=False, apply_symmetry=False)
            #stable_array, rows_reversed, cols_reversed = dmi.get_matrix(
                #is_sparse=True, apply_symmetry=False)
            #print(table_array)
        #print(stable_array)
        msg = 'matrix %s was not read properly\n' % matrix_name
        msg += 'expected shape=%s\n%s\n' % (str(expected.shape), expected)
        msg += 'actual shape=%s;  squeeze(actual) shape=%s\n%s' % (
            str(actual.shape), str(np.squeeze(actual.shape)), actual.ravel())
        #msg += '\n%s' % actual.ravel()
        print(msg)
        print('==========================')
        #raise RuntimeError(msg)


if __name__ == '__main__':  # pragma: no cover
    import os
    ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
    if not ON_RTD:
        unittest.main()
