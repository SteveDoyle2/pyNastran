import os
import unittest
from numpy import unique

import pyNastran
test_path = pyNastran.__path__[0]

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2, FatalError
from pyNastran.op2.test.test_op2 import run_op2
from pyNastran.bdf.test.bdf_unit_tests import Tester
from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealPlateBilinearForce, RealPlateForceArray, RealPlateForce

class TestOP2(Tester):
    def _spike(self):
        op2 = OP2()
        op2.set_results('solidStress.oxx')
        op2.read_op2(op2_filename, vectorized=False)

    def test_set_results(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        op2_filename = os.path.join(folder, 'solid_bending', 'solid_bending.op2')
        op2 = OP2(debug=False)
        op2.set_results('stress')
        op2.read_op2(op2_filename, vectorized=False)
        self.assertEqual(len(op2.cpenta_stress), 0), len(op2.cpenta_stress)
        self.assertEqual(len(op2.chexa_stress), 0), len(op2.chexa_stress)
        self.assertEqual(len(op2.ctetra_stress), 1), len(op2.ctetra_stress)
        self.assertEqual(len(op2.displacements), 0), len(op2.displacements)

        op2 = OP2(debug=False)
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2_filename, vectorized=False)
        self.assertEqual(len(op2.cpenta_stress), 0), len(op2.cpenta_stress)
        self.assertEqual(len(op2.chexa_stress), 0), len(op2.chexa_stress)
        self.assertEqual(len(op2.ctetra_stress), 1), len(op2.ctetra_stress)
        self.assertEqual(len(op2.displacements), 1), len(op2.displacements)

        op2 = OP2(debug=False)
        op2.set_results('stress')
        op2.read_op2(op2_filename, vectorized=True)
        self.assertEqual(len(op2.cpenta_stress), 0), len(op2.cpenta_stress)
        self.assertEqual(len(op2.chexa_stress), 0), len(op2.chexa_stress)
        self.assertEqual(len(op2.ctetra_stress), 1), len(op2.ctetra_stress)
        self.assertEqual(len(op2.displacements), 0), len(op2.displacements)

        op2 = OP2(debug=False)
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2_filename, vectorized=True)
        self.assertEqual(len(op2.cpenta_stress), 0), len(op2.cpenta_stress)
        self.assertEqual(len(op2.chexa_stress), 0), len(op2.chexa_stress)
        self.assertEqual(len(op2.ctetra_stress), 1), len(op2.ctetra_stress)
        self.assertEqual(len(op2.displacements), 1), len(op2.displacements)

    def test_op2_01(self):
        op2_filename = os.path.join('solid_bending', 'solid_bending.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2_filename)
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, is_vector=False,
                debug=debug, stopOnFailure=True, binary_debug=True)
        assert os.path.exists('debug.out'), os.listdir('.')
        os.remove('debug.out')

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, is_vector=True,
                debug=debug, stopOnFailure=True, binary_debug=True)
        assert os.path.exists('debug.out'), os.listdir('.')
        os.remove('debug.out')

    def test_op2_02(self):
        op2_filename = os.path.join('plate_py', 'plate_py.op2')
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        op2file = os.path.join(folder, op2_filename)
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, is_vector=False,
                debug=debug, stopOnFailure=True)

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, iSubcases=[],
                write_f06=write_f06, is_vector=True,
                debug=debug, stopOnFailure=True)

    def test_op2_eids_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'sol_101_elements', 'static_solid_comp_bar.bdf')
        op2_filename = os.path.join(folder, 'sol_101_elements', 'static_solid_comp_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        op2file = os.path.join(folder, op2_filename)
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename)

        debug = False
        debug_file = 'debug.out'
        op2v = OP2(debug=debug, debug_file=debug_file)
        op2v.read_op2(op2_filename, vectorized=True)

        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename, vectorized=False)
        assert os.path.exists('debug.out'), os.listdir('.')

        self._verify_ids(bdf, op2, vectorized=False, isubcase=1)
        self._verify_ids(bdf, op2v, vectorized=True, isubcase=1)

    def _verify_ids(self, bdf, op2, vectorized=True, isubcase=1):
        types = ['CQUAD4', 'CTRIA3', 'CHEXA', 'CPENTA', 'CTETRA', 'CROD', 'CONROD', 'CTUBE']
        out = bdf.get_card_ids_by_card_types(types)

        card_type = 'CQUAD4'
        if op2.cquad4_stress:
            case = op2.cquad4_stress[isubcase]
            eids = unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_strain:
            case = op2.cquad4_strain[isubcase]
            eids = unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_strain:
            case = op2.cquad4_composite_strain[isubcase]
            if vectorized:
                eids = unique(case.element_layer[:, 0])
            else:
                eids = unique(case.angle.keys())
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_stress:
            case = op2.cquad4_composite_stress[isubcase]
            if vectorized:
                eids = unique(case.element_layer[:, 0])
            else:
                eids = unique(case.angle.keys())
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_force:
            case = op2.cquad4_force[isubcase]
            if isinstance(case, RealPlateBilinearForce):
                eids = unique(case.tx.keys())
            else:
                assert isinstance(case, RealPlateBilinearForce), 'update this...'

            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)


        card_type = 'CTRIA3'
        if op2.ctria3_stress:
            case = op2.ctria3_stress[isubcase]
            eids = unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_strain:
            case = op2.ctria3_strain[isubcase]
            eids = unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_strain:
            case = op2.ctria3_composite_strain[isubcase]
            if vectorized:
                eids = unique(case.element_layer[:, 0])
            else:
                eids = unique(case.angle.keys())
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_stress:
            case = op2.ctria3_composite_stress[isubcase]
            if vectorized:
                eids = unique(case.element_layer[:, 0])
            else:
                eids = unique(case.angle.keys())
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_force:
            case = op2.ctria3_force[isubcase]
            if isinstance(case, RealPlateForceArray):
                eids = unique(case.element)
            else:
                assert isinstance(case, RealPlateForce), case
                eids = unique(case.tx.keys())
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)




    def test_op2_dmi(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'matrix', 'matrix.dat')
        op2_filename = os.path.join(folder, 'matrix', 'mymatrix.op2')
        matrices = {
            'A' : True,
            'B' : False,
            'ATB' : False,
            'BTA' : False,
            'MYDOF' : True,
        }
        from pyNastran.bdf.bdf import BDF
        model = BDF()
        model.read_bdf(bdf_filename)

        dmi_a = model.dmis['A']
        a, rows_reversed, cols_reversed = dmi_a.get_matrix(is_sparse=False, apply_symmetry=False)
        print('model.dmi.A =\n%s' % dmi_a)
        print('model.dmi.A =\n%s' % str(a))
        #return
        op2 = OP2()
        op2.set_additional_matrices_to_read(matrices)
        try:
            op2.read_op2(op2_filename)
            raise RuntimeError('this is wrong...')
        except FatalError:
            # the OP2 doesn't have a trailing zero marker
            pass

        from numpy import dot, array, array_equal
        # M rows, Ncols
        A = array([
            [1., 0.],
            [3., 6.],
            [5., 0.],
            [0., 8.],
        ], dtype='float32')
        B = A
        mydof = array([
            -1.0, 1.0, 1.0, -1.0, 1.0,
            2.0, -1.0, 1.0, 3.0, -1.0, 1.0, 4.0, -1.0,
            1.0, 5.0, -1.0, 1.0, 6.0, -1.0, 2.0, 1.0,
            -1.0, 2.0, 2.0, -1.0, 2.0, 3.0, -1.0, 2.0,
            4.0, -1.0, 2.0, 5.0, -1.0, 2.0, 6.0,
        ])
        BTA = dot(B.T, A)
        ATB = dot(A.T, B)

        expecteds = [A, ATB, B, BTA, mydof]
        matrix_names = sorted(matrices.keys())

        for table_name, expected in zip(matrix_names, expecteds):
            assert table_name in op2.matrices, table_name


            actual = op2.matrices[table_name].data
            if not array_equal(expected, actual):
                if table_name in model.dmis:
                    dmi = model.dmis[table_name]
                    table_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=False, apply_symmetry=False)
                    #stable_array, rows_reversed, cols_reversed = dmi.get_matrix(is_sparse=True, apply_symmetry=False)
                    print(table_array)
                #print(stable_array)
                msg = 'matrix %s was not read properly\n' % table_name
                msg += 'expected\n%s\n' % expected
                msg += 'actual\n%s' % actual
                print(msg)
                print('==========================')
                #raise RuntimeError(msg)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
