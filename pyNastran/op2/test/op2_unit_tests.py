from __future__ import print_function
import os
import unittest
import warnings

from six import iteritems
import numpy as np
from numpy import dot, array_equal
try:
    import pandas
    is_pandas = True
    # per http://stackoverflow.com/questions/35175949/ignore-pandas-warnings
    # doesn't work...
    #warnings.filterwarnings('ignore', '.*unorderable dtypes; returning scalar but in the future this will be an error.*')
except ImportError:
    is_pandas = False

import pyNastran

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2, FatalError, read_op2
from pyNastran.op2.op2_interface.op2_common import get_scode_word
from pyNastran.op2.op2_geom import read_op2_geom
from pyNastran.op2.test.test_op2 import run_op2

from pyNastran.bdf.test.bdf_unit_tests import Tester
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealPlateBilinearForceArray, RealPlateForceArray)
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
from pyNastran.op2.export_to_vtk import export_to_vtk_filename
from pyNastran.op2.vector_utils import filter1d
from pyNastran.utils.log import SimpleLogger

test_path = pyNastran.__path__[0]
model_path = os.path.abspath(os.path.join(test_path, '..', 'models'))
from pyNastran.op2.tables.geom.geom4 import read_rbe3s_from_idata_fdata

class TestOP2(Tester):
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_rbe3(self):
        """
        data = [99           99 123456 1.0    123    44    45  48  49  -1    -3]
        data = [61           71 123456 1.0    123    70    75  77      -1    -3
                62           71 123456 1.0    123    58    59  72      -1    -3]
        data = [1001100 1001100 123456 1.0 123456 10011 10002          -1 -2 -3
                1002500 1002500 123456 1.0 123456 10025 10020          -1 -2 -3]
                eid     refg    refc   wt  c      g     ...
        """
        model = BDF()
        data = [99, 99, 123456, 1.0, 123, 44, 45, 48, 49, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 1, rbes

        data = [61, 71, 123456, 1.0, 123, 70, 75, 77, -1, -3,
                62, 71, 123456, 1.0, 123, 58, 59, 72, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 2, rbes

        data = [11, 10, 123456, 1.0, 456, 111, 100, -1, -2, -3,
                25, 50, 123456, 1.0, 456, 125, 103, -1, -2, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 0, rbes  ## TODO: not supported b/c the -2; should be 2

    def test_filter1d(self):
        """tests filtering small values out of arrays"""
        a = np.array([1., 2., 0.1])
        i = filter1d(a, zero_tol=0.5)
        res = np.array([0, 1])
        self.assertTrue(np.array_equal(i, res), 'A i=%s res=%s' % (i, res))

        a = np.array([1., 2., 0.1])
        b = np.array([1., -0.1, 0.1])
        res = np.array([0, 1])
        i = filter1d(a, b, zero_tol=0.5)
        self.assertTrue(np.array_equal(i, res), 'B i=%s res=%s' % (i, res))

        a = np.array([1., 2., 0.1])
        b = np.array([1., -0.1, 0.1])
        i = filter1d(a, b, zero_tol=1.1)
        res = np.array([1])
        self.assertTrue(np.array_equal(i, res), 'C i=%s res=%s' % (i, res))

    def test_ibulk(self):
        """this test will fail if IBULK talble doesn't work"""
        bdf_filename = os.path.abspath(os.path.join(
            test_path, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        f06_filename = os.path.abspath(os.path.join(
            test_path, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.test_op2.f06'))
        op2_filename = os.path.abspath(os.path.join(
            test_path, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        op2 = read_op2(op2_filename, debug=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_beam_modes(self):
        """tests the Eigenvalue table reading"""
        f06_filename = os.path.abspath(os.path.join(
            model_path, 'beam_modes', 'model1_sim1-solution_1.test_op2.f06'))
        op2_filename_m1 = os.path.abspath(os.path.join(
            model_path, 'beam_modes', 'beam_modes_m1.op2'))
        op2_filename_m2 = os.path.abspath(os.path.join(
            model_path, 'beam_modes', 'beam_modes_m2.op2'))
        op2_1 = read_op2(op2_filename_m1, debug=False)
        op2_2 = read_op2(op2_filename_m2, debug=False)
        op2_1.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_set_results(self):
        op2_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')
        f06_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.test_op2.f06')

        op2 = OP2(debug=False)
        op2.set_results('stress')
        op2.read_op2(op2_filename)
        self.assertEqual(len(op2.cpenta_stress), 0, len(op2.cpenta_stress))
        self.assertEqual(len(op2.chexa_stress), 0, len(op2.chexa_stress))
        self.assertEqual(len(op2.ctetra_stress), 1, len(op2.ctetra_stress))
        self.assertEqual(len(op2.displacements), 0, len(op2.displacements))

        op2 = OP2(debug=False)
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2_filename)
        self.assertEqual(len(op2.cpenta_stress), 0, len(op2.cpenta_stress))
        self.assertEqual(len(op2.chexa_stress), 0, len(op2.chexa_stress))
        self.assertEqual(len(op2.ctetra_stress), 1, len(op2.ctetra_stress))
        self.assertEqual(len(op2.displacements), 1, len(op2.displacements))
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_solid_bending_01(self):
        folder = os.path.join(model_path, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        f06_filename = os.path.join(folder, 'solid_bending.test_op2.f06')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)

        read_op2(op2_filename)
        run_op2(op2_filename, write_bdf=write_bdf, subcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        assert os.path.exists(debug_file), os.listdir(folder)

        make_geom = False
        write_bdf = False
        write_f06 = True
        op2 = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)[0]
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_solid_bending_02(self):
        folder = os.path.join(model_path, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2 = read_op2(op2_filename, debug=False)

    def test_op2_solid_bending_02_geom(self):
        folder = os.path.join(model_path, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)

    def _test_op2_solid_bending_03(self):
        """tests basic op2 writing"""
        folder = os.path.join(model_path, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2_filename_debug = os.path.join(folder, 'solid_bending.debug.out')
        op2_filename_out = os.path.join(folder, 'solid_bending_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'solid_bending_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug)
        from pyNastran.op2.dev.op2_writer import OP2Writer
        op2w = OP2Writer()
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False,
                       delete_objects=True)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)


    def test_op2_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'static_solid_shell_bar.test_op2.f06')
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_mode_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        subcases = [1]
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        op2.get_op2_stats(short=False)
        op2.get_op2_stats(short=True)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_buckling_solid_shell_bar_01_geom(self):
        """single subcase buckling"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling_solid_shell_bar.op2')
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases)
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)

        f06_filename = os.path.join(folder, 'buckling_solid_shell_bar.test_op2_sort2.f06')
        op2.write_f06(f06_filename, is_mag_phase=False, is_sort1=False,
                      #delete_objects=True,
                      end_flag=False, quiet=True,
                      repr_check=False, close=True)

        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_buckling_solid_shell_bar_02_geom(self):
        """multi subcase buckling"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling2_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, debug=False)
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases)
        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 0, len(op2.eigenvectors)
        str(op2.iSubcaseNameMap)
        str(op2.displacements[1].subtitle)
        str(op2.displacements[1].label)

        subcases = 2
        op2, is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=True)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = 2
        op2, is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = [1, 2]
        op2, is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_transient_solid_shell_bar_01_geom(self):
        """transient test"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'transient_solid_shell_bar.test_op2.f06')
        op2, is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=False, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_frequency_solid_shell_bar_01_geom(self):
        """frequency test"""
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')
        op2 = read_op2_geom(op2_filename, debug=False)
        op2, is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=False, write_op2=False, write_xlsx=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_transfer_function_01(self):
        """tests the transfer function cards work"""
        folder = os.path.join(model_path, 'transfer_function')
        #bdf_filename = os.path.join(folder, 'actuator_tf_modeling.bdf')
        op2_filename = os.path.join(folder, 'actuator_tf_modeling.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')

        op2 = read_op2_geom(op2_filename, debug=False)

        debug = False
        write_bdf = True
        write_f06 = True
        make_geom = True
        op2, is_passed = run_op2(
            op2_filename, write_bdf=write_bdf, make_geom=make_geom, subcases=[],
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

        #fem1, fem2, diff_cards = self.run_bdf(folder, bdf_filename)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #for fem in [fem1, fem2]:
            #assert fem.card_count['CONM2'] == 3, fem.card_count
            #assert fem.card_count['SPOINT'] == 1, fem.card_count
            #assert fem.card_count['EPOINT'] == 1, fem.card_count
            #assert fem.card_count['PARAM'] == 1, fem.card_count
            #assert fem.card_count['CELAS2'] == 2, fem.card_count
            #assert fem.card_count['GRID'] == 3, fem.card_count
            #assert fem.card_count['EIGR'] == 1, fem.card_count
            #assert fem.card_count['EIGC'] == 1, fem.card_count
            #assert fem.card_count['MPC'] == 1, fem.card_count
            #assert fem.card_count['TF'] == 2, fem.card_count

    def test_monpnt3(self):
        """creates the MONPNT3 table"""
        folder = os.path.join(model_path, 'aero', 'monpnt3')
        op2_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.op2')
        f06_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.test_op2.f06')
        op2 = read_op2(op2_filename, debug=False)
        monitor3 = op2.monitor3
        assert len(monitor3.frequencies) == 11, monitor3
        str(monitor3)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_nastran_2005r3b(self):
        """Nastran2005r3 bug"""
        folder = os.path.join(model_path, 'modele_petite_zone')
        op2_filename = os.path.join(folder, 'modele_petite_zone.op2')
        f06_filename = os.path.join(folder, 'modele_petite_zone.test_op2.f06')
        op2 = read_op2(op2_filename, debug=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_solid_shell_bar_01(self):
        op2_filename = os.path.join('static_solid_shell_bar.op2')
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        rod_force = op2.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (1, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (1, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (1, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (1, 1, 15), cbar_stress.data.shape

        cbeam_force = op2.cbeam_force[isubcase]
        assert cbeam_force.nelements == 1, cbeam_force.nelements
        assert cbeam_force.data.shape == (1, 2, 8), cbeam_force.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # wrong
        assert cbeam_stress.data.shape == (1, 2, 8), cbeam_stress.data.shape

        cquad4_force = op2.cquad4_force[isubcase]
        assert cquad4_force.nelements == 4, cquad4_force.nelements
        assert cquad4_force.data.shape == (1, 20, 8), cquad4_force.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements # TODO: should this be 4; yes by actual count...
        assert cquad4_stress.data.shape == (1, 20, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert cquad4_stress.is_von_mises(), cquad4_stress

        ctria3_force = op2.ctria3_force[isubcase]
        assert ctria3_force.nelements == 8, ctria3_force.nelements
        assert ctria3_force.data.shape == (1, 8, 8), ctria3_force.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert ctria3_stress.is_von_mises(), ctria3_stress

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 10, 10), ctetra_stress.data.shape
        assert ctetra_stress.is_von_mises(), ctetra_stress

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (1, 14, 10), cpenta_stress.data.shape
        assert cpenta_stress.is_von_mises(), cpenta_stress

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (1, 9, 10), chexa_stress.data.shape
        assert chexa_stress.is_von_mises(), chexa_stress

        if is_pandas:
            rod_force.build_dataframe()
            rod_stress.build_dataframe()
            cbar_force.build_dataframe()
            cbar_stress.build_dataframe()
            cbeam_force.build_dataframe()
            cbeam_stress.build_dataframe()
            cquad4_force.build_dataframe()
            cquad4_stress.build_dataframe()
            ctria3_force.build_dataframe()
            ctria3_stress.build_dataframe()
            ctetra_stress.build_dataframe()
            cpenta_stress.build_dataframe()
            chexa_stress.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_01_export(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'static_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)
        os.remove(vtk_filename)

    def test_op2_solid_shell_bar_01_straincurvature(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert ctria3_stress.is_von_mises(), ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert cquad4_stress.is_von_mises(), cquad4_stress

        ctria3_strain = op2.ctria3_strain[isubcase]
        sword = get_scode_word(ctria3_strain.s_code, ctria3_strain.stress_bits)
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance(), '%s\n%s' % (ctria3_strain, sword)
        assert ctria3_strain.is_von_mises(), '%s\n%s' % (ctria3_strain, sword)

        cquad4_strain = op2.cquad4_strain[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance(), cquad4_strain
        assert cquad4_strain.is_von_mises(), cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert ctria3_stress.is_von_mises(), ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert cquad4_stress.is_von_mises(), cquad4_stress

        ctria3_strain = op2.ctria3_stress[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance(), ctria3_strain
        assert ctria3_strain.is_von_mises(), ctria3_strain

        cquad4_strain = op2.cquad4_stress[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), '%s\n%s' % (cquad4_strain.data.shape, sword)
        assert cquad4_strain.is_fiber_distance(), '%s\n%s' % (cquad4_strain, sword)
        assert cquad4_strain.is_von_mises(), '%s\n%s' % (cquad4_strain, sword)

    def test_op2_solid_shell_bar_01_straincurvature_shear(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve_shear.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve_shear.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert not ctria3_stress.is_von_mises(), ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert not cquad4_stress.is_von_mises(), cquad4_stress

        ctria3_strain = op2.ctria3_strain[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance(), ctria3_strain
        assert not ctria3_strain.is_von_mises(), ctria3_strain

        cquad4_strain = op2.cquad4_strain[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance(), cquad4_strain
        assert not cquad4_strain.is_von_mises(), cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance_shear(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist_shear.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist_shear.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance(), ctria3_stress
        assert ctria3_stress.is_von_mises() == False, ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance(), cquad4_stress
        assert cquad4_stress.is_von_mises() == False, cquad4_stress

        ctria3_strain = op2.ctria3_stress[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance(), ctria3_strain
        assert ctria3_strain.is_von_mises() == False, ctria3_strain

        cquad4_strain = op2.cquad4_stress[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert cquad4_strain.is_fiber_distance(), cquad4_strain
        assert cquad4_strain.is_von_mises() == False, cquad4_strain

    def test_op2_solid_shell_bar_mode(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        rod_force = op2.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (3, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (3, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        str(cbar_force.data_frame)
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (3, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 3, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (3, 1, 15), cbar_stress.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (3, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (3, 20, 8), cquad4_stress.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (3, 8, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (3, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (3, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (3, 9, 10), chexa_stress.data.shape

        if is_pandas:
            cbar_force.build_dataframe()
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_mode_export(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'mode_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'mode_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)
        os.remove(vtk_filename)

    def test_op2_solid_shell_bar_buckling(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcases = [(1, 1, 1, 0, 'DEFAULT1'), (1, 8, 1, 0, 'DEFAULT1')]
        isubcase = isubcases[1]

        try:
            rod_force = op2.crod_force[isubcase]
        except KeyError:
            msg = 'isubcase=%s was not found\nkeys=%s' % (isubcase, op2.crod_force.keys())
            raise KeyError(msg)
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (4, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (4, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (4, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 4, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (4, 1, 15), cbar_stress.data.shape

        cbeam_stress = op2.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (4, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (4, 20, 8), cquad4_stress.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (4, 8, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (4, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (4, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (4, 9, 10), chexa_stress.data.shape

        if is_pandas:
            cbar_force.build_dataframe()
            str(cbar_force.data_frame)

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    #@unittest.expectedFailure
    def test_op2_solid_shell_bar_freq(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1
        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        rod_force = op2.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape

        #if op2.is_nx:
            #return
        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (7, 2, 2), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (7, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (7, 1, 9), cbar_stress.data.shape

        #print(op2.cbeam_stress.keys())
        #cbeam_stress = op2.cbeam_stress[isubcase]
        #assert cbeam_stress.nelements == 2, cbeam_stress.nelements
        #assert cbeam_stress.data.shape == (7, 2, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: wrong
        assert cquad4_stress.data.shape == (7, 40, 3), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements # TODO: wrong
        assert ctria3_stress.data.shape == (7, 32, 3), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (7, 10, 6), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (7, 14, 6), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (7, 9, 6), chexa_stress.data.shape

        grid_point_forces = op2.grid_point_forces[isubcase]
        #print(grid_point_forces._ntotals)
        assert grid_point_forces.ntotal == 106, grid_point_forces.ntotal
        assert grid_point_forces.data.shape == (7, 106, 6), grid_point_forces.data.shape

        if is_pandas:
            rod_force.build_dataframe()
            rod_stress.build_dataframe()
            cbar_force.build_dataframe()
            cbar_stress.build_dataframe()
            #cquad4_stress.build_dataframe()
            ctria3_stress.build_dataframe()
            ctetra_stress.build_dataframe()
            cpenta_stress.build_dataframe()
            chexa_stress.build_dataframe()
            grid_point_forces.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_freq_export(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'freq_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'freq_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)
        os.remove(vtk_filename)

    def test_op2_solid_shell_bar_transient(self):
        """
        MSC 2005r2 Tables : GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPTS, DYNAMICS, DIT
                            OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        NX 10 Tables : PVT0, CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                       DYNAMICS, BGPDTS, EQEXINS, DIT,
                       OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        """
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=debug)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1
        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        rod_force = op2.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (21, 2, 2), rod_force.data.shape

        rod_stress = op2.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (21, 2, 4), rod_stress.data.shape

        cbar_force = op2.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (21, 1, 8), cbar_force.data.shape

        cbar_stress = op2.cbar_stress[isubcase]
        assert cbar_stress.nelements == 21, cbar_stress.nelements # 1-wrong
        assert cbar_stress.data.shape == (21, 1, 15), cbar_stress.data.shape

        #print(op2.cbeam_stress.keys())
        # cbeam_stress = op2.cbeam_stress[isubcase]
        # assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        # assert cbeam_stress.data.shape == (7, 11, 8), cbeam_stress.data.shape

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 40, cquad4_stress.nelements # TODO: (840-wrong, 40-correct)
        assert cquad4_stress.data.shape == (21, 40, 8), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 16, ctria3_stress.nelements # TODO: 8-wrong
        assert ctria3_stress.data.shape == (21, 16, 8), ctria3_stress.data.shape

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (21, 10, 10), ctetra_stress.data.shape

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (21, 14, 10), cpenta_stress.data.shape

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (21, 9, 10), chexa_stress.data.shape

        grid_point_forces = op2.grid_point_forces[isubcase]
        #print(grid_point_forces._ntotals)
        assert grid_point_forces.ntotal == 130, grid_point_forces.ntotal
        assert grid_point_forces.data.shape == (21, 130, 6), grid_point_forces.data.shape

        if is_pandas:
            rod_force.build_dataframe()
            rod_stress.build_dataframe()
            cbar_force.build_dataframe()
            cbar_stress.build_dataframe()
            #cquad4_stress.build_dataframe()
            ctria3_stress.build_dataframe()
            ctetra_stress.build_dataframe()
            cpenta_stress.build_dataframe()
            chexa_stress.build_dataframe()
            grid_point_forces.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_optistruct_01(self):
        """
        Optistruct 2012 Tables : CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                                OUGV1, OES1X
        """
        op2_filename = os.path.join(model_path, 'optistruct', 'hm14.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=debug)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1
        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        #assert len(op2.rod_force) == 0
        assert len(op2.crod_stress) == 0

        assert len(op2.cbar_force) == 0
        assert len(op2.cbar_stress) == 0
        assert len(op2.cbeam_stress) == 0

        assert len(op2.cquad4_stress) == 0
        assert len(op2.ctria3_stress) == 0

        ctetra_stress = op2.ctetra_stress[isubcase]
        if is_pandas:
            ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 3951, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 19755, 10), ctetra_stress.data.shape

        assert len(op2.cpenta_stress) == 0
        assert len(op2.chexa_stress) == 0

        assert len(op2.grid_point_forces) == 0
        os.remove(debug_file)

    def test_op2_plate_py_01(self):
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        op2_filename = os.path.join(model_path, 'plate_py', 'plate_py.op2')
        read_op2(op2_filename)
        run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, quiet=True)

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, quiet=True)

    def test_op2_good_sine_01(self):
        op2_filename = os.path.join(model_path, 'freq_sine', 'good_sine.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        read_op2(op2_filename)
        op2i, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                  write_f06=write_f06,
                                  debug=debug, stop_on_failure=True,
                                  quiet=True)

        nids = [5]
        isubcase = 103
        try:
            acc = op2i.accelerations[isubcase]
        except KeyError:
            raise KeyError('getting accelerations; isubcase=%s; keys=%s' % (
                isubcase, op2i.accelerations.keys()))

        with self.assertRaises(AssertionError):
            # no index 0; fortran 1-based
            acc.extract_xyplot(nids, 0, 'real')

        if is_pandas:
            acc.build_dataframe()
        accx = acc.extract_xyplot(nids, 1, 'real')
        accxi = acc.extract_xyplot(nids, 1, 'imag')
        #print(accx)
        #print(accxi)
        #make_geom = False
        #write_bdf = False
        #write_f06 = True
        #run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                #write_f06=write_f06,
                #debug=debug, stopOnFailure=True)

    def test_op2_good_sine_02(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'freq_sine', 'good_sine.dat')
        op2_filename = os.path.join(folder, 'freq_sine', 'good_sine.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename)

        debug = False
        debug_file = 'debug.out'

        read_op2(op2_filename)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')

        self._verify_ids(bdf, op2, isubcase=1)
        os.remove(debug_file)

    def test_op2_bcell_01(self):
        folder = os.path.abspath(os.path.join(test_path, '..', 'models'))
        bdf_filename = os.path.join(folder, 'other', 'bcell9p0.bdf')
        op2_filename = os.path.join(folder, 'other', 'bcell9p0.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename, xref=False)

        debug = False
        debug_file = 'debug.out'
        read_op2(op2_filename)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')
        os.remove(debug_file)

        self._verify_ids(bdf, op2, isubcase=1)

        msg = ''
        for isubcase, keys in sorted(iteritems(op2.subcase_key)):
            if len(keys) != 1:
                msg += 'isubcase=%s keys=%s len(keys) != 1\n' % (isubcase, keys)
                if len(keys) == 0:
                    continue
            if isubcase != keys[0]:
                msg += 'isubcase=%s != key0=%s keys=%s\n' % (isubcase, keys[0], keys)
        if msg:
            assert msg == '', msg
        op2.write_f06('junk.f06', quiet=True)
        os.remove('junk.f06')

    def test_op2_cbush_01(self):
        op2_filename = os.path.join(model_path, 'cbush', 'cbush.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model, ext = os.path.splitext(op2_filename)
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename)
        op2, is_passed = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf, subcases=[],
                                 write_f06=write_f06,
                                 debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)
        isubcase = 1

        cbush_stress = op2.cbush_stress[isubcase]
        assert cbush_stress.nelements == 1, cbush_stress.nelements
        assert cbush_stress.data.shape == (1, 1, 6), cbush_stress.data.shape

        cbush_strain = op2.cbush_strain[isubcase]
        assert cbush_strain.nelements == 1, cbush_strain.nelements
        assert cbush_strain.data.shape == (1, 1, 6), cbush_strain.data.shape

        cbush_force = op2.cbush_force[isubcase]
        assert cbush_force.nelements == 1, cbush_force.nelements
        assert cbush_force.data.shape == (1, 1, 6), cbush_force.data.shape
        if is_pandas:
            cbush_stress.build_dataframe()
            cbush_strain.build_dataframe()
            cbush_force.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def _verify_ids(self, bdf, op2, isubcase=1):
        types = ['CQUAD4', 'CTRIA3', 'CHEXA', 'CPENTA', 'CTETRA', 'CROD', 'CONROD', 'CTUBE']
        out = bdf.get_card_ids_by_card_types(types)

        card_type = 'CQUAD4'
        if op2.cquad4_stress:
            try:
                case = op2.cquad4_stress[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_stress; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_stress.keys()))
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_strain:
            try:
                case = op2.cquad4_strain[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_strain; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_strain.keys()))
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_strain:
            try:
                case = op2.cquad4_composite_strain[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_composite_strain; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_composite_strain.keys()))
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_composite_stress:
            try:
                case = op2.cquad4_composite_stress[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_composite_stress; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_composite_stress.keys()))

            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.cquad4_force:
            try:
                case = op2.cquad4_force[isubcase]
            except KeyError:
                raise KeyError('getting cquad4_force; isubcase=%s; keys=%s' % (
                    isubcase, op2.cquad4_force.keys()))

            eids = np.unique(case.element)
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)


        card_type = 'CTRIA3'
        if op2.ctria3_stress:
            case = op2.ctria3_stress[isubcase]
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_strain:
            case = op2.ctria3_strain[isubcase]
            eids = np.unique(case.element_node[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_strain:
            case = op2.ctria3_composite_strain[isubcase]
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_composite_stress:
            case = op2.ctria3_composite_stress[isubcase]
            eids = np.unique(case.element_layer[:, 0])
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
        if op2.ctria3_force:
            case = op2.ctria3_force[isubcase]
            eids = np.unique(case.element)
            for eid in eids:
                assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)

    @unittest.expectedFailure
    def test_set_times_01(self):
        """specify the modes to extract"""
        model = OP2(debug=True, log=None, debug_file=None, mode='msc')

        # subcase : times
        times = {1 : [2, 3]}
        model.set_transient_times(times)
        op2_filename = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements',
                                                    'mode_solid_shell_bar.op2'))
        model.read_op2(op2_filename, combine=True, build_dataframe=None,
                       skip_undefined_matrices=False,
                       encoding=None)
        isubcase = 1
        #print(model.get_op2_stats(short=False))
        eigenvector = model.eigenvectors[isubcase]
        #print(eigenvector)
        assert len(eigenvector.modes) == 2, eigenvector.modes


class TestGridPointForces(unittest.TestCase):
    """various grid point force tests"""

    def test_gpforce_01(self):
        nids = np.array([1, 2, 3])
        xyz_cid0 = np.array([
            [1., 1., 1.],
            [4., 2., 5.],
            [3., 3., 3.],
        ])
        data_code = {
            'nonlinear_factor' : None,
            'sort_bits' : [0, 0, 0],
            'analysis_code' : 1,
            'is_msc' : True,
            'format_code' : 1,
            'table_code' : 1,
            'data_names' : 'cat',
            'device_code' : 1,
            #'tcode' : 1,
        }
        is_sort1 = True
        isubcase = 1
        dt = 0.0
        gpforce = RealGridPointForcesArray(data_code, is_sort1, isubcase, dt)
        gpforce.ntimes = 1
        gpforce.ntotal = 3
        gpforce._ntotals = [3]

        gpforce.build()
        gpforce.data[0, :, :] = np.array([
            [3., 7., 11., 0., 0., 0.,], # fx, fy, fz, mx, my, mz
            [3., 7., 11., 0., 0., 0.,],
            [3., 7., 11., 0., 0., 0.,],
        ])
        gpforce.node_element[0, :, :] = np.array([
            [1, 1],
            [2, 1],
            [3, 1],
        ])
        op2 = OP2()
        summation_point = [0., 0., 0.]
        i_transform = None
        nid_cd = np.array([
            [1, 0],
            [2, 0],
            [3, 0],
        ])
        from pyNastran.bdf.bdf import CORD2R
        coord_out = CORD2R(cid=0)
        coords = {0 : coord_out}

        #eids = [1]
        #nids = [1]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #i_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #logger=op2.log)

        #print('------------')
        #eids = [1]
        #nids = [2]
        #gpforce.extract_interface_loads(
            #nids, eids, coord_out, coords, nid_cd,
            #i_transform,
            #xyz_cid0,
            #summation_point,
            #itime=0,
            #debug=True,
            #logger=op2.log)
        print('------------')

        eids = [1]
        nids = [1, 2]
        gpforce.extract_interface_loads(
            nids, eids, coord_out, coords, nid_cd,
            i_transform,
            xyz_cid0,
            summation_point,
            itime=0,
            debug=True,
            logger=op2.log)
        #print(gpforce)

    def test_op2_solid_shell_bar_01_gpforce(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2 = read_op2_geom(op2_filename, xref=False, debug=False)

        i_transform, beta_transforms = op2.get_displacement_index_transforms()
        op2.transform_displacements_to_global(i_transform, op2.coords)

        gpforce = op2.grid_point_forces[1]

        #bdf_filename = os.path.join(folder, 'solid_shell_bar_xyz.bdf')
        #model = BDF(debug=False)
        #model.read_bdf(bdf_filename, xref=True)

        op2.cross_reference(xref_elements=False,
                            xref_nodes_with_elements=False,
                            xref_properties=False,
                            xref_masses=False,
                            xref_materials=False,
                            xref_loads=False,
                            xref_constraints=False,
                            xref_aero=False,
                            xref_sets=False,
                            xref_optimization=False)
        xyz_cid0 = op2.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2.nodes))])
        coords = op2.coords

        data = _get_gpforce_data()
        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            if cid not in coords:
                continue
            #op2.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            coord_out = coords[cid]
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            #op2.log.debug('***********')
            #op2.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            #op2.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s' % (
                case, total_force_local_expected, total_force_local)
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.005), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s' % (
                case, total_moment_local_expected, total_moment_local)
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    def test_op2_solid_shell_bar_01_gpforce_xyz(self):
        folder = os.path.join(model_path, 'sol_101_elements')
        #bdf_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.bdf')
        op2_filename1 = os.path.join(folder, 'static_solid_shell_bar_xyz.op2')
        op2_1 = read_op2_geom(op2_filename1, xref=False, debug=False)
        #print("disp_orig =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        nids_all, nids_transform_1, i_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(i_transform_1, op2_1.coords)
        op2_1.transform_gpforce_to_global(nids_all, nids_transform_1, i_transform_1, op2_1.coords)
        #print("disp_new =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        gpforce = op2_1.grid_point_forces[1]
        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)
        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2_1.nodes))])

        #bdf_filename2 = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, i_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(i_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(nids_all, nids_transform_2, i_transform_2, op2_2.coords)

        #print("disp_goal =\n", op2_2.displacements[1].data[0, :2, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, :2, :])

        msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04), msg
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1]), msg
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1]), msg
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123), msg

        #-------------------------------------------------
        return
        data = _get_gpforce_data()
        coords = op2_1.coords
        used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        for datai in data:
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test' + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, coords,
                nid_cd, i_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    @unittest.expectedFailure
    def test_op2_solid_shell_bar_01_gpforce_radial_global_cd(self):
        warning_log = SimpleLogger(level='warning')
        debug_log = SimpleLogger(level='debug')
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename1 = os.path.join(folder, 'static_solid_shell_bar_global_radial_cd.op2')
        op2_1 = read_op2_geom(op2_filename1, xref=False, log=warning_log)
        op2_1.log = debug_log

        #print("disp_orig =\n", op2_1.displacements[1].data[0, :4, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, 2:8, :])

        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)

        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2_1.nodes))])
        #-------------------------------------------------
        #coords = op2_1.coords
        #used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))

        nids_all, nids_transform_1, i_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(i_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        op2_1.transform_gpforce_to_global(nids_all, nids_transform_1, i_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        #print('stuff...')
        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))
        #print(disp.data[0, :, :3])

        #print("disp_new =\n", op2_1.displacements[1].data[0, :4, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, 2:8, :])

        #-----------------------------------------------------------------------
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, i_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(i_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(nids_all, nids_transform_2, i_transform_2, op2_2.coords, xyz_cid0=xyz_cid0)

        #print("disp_goal =\n", op2_2.displacements[1].data[0, :4, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, 2:8, :])

        #return
        #msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            #op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        #assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        #msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            #op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        from six import StringIO
        csv_file = StringIO()
        op2_1.spc_forces[1].write_csv(csv_file)
        op2_1.grid_point_forces[1].write_csv(csv_file)
        #print(csv_file.getvalue())

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04), msg
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1]), msg
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1]), msg
        #print('op2_2.grid_point_forces[1]\n', op2_2.grid_point_forces[1].data)
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123), msg
        return
        #-----------------------------------------------------------------------
        gpforce = op2_1.grid_point_forces[1]
        data = _get_gpforce_data()
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = op2_1.coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, op2_1.coords,
                nid_cd, i_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

    #@unittest.expectedFailure
    def test_op2_solid_shell_bar_01_gpforce_radial(self):
        warning_log = SimpleLogger(level='warning')
        debug_log = SimpleLogger(level='debug')
        folder = os.path.join(model_path, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_radial.op2')
        op2_1 = read_op2_geom(op2_filename, xref=False, log=warning_log)
        op2_1.log = debug_log

        print("disp_orig =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_orig =\n", op2_1.spc_forces[1].data[0, -3:, :])
        #print("gpf_orig =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        op2_1.cross_reference(xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        xyz_cid0 = op2_1.get_xyz_in_coord(cid=0)

        nid_cd = np.array([[nid, node.Cd()] for nid, node in sorted(iteritems(op2_1.nodes))])
        #-------------------------------------------------
        #coords = op2_1.coords
        #used_cds = np.unique(nid_cd[:, 1])
        #for cd in used_cds:
            #coord = op2_1.coords[cd]
            #print(coord)
            #print('origin = %s' % coord.origin)
            #print('beta =\n%s' % coord.beta())
            #print('-----------------------------')

        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))

        nids_all, nids_transform_1, i_transform_1 = op2_1.get_displacement_index()
        op2_1.transform_displacements_to_global(i_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        op2_1.transform_gpforce_to_global(
            nids_all, nids_transform_1, i_transform_1, op2_1.coords, xyz_cid0=xyz_cid0)
        #print('stuff...')
        #disp = op2_1.displacements[1]
        #for line in list(disp.data[0, :, :3]):
            #print('%10.4e %10.4e %10.4e' % tuple(line))
        #print(disp.data[0, :, :3])

        print("disp_new =\n", op2_1.displacements[1].data[0, :2, :])
        #print("spc_new =\n", op2_1.spc_forces[1].data[0, -3:, :])
        #print("gpf_new =\n", op2_1.grid_point_forces[1].data[0, :2, :])

        #-----------------------------------------------------------------------
        op2_filename2 = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_2 = read_op2_geom(op2_filename2, debug=False)
        nids_all, nids_transform_2, i_transform_2 = op2_2.get_displacement_index()
        op2_2.transform_displacements_to_global(i_transform_2, op2_2.coords)
        op2_2.transform_gpforce_to_global(
            nids_all, nids_transform_2, i_transform_2, op2_2.coords)

        print("disp_goal =\n", op2_2.displacements[1].data[0, :2, :])
        #print("spc_goal =\n", op2_2.spc_forces[1].data[0, -3:, :])
        #print("gpf_goal =\n", op2_2.grid_point_forces[1].data[0, :2, :])

        return
        msg = 'displacements baseline=\n%s\ndisplacements xyz=\n%s' % (
            op2_1.displacements[1].data[0, :, :], op2_2.displacements[1].data[0, :, :])
        #print(msg)
        assert op2_1.displacements[1].assert_equal(op2_2.displacements[1])

        msg = 'grid_point_forces baseline=\n%s\ngrid_point_forces xyz=\n%s' % (
            op2_1.grid_point_forces[1].data[0, :, :], op2_2.grid_point_forces[1].data[0, :, :])
        #print(msg)

        assert op2_1.spc_forces[1].assert_equal(op2_2.spc_forces[1], atol=4.4341e-04), msg
        assert op2_1.mpc_forces[1].assert_equal(op2_2.mpc_forces[1]), msg
        assert op2_1.load_vectors[1].assert_equal(op2_2.load_vectors[1]), msg
        assert op2_1.grid_point_forces[1].assert_equal(op2_2.grid_point_forces[1], atol=0.000123), msg
        #-----------------------------------------------------------------------
        gpforce = op2_1.grid_point_forces[1]
        data = _get_gpforce_data()
        for i, datai in enumerate(data):
            eids, nids, cid, summation_point, total_force_local_expected, total_moment_local_expected = datai
            coord_out = op2_1.coords[cid]
            op2_1.log.debug('*' * 30 + 'Next Test #%s' % i + '*' * 30)
            out = gpforce.extract_interface_loads(
                nids, eids,
                coord_out, op2_1.coords,
                nid_cd, i_transform_1,
                xyz_cid0, summation_point, itime=0, debug=False, logger=op2_1.log)
            total_force_global, total_moment_global, total_force_local, total_moment_local = out

            op2_1.log.debug('***********')
            op2_1.log.debug('force = %s; %s' % (total_force_global, np.linalg.norm(total_force_global)))
            op2_1.log.debug('moment = %s; %s' % (total_moment_global, np.linalg.norm(total_moment_global)))

            case = 'eids=%s nids=%s cid=%s summation_point=%s' % (
                eids, nids, cid, summation_point)
            msg = '%s\ntotal_force_local_expected=%s total_force_local=%s delta=%s' % (
                case, total_force_local_expected, total_force_local,
                np.abs(total_force_local_expected - total_force_local))
            self.assertTrue(np.allclose(total_force_local_expected, total_force_local, atol=0.2), msg), msg

            msg = '%s\ntotal_moment_local_expected=%s total_moment_local=%s delta=%s' % (
                case, total_moment_local_expected, total_moment_local,
                np.abs(total_moment_local_expected - total_moment_local))
            self.assertTrue(np.allclose(total_moment_local_expected, total_moment_local, atol=0.005), msg), msg

def _get_gpforce_data():
    data = [
        #eids, nids, cid, summation_point
        #[1], [1], 0, [0., 0., 0.],
        #[[1], [1, 2, 3, 4], 0, [0., 0., 0.], [0.0, 0.0, -10000.0], [-5000.0, 5000.0, 0.0],],  # total; good for gpforce

        # cid=0; eid=[1]; nid=[3]; sum=[0., 0., 0.] - done
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2589.95,     0.0,  26.34, -44.15, -2589.44,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2589.95, 3862.70,  26.34, -44.15, -2589.44, -2589.44, 2589.44, -70.49]  # total
        #[[1], [3], 0, [0., 0., 0.], [26.34, -44.15, -2589.44], [-2589.44, 2589.44, -70.49],], # good for gpforce; failing for xyz (cid 11)

        # cid=0; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #                            fx      fy       fz       mx       my       mz
        [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment; good for gpforce; failing for xyz (cid 11)

        # cid=0/1/2/3; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        [[1], [1], 0, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 1, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 2, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment
        [[1], [1], 3, [0., 0., 0.], [-37.18, 32.00, -2589.44], [0.0, 0.0, 0.0],],  # only 1 line b/c no moment

        # cid=1; eid=[1]; nid=[1]; sum=[0., 0., 0.]
        #               fmag     mmag      fx      fy       fz       mx       my       mz
        # F1      = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # ith
        # F1Total = [2589.90,     0.0,1853.64, 567.74, 1717.35,     0.0,     0.0,     0.0]  # total
        [[1], [1], 11, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce
        [[1], [1], 12, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce
        [[1], [1], 13, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0.0, 0.0, 0.0],], # good; failing for gpforce

        # cid=1; eid=[1]; nid=[2]; sum=[0., 0., 0.]
        #               fmag     mmag       fx      fy       fz       mx       my       mz
        # F2      = [2411.67,     0.0, 1710.67, 634.80, 1577.03,     0.0,     0.0,     0.0]  # ith
        # F2Total = [2411.67, 2410.58, 1710.67, 634.80, 1577.03, 1698.38, -570.22, -1612.84] # total
        [[1], [2], 11, [0., 0., 0.], [1710.67, 634.60, 1577.03], [1698.38, -570.22, -1612.84],],

        # cid=1; eid=[1]; nid=[3]; sum=[0., 0., 0.]
        #           fmag          mmag     fx       fy       fz       mx        my       mz
        # F3      = [2589.95,     0.0, 1799.79, 645.58, 1746.94,     0.0,      0.0,     0.0]  # ith
        # F3Total = [2589.95, 3862.70, 1799.79, 645.58, 1746.94, 1880.85, -3035.07, -816.15]  # total
        [[1], [3], 11, [0., 0., 0.], [1799.79, 645.58, 1746.94], [1880.85, -3035.07, -816.15]],

        #[[1], [1], 11, [0., 0., 0.], [1853.64, 567.74, 1717.35], [0., 0., 0.],],
        #[[1], [1], 12, [0., 0., 0.], [1938.05, -47.57, 1717.35], [0., 0., 0.],],
        #[[1], [1], 13, [0., 0., 0.], [2069.00, 1557.11, -47.57], [0., 0., 0.],],
    ]
    return data


if __name__ == '__main__':  # pragma: no cover
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
    if not on_rtd:
        unittest.main()
