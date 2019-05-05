from __future__ import print_function
import os
import unittest
import numpy as np
from cpylog import get_logger

try:
    import pandas  # pylint: disable=unused-import
    IS_PANDAS = True
    # per http://stackoverflow.com/questions/35175949/ignore-pandas-warnings
    # doesn't work...
    #warnings.filterwarnings(
        #'ignore',
        #'.*unorderable dtypes; returning scalar but in the future this will be an error.*')
except ImportError:
    IS_PANDAS = False

try:
    import h5py  # pylint: disable=unused-import
    IS_H5PY = True
except ImportError:  # pragma: no cover
    IS_H5PY = False

IS_TRANSIENT_PANDAS = False
if IS_PANDAS and (np.lib.NumpyVersion(np.__version__) < '1.13.0'):
    IS_TRANSIENT_PANDAS = True

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.op2.op2_interface.op2_common import get_scode_word
from pyNastran.op2.op2_geom import OP2Geom, read_op2_geom
from pyNastran.op2.test.test_op2 import run_op2, main as test_op2

from pyNastran.bdf.test.bdf_unit_tests import Tester
#from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    #RealPlateBilinearForceArray, RealPlateForceArray)
#from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
from pyNastran.op2.export_to_vtk import export_to_vtk_filename
from pyNastran.op2.vector_utils import filter1d, abs_max_min_global, abs_max_min_vector
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
from pyNastran.femutils.test.utils import is_array_close

from pyNastran.op2.tables.geom.geom4 import _read_spcadd_mpcadd

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))


class TestOP2(Tester):
    """various OP2 tests"""
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_cd_displacement(self):
        log = get_logger(level='debug')
        data_code = {
            'device_code' : 1,
            'analysis_code' : 1,
            'table_code' : 1,
            'nonlinear_factor' : None,
            'sort_bits' : [0, 0, 0],
            'sort_method' : 1,
            'is_msc' : True,
            'format_code' : 1,
            'data_names' : [],
            'tCode' : 1,
            'table_name' : 'OUGV1',
            '_encoding' : 'utf-8',
        }

        bdf_model = BDF(log=log)
        bdf_model.add_grid(1, [0., 0., 0.], cd=0)
        bdf_model.add_grid(2, [0., 0., 0.], cd=1)
        bdf_model.add_grid(3, [0., 0., 0.], cp=2, cd=2)

        bdf_model.add_grid(11, [1., 0., 0.], cd=0)
        bdf_model.add_grid(12, [1., 0., 0.], cd=1)
        bdf_model.add_grid(13, [1., 0., 0.], cp=2, cd=2)

        #bdf_model.add_grid(21, [1., 0., 0.], cd=0)
        #bdf_model.add_grid(22, [1., 0., 0.], cd=1)
        bdf_model.add_grid(23, [1., 90., 0.], cp=2, cd=2)
        bdf_model.add_grid(24, [0., 1., 0.], cp=0, cd=2)
        bdf_model.add_grid(25, [-1., 0., 0.], cp=0, cd=2)

        #bdf_model.add_grid(31, [1., 0., 0.], cp=3, cd=3)  # [0,1,0]
        #bdf_model.add_grid(32, [1., 90., 0.], cp=3, cd=3) # [0,-1,0]

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        bdf_model.add_cord2r(1, origin, zaxis, xzplane, rid=0, comment='')
        unused_coord = bdf_model.add_cord2c(2, origin, zaxis, xzplane, rid=0, comment='')

        origin = [0., 0., 0.]
        zaxis = [1., 0., 0.]
        xzplane = [0., 1., 0.]
        bdf_model.add_cord2c(3, origin, zaxis, xzplane, rid=0, comment='')

        dxyz = np.array([[
            [1., 0., 0., 0., 0., 0.], # 1
            [1., 0., 0., 0., 0., 0.], # 2
            [1., 0., 0., 0., 0., 0.], # 3

            [1., 0., 0., 0., 0., 0.], # 11
            [1., 0., 0., 0., 0., 0.], # 12
            [1., 0., 0., 0., 0., 0.], # 13

            [1., 0., 0., 0., 0., 0.], # 23 - [0., 1., 0.]
            [1., 0., 0., 0., 0., 0.], # 24 - answer=same as 23
            [1., 0., 0., 0., 0., 0.], # 25 - [-1, 0., 0.]

            #[1., 0., 0., 0., 0., 0.], # 31 - [0,1,0]
            #[1., 0., 0., 0., 0., 0.], # 32 - [0,-1,0]
        ]])
        #--------------------------------------------
        #out = bdf_model.get_displacement_index_xyz_cp_cd(
            #fdtype='float64', idtype='int32', sort_ids=True)
        #out = icd_transform, icp_transform, xyz_cp, nid_cp_cd
        out = bdf_model.get_xyz_in_coord_array(
            cid=0, fdtype='float64', idtype='int32')
        unused_nid_cp_cd, xyz_cid0, unused_xyz_cp, icd_transform, unused_icp_transform = out

        op2_model = OP2(log=log)

        is_sort1 = True
        isubcase = 1
        dt = None
        disp = RealDisplacementArray(data_code, is_sort1, isubcase, dt)
        disp.data = dxyz
        op2_model.displacements[1] = disp

        op2_model.transform_displacements_to_global(
            icd_transform, bdf_model.coords, xyz_cid0=xyz_cid0, debug=True)


        # we're working in a 2D plane
        icd2 = icd_transform[2]
        unused_dispi_cd2 = op2_model.displacements[1].data[0, icd2, :2]
        #op2_model.log.info("dispi2:\n%s" % dispi_cd2)

        dispi = op2_model.displacements[1].data[0, :, :2]
        expected_disp = np.array([
            [1., 0.,], # 1
            [1., 0.,], # 2
            [1., 0.,], # 3

            [1., 0.,], # 11
            [1., 0.,], # 12
            [1., 0.,], # 13

            [0., 1.,], # 23
            [0., 1.,], # 24
            [-1., 0.,], # 25

            #[0., 1.,], # 31
            #[0., -1.,], # 32
        ])
        assert is_array_close(dispi, expected_disp)
        #print(is_array_close(dispi, expected_disp))
        #print(dispi)

        ## TODO: fix the thetad in the cid=3 coordinates (nid=33,34)

    def test_generalized_tables(self):
        """tests that set_additional_generalized_tables_to_read overwrites the GEOM1S class"""
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        model = OP2Geom()
        model.read_op2(op2_filename=op2_filename, combine=True, build_dataframe=None,
                       skip_undefined_matrices=False,
                       encoding=None)

        def read_some_table(self):
            """crashes"""
            raise NotImplementedError('read_some_table')

        model2 = OP2Geom()
        tables = {
            b'GEOM1S' : read_some_table,
        }
        model2.set_additional_generalized_tables_to_read(tables)
        with self.assertRaises(NotImplementedError):
            model2.read_op2(op2_filename=op2_filename, combine=True, build_dataframe=None,
                            skip_undefined_matrices=False,
                            encoding=None)


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

    def test_abs_max_min_global(self):
        #print(iformat('4si3f', 2))
        str(abs_max_min_global([0.0, 2.0, 1.0]))
        str(abs_max_min_global([0.0, 2.0, -1.0]))
        str(abs_max_min_global([0.0, 2.0, -3.0]))
        str(abs_max_min_global(np.array([0.0, 2.0, -3.0])))
        str(abs_max_min_global([1.0]))

        # gets the global max/min value
        str(abs_max_min_global([
            [0.0, 2.0, -3.0],
            [0.0, 2.0, -4.0],
        ]))
        str(abs_max_min_global(np.array([
            [0.0, 2.0, -3.0],
            [0.0, 2.0, -4.0],
        ])))

    def test_abs_max_min_vector(self):
        str(abs_max_min_vector(np.array([
            [0.0, 2.0, 1.0],
            [0.0, 2.0, -1.0],
            [0.0, 2.0, -3.0],
        ])))

        str(abs_max_min_vector([
            [0.0, 2.0, 1.0],
            [0.0, 2.0, -1.0],
            [0.0, 2.0, -3.0],
            [0.0, 2.0, 4.0],
        ]))
        str(abs_max_min_vector(np.array([
            [0.0, 2.0, 1.0],
            [0.0, 2.0, -1.0],
            [0.0, 2.0, -3.0],
            [0.0, 2.0, 4.0],
        ])))

        str(abs_max_min_vector(np.array([
            [3.0, 2.0, -3.0],
            [-3.0, 2.0, 3.0],
        ])))

        # not an array
        #print(abs_max_min([
            #[0.0, 2.0, 1.0],
            #[0.0, 2.0, -1.0],
            #[0.0, 2.0, -3.0],
            #[0.0, 2.0, 4.0],
        #]))

    def test_ibulk(self):
        """this test will fail if IBULK talble doesn't work"""
        unused_bdf_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        f06_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.test_op2.f06'))
        op2_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        op2 = read_op2_geom(op2_filename, xref=False, debug=False, debug_file='temp.debug')
        op2.write_f06(f06_filename)
        os.remove(f06_filename)
        os.remove('temp.debug')

    def test_beam_modes(self):
        """tests the eigenvalue table reading"""
        f06_filename = os.path.abspath(os.path.join(
            MODEL_PATH, 'beam_modes', 'model1_sim1-solution_1.test_op2.f06'))
        op2_filename_m1 = os.path.abspath(os.path.join(
            MODEL_PATH, 'beam_modes', 'beam_modes_m1.op2'))
        op2_filename_m2 = os.path.abspath(os.path.join(
            MODEL_PATH, 'beam_modes', 'beam_modes_m2.op2'))
        op2_1 = read_op2(op2_filename_m1, debug=False)
        unused_op2_2 = read_op2_geom(op2_filename_m2, debug=False, debug_file='temp.debug')
        op2_1.write_f06(f06_filename)
        os.remove(f06_filename)
        os.remove('temp.debug')

    def test_bdf_op2_elements_01(self):
        """tests a large number of elements and results in SOL 101"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True, check_memory=False,
            stop_on_failure=True, dev=False)
        with self.assertRaises(NotImplementedError):
            op2.save()

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True, check_memory=False,
            stop_on_failure=True, dev=False)
        op2.save()

    def test_bdf_op2_elements_02(self):
        """tests a large number of elements and results in SOL 103-modes"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        print(bdf_filename)
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)
        run_op2(op2_filename, make_geom=True, write_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False)


    def test_bdf_op2_elements_03(self):
        """tests a large number of elements and results in SOL 108-freq"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_04(self):
        """tests a large number of elements and results in SOL 108-freq"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        build_pandas = IS_TRANSIENT_PANDAS
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False, build_pandas=build_pandas)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_05(self):
        """tests a large number of elements and results in SOL 106-loadstep"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        build_pandas = IS_TRANSIENT_PANDAS
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False, build_pandas=build_pandas)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_06(self):
        """tests a large number of elements and results in SOL 107-complex modes"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_post_minus4(self):
        """tests a large number of elements and results in SOL 107-complex modes"""
        unused_bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements_post4.op2')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements_post4.op2')
        #fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #read_op2(op2_filename=op2_filename, combine=True, subcases=None,
                 #exclude_results=None, include_results=None,
                 #log=None, debug=True, debug_file=None,
                 #build_dataframe=None,
                 #skip_undefined_matrices=True, mode='msc',
                 #encoding=None)
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False, post=-4)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_01(self):
        """checks time_thermal_elements.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_02(self):
        """checks hd15306.bdf"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        #op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')

        bdf_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False,
                build_pandas=False)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_03(self):
        """checks time_thermal_elements.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False,
                build_pandas=False)

    def test_bdf_op2_other_01(self):
        """checks ofprand1.bdf"""
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        #op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')

        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'ofprand1.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ofprand1.op2')
        #fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False, xref_safe=True,
                write_f06=False, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True, check_memory=False,
                stop_on_failure=True, dev=False,
                build_pandas=False)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_set_results(self):
        """tests setting only a subset of results"""
        op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.op2')
        f06_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.test_op2.f06')

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
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        f06_filename = os.path.join(folder, 'solid_bending.test_op2.f06')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)

        read_op2(op2_filename, debug=False)
        run_op2(op2_filename, write_bdf=write_bdf,
                write_f06=write_f06,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
                load_as_h5=False)
        assert os.path.exists(debug_file), os.listdir(folder)

        make_geom = False
        write_bdf = False
        write_f06 = True
        build_pandas = True
        op2 = run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf,
                      write_f06=write_f06,
                      debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
                      build_pandas=build_pandas)[0]
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    @unittest.skipIf(not IS_H5PY, "No h5py")
    def test_op2_solid_bending_02(self):
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2 = OP2(debug=False, log=None, debug_file=None, mode=None)
        op2.load_as_h5 = True
        op2.read_op2(op2_filename=op2_filename, combine=True,
                     build_dataframe=None, skip_undefined_matrices=False,
                     encoding=None)
        #op2 = read_op2(op2_filename, debug=False)
        del op2

    def test_op2_solid_bending_02_geom(self):
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        hdf5_filename = os.path.join(folder, 'solid_bending.h5')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False, write_hdf5=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False, build_pandas=True)
        if IS_PANDAS:
            assert op2.displacements[1].data_frame is not None

        op2.print_subcase_key()
        if IS_H5PY:
            op2.export_hdf5_filename(hdf5_filename)
            op2b = OP2(debug=False)
            op2b.load_hdf5_filename(hdf5_filename, combine=True)
            op2b.print_subcase_key()

    def test_op2_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'static_solid_shell_bar.test_op2.f06')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False, build_pandas=True)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_mode_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        subcases = [1]
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False,
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling_solid_shell_bar.op2')
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False,
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling2_solid_shell_bar.op2')
        unused_op2 = read_op2_geom(op2_filename, debug=False)
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases)
        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 0, len(op2.eigenvectors)
        str(op2.isubcase_name_map)
        str(op2.displacements[1].subtitle)
        str(op2.displacements[1].label)

        subcases = 2
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=True)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = 2
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = [1, 2]
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        assert len(op2.displacements) == 1, len(op2.displacements).v
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_transient_solid_shell_bar_01_geom(self):
        """transient test"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'transient_solid_shell_bar.test_op2.f06')
        build_pandas = IS_TRANSIENT_PANDAS
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=False, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False, build_pandas=build_pandas)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_frequency_solid_shell_bar_01_geom(self):
        """frequency test"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')
        unused_op2 = read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=False, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, check_memory=False, stop_on_failure=True,
            dev=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_transfer_function_01(self):
        """tests the transfer function cards work"""
        folder = os.path.join(MODEL_PATH, 'transfer_function')
        #bdf_filename = os.path.join(folder, 'actuator_tf_modeling.bdf')
        op2_filename = os.path.join(folder, 'actuator_tf_modeling.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')

        unused_op2 = read_op2_geom(op2_filename, debug=False)

        debug = False
        write_bdf = True
        write_f06 = True
        make_geom = True
        op2, unused_is_passed = run_op2(
            op2_filename, write_bdf=write_bdf, make_geom=make_geom,
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
        folder = os.path.join(MODEL_PATH, 'aero', 'monpnt3')
        op2_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.op2')
        f06_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.test_op2.f06')
        op2 = read_op2(op2_filename, debug=False, debug_file='temp.debug')
        monitor3 = op2.monitor3
        assert len(monitor3.frequencies) == 11, monitor3
        str(monitor3)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)
        os.remove('temp.debug')

    def test_op2_nastran_2005r3b(self):
        """Nastran 2005r3 bug"""
        folder = os.path.join(MODEL_PATH, 'modele_petite_zone')
        op2_filename = os.path.join(folder, 'modele_petite_zone.op2')
        f06_filename = os.path.join(folder, 'modele_petite_zone.test_op2.f06')
        op2 = read_op2_geom(op2_filename, debug=False)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_solid_shell_bar_01(self):
        """tests sol_101_elements/static_solid_shell_bar.op2"""
        op2_filename = os.path.join('static_solid_shell_bar.op2')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, op2_filename)
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
            load_as_h5=False)

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
        # TODO: should this be 4; yes by actual count...
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (1, 20, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_force = op2.ctria3_force[isubcase]
        assert ctria3_force.nelements == 8, ctria3_force.nelements
        assert ctria3_force.data.shape == (1, 8, 8), ctria3_force.data.shape

        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        ctetra_stress = op2.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 10, 10), ctetra_stress.data.shape
        assert ctetra_stress.is_von_mises, ctetra_stress

        cpenta_stress = op2.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (1, 14, 10), cpenta_stress.data.shape
        assert cpenta_stress.is_von_mises, cpenta_stress

        chexa_stress = op2.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (1, 9, 10), chexa_stress.data.shape
        assert chexa_stress.is_von_mises, chexa_stress

        if IS_PANDAS:
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'static_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)
        os.remove(vtk_filename)

    def test_op2_solid_shell_bar_01_straincurvature(self):
        """tests sol_101_elements/static_solid_shell_bar_straincurve.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        unused_bdf_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = op2.ctria3_strain[isubcase]
        sword = get_scode_word(ctria3_strain.s_code, ctria3_strain.stress_bits)
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance, '%s\n%s' % (ctria3_strain, sword)
        assert ctria3_strain.is_von_mises, '%s\n%s' % (ctria3_strain, sword)

        cquad4_strain = op2.cquad4_strain[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance, cquad4_strain
        assert cquad4_strain.is_von_mises, cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance(self):
        """tests sol_101_elements/static_solid_shell_bar_fiberdist.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = op2.ctria3_stress[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance, ctria3_strain
        assert ctria3_strain.is_von_mises, ctria3_strain

        cquad4_strain = op2.cquad4_stress[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), '%s\n%s' % (cquad4_strain.data.shape, sword)
        assert cquad4_strain.is_fiber_distance, '%s\n%s' % (cquad4_strain, sword)
        assert cquad4_strain.is_von_mises, '%s\n%s' % (cquad4_strain, sword)

    def test_op2_solid_shell_bar_01_straincurvature_shear(self):
        """tests sol_101_elements/static_solid_shell_bar_straincurve_shear.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve_shear.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_straincurve_shear.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert not ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert not cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = op2.ctria3_strain[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance, ctria3_strain
        assert not ctria3_strain.is_von_mises, ctria3_strain

        cquad4_strain = op2.cquad4_strain[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance, cquad4_strain
        assert not cquad4_strain.is_von_mises, cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance_shear(self):
        """tests sol_101_elements/static_solid_shell_bar_fiberdist_shear.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        #bdf_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist_shear.bdf')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar_fiberdist_shear.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcase = 1
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises is False, ctria3_stress

        cquad4_stress = op2.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises is False, cquad4_stress

        ctria3_strain = op2.ctria3_stress[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance, ctria3_strain
        assert ctria3_strain.is_von_mises is False, ctria3_strain

        cquad4_strain = op2.cquad4_stress[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert cquad4_strain.is_fiber_distance, cquad4_strain
        assert cquad4_strain.is_von_mises is False, cquad4_strain

    def test_op2_solid_shell_bar_mode(self):
        """tests sol_101_elements/mode_solid_shell_bar.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
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

        if IS_PANDAS:
            cbar_force.build_dataframe()
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_mode_export(self):
        """tests sol_101_elements/mode_solid_shell_bar.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        bdf_filename = os.path.join(folder, 'mode_solid_shell_bar.bdf')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        vtk_filename = os.path.join(folder, 'mode_solid_shell_bar.vtk')
        export_to_vtk_filename(bdf_filename, op2_filename, vtk_filename)
        os.remove(vtk_filename)

    def test_op2_solid_shell_bar_buckling(self):
        """tests sol_101_elements/buckling_solid_shell_bar.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True)

        isubcases = [(1, 1, 1, 0, 0, '', ''), (1, 8, 1, 0, 0, '', '')]
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

        if IS_PANDAS:
            cbar_force.build_dataframe()
            str(cbar_force.data_frame)

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_freq(self):
        """tests sol_101_elements/freq_solid_shell_bar.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, debug=False)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
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
        assert ctria3_stress.data.shape == (7, 16, 3), ctria3_stress.data.shape

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

        if IS_PANDAS:
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
        """tests sol_101_elements/freq_solid_shell_bar.op2"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
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
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)

        build_pandas = IS_TRANSIENT_PANDAS
        read_op2(op2_filename, debug=debug)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
            build_pandas=build_pandas)
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
        assert cquad4_stress.nelements == 40, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (21, 40, 8), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = op2.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 16, ctria3_stress.nelements  # TODO: 8-wrong
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

        if IS_TRANSIENT_PANDAS:
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

    def _test_op2_autodesk_1(self):
        """tests an Autodesk Nastran example"""
        op2_filename = os.path.join(PKG_PATH, 'op2', 'test', 'examples',
                                    'autodesk', 'aa8lzviq9.op2')
        log = get_logger(level='warning')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, write_f06=False,
            log=log, stop_on_failure=True, binary_debug=True, quiet=True,
            post=-4)

        assert len(op2.displacements) == 1
        assert len(op2.spc_forces) == 1
        assert len(op2.ctetra_stress) == 1

        isubcase = 1
        ctetra_stress = op2.ctetra_stress[isubcase]
        if IS_PANDAS:
            ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 810, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 810*5, 10), ctetra_stress.data.shape

        assert len(op2.cpenta_stress) == 0
        assert len(op2.chexa_stress) == 0
        assert len(op2.grid_point_forces) == 0

    def test_op2_optistruct_1(self):
        """
        Optistruct 2012 Tables : CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                                OUGV1, OES1X
        """
        op2_filename = os.path.join(MODEL_PATH, 'optistruct', 'hm14.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        log = get_logger(level='warning')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2_geom(op2_filename, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            log=log, stop_on_failure=True, binary_debug=True, quiet=True)
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
        if IS_PANDAS:
            ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 3951, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 19755, 10), ctetra_stress.data.shape

        assert len(op2.cpenta_stress) == 0
        assert len(op2.chexa_stress) == 0

        assert len(op2.grid_point_forces) == 0
        os.remove(debug_file)

    def test_op2_plate_py_01(self):
        """tests plate_py/plate_py.op2"""
        make_geom = False
        write_bdf = False
        write_f06 = False
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.op2')
        read_op2(op2_filename, log=log)
        run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf,
                write_f06=write_f06,
                log=log, stop_on_failure=True, quiet=True)

        make_geom = False
        write_bdf = False
        write_f06 = True
        run_op2(op2_filename, make_geom=make_geom, write_bdf=write_bdf,
                write_f06=write_f06,
                log=log, stop_on_failure=True, quiet=True)

        argv = ['test_op2', op2_filename, '-tgc', '--quiet', '--safe']
        test_op2(argv)

    def test_op2_good_sine_01(self):
        """tests freq_sine/good_sine.op2"""
        op2_filename = os.path.join(MODEL_PATH, 'freq_sine', 'good_sine.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        log = get_logger(level='warning')
        read_op2(op2_filename, log=log)
        build_pandas = False # IS_TRANSIENT_PANDAS #or PY3
        op2i, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            log=log, stop_on_failure=True,
            quiet=True, build_pandas=build_pandas)

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

        #if IS_PANDAS and not PY3:
            #acc.build_dataframe()
        unused_accx = acc.extract_xyplot(nids, 1, 'real')
        unused_accxi = acc.extract_xyplot(nids, 1, 'imag')
        #print(accx)
        #print(accxi)
        #make_geom = False
        #write_bdf = False
        #write_f06 = True
        #run_op2(op2file, make_geom=make_geom, write_bdf=write_bdf,
                #write_f06=write_f06,
                #debug=debug, stopOnFailure=True)

    def test_op2_good_sine_02(self):
        """tests freq_sine/good_sine.op2"""
        folder = os.path.abspath(os.path.join(MODEL_PATH))
        bdf_filename = os.path.join(folder, 'freq_sine', 'good_sine.dat')
        op2_filename = os.path.join(folder, 'freq_sine', 'good_sine.op2')
        #make_geom = False
        #write_bdf = False
        #write_f06 = True
        log = get_logger(level='warning')
        bdf = BDF(log=log)
        bdf.read_bdf(bdf_filename)

        debug = False
        debug_file = 'debug.out'

        read_op2(op2_filename, log=log)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')

        _verify_ids(bdf, op2, isubcase=1)
        os.remove(debug_file)

    def test_op2_bcell_01(self):
        """tests other/bcell9p0.op2"""
        folder = os.path.abspath(os.path.join(MODEL_PATH))
        bdf_filename = os.path.join(folder, 'other', 'bcell9p0.bdf')
        op2_filename = os.path.join(folder, 'other', 'bcell9p0.op2')
        #make_geom = False
        #write_bdf = False
        #write_f06 = True
        log = get_logger(level='warning')
        bdf = BDF(log=log)
        bdf.read_bdf(bdf_filename, xref=False)

        debug = False
        debug_file = 'debug.out'
        read_op2(op2_filename, debug=False)
        op2 = OP2(debug=debug, debug_file=debug_file)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')
        os.remove(debug_file)

        _verify_ids(bdf, op2, isubcase=1)

        msg = ''
        for isubcase, keys in sorted(op2.subcase_key.items()):
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
        """tests cbush/cbush.op2"""
        op2_filename = os.path.join(MODEL_PATH, 'cbush', 'cbush.op2')
        make_geom = False
        write_bdf = False
        write_f06 = True
        log = get_logger(level='warning')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06, compare=True,
            log=log, stop_on_failure=True, binary_debug=True, quiet=True)
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
        if IS_PANDAS:
            op2.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(os.path.dirname(op2_filename))
        os.remove(debug_file)

    @unittest.expectedFailure
    def test_set_times_01(self):
        """specify the modes to extract"""
        model = OP2(debug=False, log=None, debug_file=None, mode='msc')

        # subcase : times
        times = {1 : [2, 3]}
        model.set_transient_times(times)
        op2_filename = os.path.abspath(os.path.join(MODEL_PATH, 'sol_101_elements',
                                                    'mode_solid_shell_bar.op2'))
        model.read_op2(op2_filename, combine=True, build_dataframe=None,
                       skip_undefined_matrices=False,
                       encoding=None)
        isubcase = 1
        #print(model.get_op2_stats(short=False))
        eigenvector = model.eigenvectors[isubcase]
        #print(eigenvector)
        assert len(eigenvector.modes) == 2, eigenvector.modes

    def test_random_ctria3(self):
        """runs a random test"""
        folder = os.path.join(MODEL_PATH, 'random')
        op2_filename = os.path.join(folder, 'random_test_bar_plus_tri.op2')
        f06_filename = os.path.join(folder, 'random_test_bar_plus_tri.test_op2.f06')
        op2 = read_op2_geom(op2_filename, debug=False)

        op2res = op2.op2_results
        assert len(op2res.psd.displacements) == 1
        assert len(op2res.rms.displacements) == 1
        assert len(op2res.crm.displacements) == 1
        assert len(op2res.no.displacements) == 1
        assert len(op2res.psd.accelerations) == 1
        assert len(op2res.rms.accelerations) == 1
        assert len(op2res.crm.accelerations) == 1
        assert len(op2res.no.accelerations) == 1
        assert len(op2res.crm.cbar_force) == 1
        assert len(op2res.psd.cbar_force) == 1
        assert len(op2res.rms.cbar_force) == 1
        assert len(op2res.no.cbar_force) == 1
        assert len(op2res.crm.cquad4_force) == 1
        assert len(op2res.psd.cquad4_force) == 1
        assert len(op2res.rms.cquad4_force) == 1
        assert len(op2res.no.cquad4_force) == 1
        assert len(op2res.crm.ctria3_force) == 1
        assert len(op2res.psd.ctria3_force) == 1
        assert len(op2res.rms.ctria3_force) == 1
        assert len(op2res.no.ctria3_force) == 1
        assert len(op2res.no.cbar_force) == 1
        assert len(op2res.no.cbar_force) == 1
        assert len(op2res.no.cbar_force) == 1
        assert len(op2.eigenvalues) == 1
        assert 'BHH' in op2.matrices
        assert 'KHH' in op2.matrices

        #psd.displacements_psd[1]
        #rms.displacements[1]
        #crm.displacements[1]
        #no.displacements[1]
        #psd.accelerations[1]
        #rms.accelerations[1]
        #crm.accelerations[1]
        #no.accelerations[1]
        #crm.cbar_force[1]
        #psd.cbar_force[1]
        #rms.cbar_force[1]
        #no.cbar_force[1]
        #eigenvalues[u'RANDOM TEST']
        #crm.cquad4_force[1]
        #psd.cquad4_force[1]
        #rms.cquad4_force[1]
        #no.cquad4_force[1]
        #crm.ctria3_force[1]
        #psd.ctria3_force[1]
        #rms.ctria3_force[1]
        #no.ctria3_force[1]
        #Matrix['BHH'];   shape=(20, 20);
        # type=numpy.matrixlib.defmatrix.matrix; dtype=float64; desc=symmetric
        #Matrix['KHH'];   shape=(20, 20);
        # type=numpy.matrixlib.defmatrix.matrix; dtype=float64; desc=symmetric

        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_random_ctria3_oesrmx1(self):
        """runs a random test"""
        folder = os.path.join(MODEL_PATH, 'random')
        op2_filename = os.path.join(folder, 'rms_tri_oesrmx1.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        unused_op2 = read_op2_geom(op2_filename, debug=False)

    def test_aero_cpmopt(self):
        """test optimization with multiple GEOM1 tables and multiple CSTM tables"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'aero')
        op2_filename = os.path.join(folder, 'cpmopt.op2')
        bdf_filename = os.path.join(folder, 'cpmopt.bdf')
        #unused_op2 = read_op2_geom(op2_filename, debug=False)
        read_bdf(bdf_filename, debug=False)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=None, write_f06=True,
            write_op2=False, write_hdf5=False, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=False, build_pandas=False, subcases=None,
            exclude=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, check_memory=False, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_ogs(self):
        """test grid_point_stresses"""
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'ogs', 'ogs.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=False, write_hdf5=True, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=False, build_pandas=True, subcases=None,
            exclude=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, check_memory=False, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_spcadd(self):
        """tests loading SPCADD/MPCADDs"""
        model = BDF()
        model.is_debug_file = False
        datai = np.array([2, 1, 10, -1], dtype='int32')
        _read_spcadd_mpcadd(model, 'SPCADD', datai)

        datai = np.array([3, 1, -1], dtype='int32')
        _read_spcadd_mpcadd(model, 'SPCADD', datai)


        datai = np.array([4, 1, 10, -1], dtype='int32')
        _read_spcadd_mpcadd(model, 'MPCADD', datai)

        datai = np.array([5, 1, -1], dtype='int32')
        _read_spcadd_mpcadd(model, 'MPCADD', datai)
        assert len(model.spcadds) == 2, model.spcadds
        assert len(model.mpcadds) == 2, model.mpcadds

def _verify_ids(bdf, op2, isubcase=1):
    """helper function for tests"""
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

if __name__ == '__main__':  # pragma: no cover
    ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
    if not ON_RTD:
        unittest.main()
