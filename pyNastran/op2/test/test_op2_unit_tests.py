"""various OP2 tests"""
import os
import unittest
from pathlib import Path

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
except ModuleNotFoundError:
    IS_PANDAS = False

try:
    import h5py  # pylint: disable=unused-import
    IS_H5PY = True
except ModuleNotFoundError:  # pragma: no cover
    IS_H5PY = False


import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, CORD2R
from pyNastran.op2.op2 import OP2, read_op2, FatalError, FortranMarkerError
from pyNastran.op2.op2_interface.op2_common import get_scode_word
from pyNastran.op2.op2_geom import OP2Geom, read_op2_geom
from pyNastran.op2.test.test_op2 import run_op2, main as test_op2

from pyNastran.bdf.test.test_bdf_unit_tests import Tester
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.bdf_interface.compare_card_content import compare_elements

#from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    #RealPlateBilinearForceArray, RealPlateForceArray)
#from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
from pyNastran.op2.vector_utils import filter1d, abs_max_min_global, abs_max_min_vector
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray
from pyNastran.femutils.test.utils import is_array_close
from pyNastran.op2.result_objects.grid_point_weight import make_grid_point_weight
from pyNastran.op2.tables.geom.geom4 import _read_spcadd_mpcadd
from pyNastran.f06.csv_writer import write_csv

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..'/ 'models').resolve()
OP2_TEST_PATH = (PKG_PATH / 'op2' / 'test' / 'examples').resolve()


class TestOP2Unit(Tester):
    """various OP2 tests"""
    def test_grid_point_weight(self):
        """tests GridPointWeight"""
        reference_point = 0
        MO = np.array([
            [9.25, -.933, 0, 0, 2., -2.967],
            [-.933, 10.75, 0., -3., 0., 6.058],
            [0., 0., 20., 5., -10., 0.],
            [0., -3., 5., 8., -2.5, -1.5],
            [2., 0., -10., -2.5, 9.5, 0.],
            [-2.967, 6.058, 0., -1.5, 0., 7.496],
        ])
        weight = make_grid_point_weight(
            reference_point, MO, approach_code=1, table_code=13,
            title='', subtitle='', label='',
            superelement_adaptivity_index='')
        str(weight)


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


class TestOP2Functions(unittest.TestCase):
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


class TestAutodeskOP2(Tester):
    """various OP2 tests"""
    def _test_op2_autodesk_1(self):
        """tests an Autodesk Nastran example"""
        op2_filename = os.path.join(PKG_PATH, 'op2', 'test', 'examples',
                                    'autodesk', 'aa8lzviq9.op2')
        log = get_logger(level='warning')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, write_f06=False,
            is_autodesk=True, log=log, stop_on_failure=True, binary_debug=True, quiet=True,
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

    def test_op2_autodesk_2(self):
        """tests an Autodesk Nastran example"""
        op2_filename = os.path.join(MODEL_PATH, 'autodesk', '9zk6b5uuo.op2')
        log = get_logger(level='warning')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, write_f06=False,
            log=log, stop_on_failure=True, binary_debug=True, quiet=True,
            is_autodesk=True, post=-4)

        stress = op2.op2_results.stress
        strain_energy = op2.op2_results.strain_energy
        assert len(op2.displacements) == 4, len(op2.displacements)
        assert len(op2.spc_forces) == 4, len(op2.spc_forces)
        assert len(stress.ctetra_stress) == 4, len(stress.ctetra_stress)
        assert len(strain_energy.ctetra_strain_energy) == 4, len(strain_energy.ctetra_strain_energy)

        isubcase = 1
        ctetra_stress = stress.ctetra_stress[isubcase]
        if IS_PANDAS:
            ctetra_stress.build_dataframe()
        nelements = 36
        assert ctetra_stress.nelements == nelements, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, nelements*5, 10), ctetra_stress.data.shape

        assert len(stress.cpenta_stress) == 0
        assert len(stress.chexa_stress) == 0
        assert len(op2.grid_point_forces) == 0

    def test_op2_autodesk_3(self):
        """tests an Autodesk Nastran example"""
        op2_filename = os.path.join(MODEL_PATH, 'autodesk', 'nonlinear_beam.op2')
        log = get_logger(level='warning')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, write_f06=True,
            log=log, stop_on_failure=True, binary_debug=True, quiet=True,
            is_autodesk=True, post=-4)

        stress = op2.op2_results.stress
        strain_energy = op2.op2_results.strain_energy
        assert len(op2.displacements) == 4, len(op2.displacements)
        assert len(op2.spc_forces) == 4, len(op2.spc_forces)
        assert len(stress.ctetra_stress) == 4, len(stress.ctetra_stress)
        assert len(op2.nonlinear_ctetra_stress_strain) == 4, len(op2.nonlinear_ctetra_stress_strain)

        #isubcase = 1
        #ctetra_stress = op2.ctetra_stress[isubcase]
        #if IS_PANDAS:
            #ctetra_stress.build_dataframe()
        #nelements = 36
        #assert ctetra_stress.nelements == nelements, ctetra_stress.nelements
        #assert ctetra_stress.data.shape == (1, nelements*5, 10), ctetra_stress.data.shape

        assert len(strain_energy.ctetra_strain_energy) == 0, len(strain_energy.ctetra_strain_energy)
        assert len(stress.cpenta_stress) == 0, len(stress.cpenta_stress)
        assert len(stress.chexa_stress) == 0, len(stress.chexa_stress)
        assert len(op2.grid_point_forces) == 4, len(op2.grid_point_forces)

class TestOptistructOP2(Tester):
    """various OP2 tests"""
    def test_op2_optistruct_1(self):
        """
        Optistruct 2012 Tables : CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                                OUGV1, OES1X
        """
        op2_filename = os.path.join(MODEL_PATH, 'optistruct', 'hm14.op2')
        make_geom = True
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

        stress = op2.op2_results.stress
        strain = op2.op2_results.strain
        force = op2.op2_results.force

        # rod_force = force.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        #assert len(op2.rod_force) == 0
        assert len(stress.crod_stress) == 0

        assert len(force.cbar_force) == 0
        assert len(stress.cbar_stress) == 0
        assert len(stress.cbeam_stress) == 0

        assert len(stress.cquad4_stress) == 0
        assert len(stress.ctria3_stress) == 0

        ctetra_stress = stress.ctetra_stress[isubcase]
        if IS_PANDAS:
            ctetra_stress.build_dataframe()
        assert ctetra_stress.nelements == 3951, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 19755, 10), ctetra_stress.data.shape

        assert len(stress.cpenta_stress) == 0
        assert len(stress.chexa_stress) == 0

        assert len(op2.grid_point_forces) == 0
        os.remove(debug_file)


class TestNASA95OP2(Tester):
    """various OP2 tests"""
    def test_op2_nasa_nastran_01(self):
        """checks sdr11se_s2dc.bdf, which tests ComplexCBushStressArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'nasa_nastran', 'balsa_wingbox.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nasa_nastran', 'balsa_wingbox.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, is_nasa95=True, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)


class TestSATKOP2(Tester):
    """various OP2 tests"""
    def test_bdf_op2_satk_1(self):
        """checks pn_mwe_s-solution_1.dat, which tests dynamics matrices"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'satk', 'pn_mwe_s-solution_1.dat')
        op2_filename = os.path.join(MODEL_PATH, 'satk', 'pn_mwe_s-solution_1.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_satk_2(self):
        """
        checks pn_mwe_s-response_dynamics_1-random_base_excitation_1.op2,
        which tests OUGPK1, OEFPK1 RMS results
        """
        log = get_logger(level='info')
        #bdf_filename = os.path.join(MODEL_PATH, 'satk', 'pn_mwe_s-solution_1.dat')
        op2_filename = os.path.join(MODEL_PATH, 'satk',
                                    'pn_mwe_s-response_dynamics_1-random_base_excitation_1.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=False, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        assert is_passed, is_passed
        str(op2.get_op2_stats())

    def test_bdf_op2_satk_3(self):
        """checks pn_mwe_s-sol_111.dat, which tests PSD tables"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'cbush_psd_bug', 'pn_mwe_s-sol_111.dat')
        op2_filename = os.path.join(MODEL_PATH, 'cbush_psd_bug', 'pn_mwe_s-sol_111.op2')
        op2_filename2 = os.path.join(MODEL_PATH, 'cbush_psd_bug', 'pn_mwe_s-rd_satk-rbe.op2')

        #  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        save_load_deck(model, run_save_load=True, run_renumber=False, run_save_load_hdf5=False)

        log = get_logger(level='warning')
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=False, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        assert is_passed, is_passed
        #op2_results.psds.accelerations:
            # (subtitle, analysis_code, stress_strain_flag, node, dof)
            #('SUBCASE - MODAL FREQUENCY 1', 9, 4)

        #op2_results.psds.force:
            # (subtitle, analysis_code, stress_strain_flag, node, dof)
            #('SUBCASE - MODAL FREQUENCY 1', 2, 2)
            #('SUBCASE - MODAL FREQUENCY 1', 8, 2)

        #op2_results.eqexin: EQEXIN(nid, ndof, doftype); nnodes=7
        #op2_results.gpdt: GPDT(nid_cp_cd_ps, xyz); nnodes=7
        #op2_results.bgpdt: BGPDT(cd, xyz); nnodes=7
        assert len(op2.op2_results.force.cbush_force[1].modes) == 8
            #type=RealCBushForceArray ntimes=8 nelements=1
            #data: [ntimes, nnodes, 6] where 6=[fx, fy, fz, mx, my, mz]
            #data.shape = (8, 1, 6)
            #element.shape = (1,)
            #element name: CBUSH
            #sort1
            #modes = [1 2 3 4 5 6 7 8]; dtype=int32
            #eigns = [1.607e+05 1.607e+05 2.224e+07 1.204e+09 1.204e+09 1.102e+10 1.102e+10
                   #2.544e+10]; dtype=float64
            #cycles = [   63.804    63.804   750.561  5522.769  5522.769 16707.479 16707.479
                       #25385.576]; dtype=float64

        assert len(op2.op2_results.force.cbar_force[1].modes) == 8
            #type=RealCBarForceArray ntimes=8 nelements=5; table_name='OEF1'
            #data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
            #data.shape = (8, 5, 8)
            #element.shape = (5,)
            #element name: CBAR-34
            #sort1
            #modes = [1 2 3 4 5 6 7 8]; dtype=int32
            #eigns = [1.607e+05 1.607e+05 2.224e+07 1.204e+09 1.204e+09 1.102e+10 1.102e+10
                   #2.544e+10]; dtype=float64
            #cycles = [   63.804    63.804   750.561  5522.769  5522.769 16707.479 16707.479
                       #25385.576]; dtype=float64

        assert len(op2.op2_results.psd.accelerations[(3, 5, 2, 0, 0, '', 'RANDOM  103')].freqs) == 127
            #isubcase = 3
            #type=RealAccelerationArray ntimes=127 nnodes=7, table_name=OUGPSD1
            #data: [t1, t2, t3, r1, r2, r3] shape=[127, 7, 6] dtype=float32
            #node_gridtype.shape = (7, 2)
            #sort1
            #freqs = [  20.      20.943   21.93  ... 1861.173 1909.985 2000.   ]; dtype=float32

        assert len(op2.op2_results.psd.cbar_force[(3, 5, 2, 0, 0, '', 'RANDOM  103')].freqs) == 127
            #type=RealCBarForceArray ntimes=127 nelements=5; table_name='OEFPSD1'
            #data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
            #data.shape = (127, 10, 8)
            #element.shape = (5,)
            #element name: CBAR-34
            #sort1
            #freqs = [  20.      20.943   21.93  ... 1861.173 1909.985 2000.   ]; dtype=float32

        assert len(op2.op2_results.psd.cbush_force[(3, 5, 2, 0, 0, '', 'RANDOM  103')].freqs) == 127
            #type=RealCBushForceArray ntimes=1 nelements=127
            #data: [ntimes, nnodes, 6] where 6=[fx, fy, fz, mx, my, mz]
            #data.shape = (127, 1, 6)
            #element.shape = (1,)
            #element name: CBUSH
            #sort1
            #freqs = [  20.      20.943   21.93  ... 1861.173 1909.985 2000.   ]; dtype=float32

        assert len(op2.op2_results.rms.accelerations[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #isubcase = 3
            #type=RealAccelerationArray ntimes=1 nnodes=7, table_name=OUGRMS1
            #data: [t1, t2, t3, r1, r2, r3] shape=[1, 7, 6] dtype=float32
            #node_gridtype.shape = (7, 2)
            #sort1
            #freqs = [0.]; dtype=float32

        assert len(op2.op2_results.rms.cbar_force[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #type=RealCBarForceArray ntimes=1 nelements=5; table_name='OEFRMS1'
            #data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
            #data.shape = (1, 5, 8)
            #element.shape = (5,)
            #element name: CBAR-34
            #sort1
            #freqs = [0.]; dtype=float32

        assert len(op2.op2_results.rms.cbush_force[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #type=RealCBushForceArray ntimes=1 nelements=1
            #data: [ntimes, nnodes, 6] where 6=[fx, fy, fz, mx, my, mz]
            #data.shape = (1, 1, 6)
            #element.shape = (1,)
            #element name: CBUSH
            #sort1
            #freqs = [0.]; dtype=float32

        assert len(op2.op2_results.no.accelerations[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #isubcase = 3
            #type=RealAccelerationArray ntimes=1 nnodes=7, table_name=OUGNO1
            #data: [t1, t2, t3, r1, r2, r3] shape=[1, 7, 6] dtype=float32
            #node_gridtype.shape = (7, 2)
            #sort1
            #freqs = [0.]; dtype=float32

        assert len(op2.op2_results.no.cbar_force[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #type=RealCBarForceArray nelements=5; table_name='OEFNO1'
            #data: [1, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
            #data.shape = (1, 5, 8)
            #element.shape = (5,)
            #element name: CBAR-34
            #sort1
            #freqs = [0.]; dtype=float32

        assert len(op2.op2_results.no.cbush_force[(3, 5, 1, 0, 0, '', 'RANDOM  103')].freqs) == 1
            #type=RealCBushForceArray nelements=1
            #data: [1, nnodes, 6] where 6=[fx, fy, fz, mx, my, mz]
            #data.shape = (1, 1, 6)
            #element.shape = (1,)
            #element name: CBUSH
            #sort1
            #freqs = [0.]; dtype=float32

        op2.eigenvalues['']
            #type=RealEigenvalues neigenvalues=8
            #title, extraction_order, eigenvalues, radians, cycles, generalized_mass, generalized_stiffness

        op2.matrices['BHH']
        op2.matrices['KHH']
            #Matrix['BHH'];      shape=(8, 8);     type=scipy.sparse._coo.coo_matrix;     dtype=float64;   desc=symmetric
            #Matrix['KHH'];      shape=(8, 8);     type=scipy.sparse._coo.coo_matrix;     dtype=complex128; desc=symmetric

        str(op2.get_op2_stats())

        op2, is_passed = run_op2(
            op2_filename2, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=False, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        assert is_passed, is_passed
        #print(op2.get_op2_stats())
        #op2_results.gpdt: GPDT(nid_cp_cd_ps, xyz); nnodes=7

        assert len(op2.op2_results.psd.displacements[(1, 5, 2, 0, 0, '', '')].freqs) == 298
          #isubcase = 1
          #type=RealDisplacementArray ntimes=298 nnodes=6, table_name=OUGPSD1
          #data: [t1, t2, t3, r1, r2, r3] shape=[298, 6, 6] dtype=float32
          #node_gridtype.shape = (6, 2)
          #sort1
          #freqs = [  20.      20.222   20.447 ... 1987.249 1993.657 2000.   ]; dtype=float32

        assert len(op2.op2_results.psd.accelerations[(1, 5, 2, 0, 0, '', '')].freqs) == 298
          #isubcase = 1
          #type=RealAccelerationArray ntimes=298 nnodes=6, table_name=OUGPSD1
          #data: [t1, t2, t3, r1, r2, r3] shape=[298, 6, 6] dtype=float32
          #node_gridtype.shape = (6, 2)
          #sort1
          #freqs = [  20.      20.222   20.447 ... 1987.249 1993.657 2000.   ]; dtype=float32

        assert len(op2.op2_results.psd.cbar_force[(1, 5, 2, 0, 0, '', '')].freqs) == 298
          #type=RealCBarForceArray ntimes=298 nelements=5; table_name='OEFPSD1'
          #data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
          #data.shape = (298, 10, 8)
          #element.shape = (5,)
          #element name: CBAR-34
          #sort1
          #freqs = [  20.      20.222   20.447 ... 1987.249 1993.657 2000.   ]; dtype=float32

        assert len(op2.op2_results.psd.cbush_force[(1, 5, 2, 0, 0, '', '')].freqs) == 298
          #type=RealCBushForceArray ntimes=1 nelements=298
          #data: [ntimes, nnodes, 6] where 6=[fx, fy, fz, mx, my, mz]
          #data.shape = (298, 1, 6)
          #element.shape = (1,)
          #element name: CBUSH
          #sort1
          #freqs = [  20.      20.222   20.447 ... 1987.249 1993.657 2000.   ]; dtype=float32

        assert len(op2.op2_results.rms.displacements[(1, 5, 1, 0, 0, '', '')].freqs) == 1
        assert len(op2.op2_results.rms.velocities[(1, 5, 1, 0, 0, '', '')].freqs) == 1
        assert len(op2.op2_results.rms.accelerations[(1, 5, 1, 0, 0, '', '')].freqs) == 1
        assert len(op2.op2_results.rms.cbar_force[(1, 5, 1, 0, 0, '', '')].freqs) == 1
        assert len(op2.op2_results.rms.cbush_force[(1, 5, 1, 0, 0, '', '')].freqs) == 1


class TestNX(Tester):
    def test_nx_flutter(self):
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'aero', 'flutter_bug', 'wing_b1.op2')
        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_nx_glue_slide_distance(self):
        """test NX 2020 version"""
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'glue', 'n401gsh01.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_nx_bolt(self):
        """test NX 2020 version"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'nx', 'bolt', 's402bolt304.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'bolt', 's402bolt304.op2')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_nx_initial_final_separation(self):
        """
        checks nx/contact_model.bdf, which tests
        initial/final contact separation distance
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'nx', 'contact_model.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'contact_model.op2')

        #  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_nx_composite_solids(self):
        """
        checks nx/composite_solids/test.bdf, which tests
        centroidal CHEXA composite stress
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'nx', 'composite_solids', 'test.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'composite_solids', 'test.op2')

        #  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_nx_composite_solids_corner(self):
        """
        checks nx/composite_solids/test_nx_corner.bdf, which tests
        corner CHEXA composite stress
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'nx', 'composite_solids', 'test_nx_corner.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'composite_solids', 'test_nx_corner.op2')

        #  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)


class TestMSC(Tester):
    def test_msc_2014(self):
        """test MSC 2014 version"""
        log = get_logger(level='warning')
        op2_filename1 = os.path.join(MODEL_PATH, 'bugs', 'msc_2014', 'sdof_crod_2014.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename1, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_msc_cfast(self):
        """test MSC 126-CFAST"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'msc', 'test_model_cfast.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'msc', 'test_model_cfast.op2')
        model = read_bdf(bdf_filename, encoding='ascii', debug=False, log=log)
        bdf_filename_out = os.path.join(MODEL_PATH, 'msc', 'test_model_cfast_out.bdf')
        model.write_bdf(bdf_filename_out)
        os.remove(bdf_filename_out)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, read_bdf=None, write_f06=False,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=True)

    def test_msc_dscmcol(self):
        """test MSC 126 DSCMCOL-matrix sensitivites"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_dscmcol', 'goland_final_test.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_dscmcol', 'goland_final_test.op2')
        model = read_bdf(bdf_filename, encoding='ascii', debug=False, log=log)
        bdf_filename_out = os.path.join(MODEL_PATH, 'bugs', 'msc_dscmcol', 'test_goland_final_test.bdf')
        model.write_bdf(bdf_filename_out)
        os.remove(bdf_filename_out)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=True)

    def test_msc_2017_failure_indices_strength_ratio(self):
        """
        checks msc/failure_indices_strength_ratio/TestStressTemp.op2, which tests
         - op2_results.strength_ratio.cquad4_composite_stress[1]
         - op2_results.failure_indices.cquad4_composite_force[1]

        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'msc', 'failure_indices_strength_ratio', 'TestStressTemp.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'msc', 'failure_indices_strength_ratio', 'TestStressTemp.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=True)
        #model.safe_cross_reference()

        save_load_deck(model, run_save_load=True)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False, write_hdf5=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_msc_2017_units(self):
        """
        checks msc/units_mass_spring_damper/units_mass_spring_damper.op2, which tests
         - UNITS table for MSC 2014
        """
        log = get_logger(level='info')
        #bdf_filename = os.path.join(MODEL_PATH, 'msc', 'units_mass_spring_damper', 'test_nx_corner.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'msc', 'units_mass_spring_damper', 'units_mass_spring_damper.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        with self.assertRaises(NotImplementedError):
            run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                    write_f06=True, write_op2=False,
                    is_mag_phase=False,
                    is_sort2=False, is_nx=None, delete_f06=True,
                    subcases=None, exclude_results=None, short_stats=False,
                    compare=False, debug=False, binary_debug=True,
                    quiet=True,
                    stop_on_failure=True, dev=False,
                    build_pandas=True, log=log)

    def test_msc_2021_cbush_rbe3(self):
        """
        checks msc/cbush_2021/cbush_test.op2, which tests
         - UNITS table for MSC 2021
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'msc', 'cbush_2021', 'cbush_test.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'msc', 'cbush_2021', 'cbush_test.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=True)
        #model.safe_cross_reference()

        save_load_deck(model, run_save_load=True)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_msc_2020_rbe2(self):
        """
        checks bugs/msc_RBE_tests/rigid_rbe2--v2020.op2, which tests
         - RBE2 alpha for MSC 2020
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_RBE_tests', 'rigid_rbe2.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_RBE_tests', 'rigid_rbe2--v2020.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=True)
        #model.safe_cross_reference()

        save_load_deck(model, run_save_load=True)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_msc_2021_rbe2(self):
        """
        checks bugs/msc_RBE_tests/rigid_rbe2--v2021.1.op2, which tests
         - RBE2 tref for MSC 2021
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_RBE_tests', 'rigid_rbe2.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'bugs', 'msc_RBE_tests', 'rigid_rbe2--v2021.1.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=True)
        #model.safe_cross_reference()

        save_load_deck(model, run_save_load=True)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_op2_nastran_2005r3b(self):
        """Nastran 2005r3 bug"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'modele_petite_zone')
        op2_filename = os.path.join(folder, 'modele_petite_zone.op2')
        f06_filename = os.path.join(folder, 'modele_petite_zone.test_op2.f06')
        op2 = read_op2_geom(op2_filename, debug=False, log=log)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)


class TestOP2Main(Tester):
    """various OP2 tests"""
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_generalized_tables(self):
        """tests that set_additional_generalized_tables_to_read overwrites the GEOM1S class"""
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        model = OP2Geom(log=log)
        model.read_op2(op2_filename=op2_filename, combine=True, build_dataframe=False,
                       skip_undefined_matrices=False,
                       encoding=None)

        def read_some_table(self):
            """crashes"""
            raise NotImplementedError('read_some_table')

        model2 = OP2Geom(log=log)
        tables = {
            b'GEOM1S' : read_some_table,
        }
        model2.set_additional_generalized_tables_to_read(tables)
        with self.assertRaises(NotImplementedError):
            model2.read_op2(op2_filename=op2_filename, combine=True, build_dataframe=False,
                            skip_undefined_matrices=False,
                            encoding=None)

    def test_ibulk(self):
        """this test will fail if IBULK talble doesn't work"""
        log = get_logger(level='warning')
        unused_bdf_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        f06_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.test_op2.f06'))
        op2_filename = os.path.abspath(os.path.join(
            PKG_PATH, 'op2', 'test', 'examples', 'ibulk', 'model1_sim1-solution_1.op2'))
        op2 = OP2Geom(log=log, debug=False, debug_file='temp.debug')
        op2.read_op2(op2_filename)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)
        os.remove('temp.debug')

    def test_beam_modes(self):
        """tests the eigenvalue table reading
        nfreq = 10
        CBAR; n=9
        CBEAM; n=1
        """
        log = get_logger(level='warning')
        dirname = os.path.abspath(os.path.join(
            MODEL_PATH, 'beam_modes'))
        f06_filename = os.path.join(dirname, 'model1_sim1-solution_1.test_op2.f06')
        op2_filename_m1 = os.path.join(dirname, 'beam_modes_m1.op2')
        op2_filename_m2 = os.path.join(dirname, 'beam_modes_m2.op2')

        op2_filename_m1_out = os.path.join(dirname, 'beam_modes_m1_out.op2')
        op2_filename_m2_out = os.path.join(dirname, 'beam_modes_m2_out.op2')
        op2_1 = read_op2(op2_filename_m1, debug=False, log=log)
        op2_2 = OP2Geom(log=log, debug=False, debug_file='temp.debug')
        op2_2.read_op2(op2_filename_m2)
        op2_1.write_f06(f06_filename)

        op2_1.write_op2(op2_filename_m1_out, skips=['grid_point_weight']) #, is_mag_phase=False)
        op2_2.write_op2(op2_filename_m2_out, skips=['grid_point_weight']) #, is_mag_phase=False)
        os.remove(f06_filename)
        os.remove('temp.debug')
        os.remove(op2_filename_m1_out)
        os.remove(op2_filename_m2_out)

    def test_bdf_op2_elements_01(self):
        """tests a large number of elements and results in SOL 101"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.csv')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
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
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        with self.assertRaises(NotImplementedError):
            op2.save()
        write_csv(op2, csv_filename)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        write_csv(op2, csv_filename)

    def test_bdf_op2_elements_02(self):
        """tests a large number of elements and results in SOL 103-modes"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
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
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        write_csv(op2, csv_filename)


    def test_bdf_op2_elements_03(self):
        """tests a large number of elements and results in SOL 108-freq"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_04(self):
        """tests a large number of elements and results in SOL 108-freq"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.csv')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.test_op2.f06')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'freq_elements2.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, build_pandas=True, log=log)

        #print(op2.get_op2_stats())
        write_csv(op2, csv_filename)

        exclude_results = ['modal_contribution.cquadr_composite_strain', 'grid_point_forces',
                           'force*', 'load_vectors']
        op2b, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
            write_f06=False, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=exclude_results, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, build_pandas=True, log=log)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_05(self):
        """tests a large number of elements and results in SOL 106-loadstep"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'loadstep_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, build_pandas=True, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_elements_06(self):
        """tests a large number of elements and results in SOL 107-complex modes"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        #csv_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        #write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def _test_op2_shock_01(self):  # pragma: no cover
        """tests a large number of elements and results in SOL 103-shock analysis"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'shock', 'shock_analysis.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'shock', 'shock_analysis.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'shock', 'shock_analysis.csv')
        op2_filename = os.path.join(MODEL_PATH, 'shock', 'shock_analysis.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=True,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_post_minus4(self):
        """tests a large number of elements and results in SOL 107-complex modes"""
        log = get_logger(level='warning')
        unused_bdf_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements_post4.op2')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements_post4.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements_post4.op2')
        #fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #read_op2(op2_filename=op2_filename, combine=True, subcases=None,
                 #exclude_results=None, include_results=None,
                 #log=None, debug=True, debug_file=None,
                 #build_dataframe=False,
                 #skip_undefined_matrices=True, mode='msc',
                 #encoding=None)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, post=-4, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_01(self):
        """checks time_thermal_elements.bdf"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_02(self):
        """checks hd15306.bdf"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        #op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')

        bdf_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.csv')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'hd15306.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        write_csv(op2, csv_filename)
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_thermal_03(self):
        """checks time_thermal_elements.bdf"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.csv')
        op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        write_csv(op2, csv_filename)

    def test_bdf_op2_thermal_04(self):
        """checks time_thermal_elements.bdf"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'thermal', 'thermal_elements2.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'thermal', 'thermal_elements2.csv')
        op2_filename = os.path.join(MODEL_PATH, 'thermal', 'thermal_elements2.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        save_load_deck(model)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        write_csv(op2, csv_filename)

    def test_bdf_op2_thermal_05(self):
        """checks htflw47.bdf"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'thermal', 'htflw47.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'thermal', 'htflw47.csv')
        op2_filename = os.path.join(MODEL_PATH, 'thermal', 'htflw47.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log)
        #save_load_deck(model)

        # make_geom=False: duplicate ids
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=True, log=log)
        write_csv(op2, csv_filename)


    def test_cbar100(self):
        """tests a CBAR-100 model"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_bar_100.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'pbarl_bar_100.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_01(self):
        """checks ofprand1.bdf which tests nonlinear elements"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        #f06_filename = os.path.join(MODEL_PATH, 'elements', 'modes_complex_elements.test_op2.f06')
        #op2_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.op2')

        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'ofprand1.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ofprand1.op2')
        #fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        include_results = None
        #include_results = 'ato.displacements'
        #include_results = 'psd.cquad8_force'
        #include_results = 'ato.cbend_force'
        op2, is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=False, xref_safe=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False,
            is_sort2=False, is_nx=None, delete_f06=True,
            subcases=None,
            exclude_results=None,
            include_results=include_results,
            short_stats=False,
            compare=True, debug=False, binary_debug=True,
            quiet=True,
            stop_on_failure=True, dev=False,
            build_pandas=False, log=log)  # TODO: enable pandas...
        if IS_PANDAS:
            op2.op2_results.stress.cbush_stress[1].build_dataframe()
        #op2 = read_op2_geom(op2_filename, debug=False)
        #op2.write_f06(f06_filename)
        #os.remove(f06_filename)

    def test_bdf_op2_other_02(self):
        """checks ac10707a.bdf, which is an acoustic problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'ac10707a.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ac10707a.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, xref=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_03(self):
        """checks ac10901a.bdf, which is an acoustic problem"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'ac10901a.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ac10901a.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, xref=False, log=log)
        #model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)  # TODO:enable pandas

    def test_bdf_op2_other_04(self):
        """checks v10111.bdf, which is an conical problem"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'v10111.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'v10111.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log)
        #model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_05(self):
        """checks ar29sadl.bdf, which is an CBUSH1D problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'ar29sadl.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ar29sadl.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model, run_renumber=False)  # excite_id

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_06(self):
        """
        checks randvar2.bdf, which is an CTRIAX problem
        CBEND - 2
        """
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'randvar2.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'randvar2.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_07(self):
        """checks randvar2.bdf, which is an CTRIAX problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'v12902.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'v12902.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_08(self):
        """checks randvar2.bdf, which is an CTRIAX problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'mne7a.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'mne7a.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model, run_convert=False)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_09(self):
        """checks sdbush10.bdf, which is a ??? problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'sdbush10.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'sdbush10.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model, run_convert=False, run_renumber=False)  # excite_id

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_10(self):
        """checks v10112.bdf, which is an ??? problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'v10112.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'v10112.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_11(self):
        """checks cbus129.bdf, which is an transient real/oes/cbush problem"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'cbus129.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'cbus129.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log)
        #model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=True,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_12(self):
        """checks api3.bdf, which is a ??? problem"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'api3.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'api3.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, run_skin_solids=False, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model, run_test_bdf=False, run_mass_properties=False)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_13(self):
        """checks v10601s.bdf, which is an superelement complex cgap_force example"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'v10601s.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'v10601s.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False, write_hdf5=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=IS_PANDAS, log=log)

    def test_bdf_op2_other_14(self):
        """checks sbuckl2a.bdf, which is an buckling elgenvalues example"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'sbuckl2a.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'sbuckl2a.op2')
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_15(self):
        """checks dbxdra2.bdf, which is an ComplexTriaxStressArray example"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'dbxdra2.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'dbxdra2.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_16(self):
        """checks phsflux4.bdf, which tests feedge"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'phsflux4.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'phsflux4.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model = BDF(debug=True, log=log, mode='msc')
        model._add_disabled_cards()
        model.read_bdf(bdf_filename=bdf_filename, validate=True, xref=False,
                       punch=False, read_includes=True, save_file_structure=False,
                       encoding=None)
        nid = 1
        sid = 601
        mag = 1.0
        fxyz = [0., 0., 1.]
        model.add_force(sid, nid, mag, fxyz)
        model.log.warning('faking GMLOAD sid=601')
        model.safe_cross_reference()
        save_load_deck(model, run_renumber=False)  # GMSPC

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_17(self):
        """checks cc508a.bdf, which tests feface, gmcurv, gmsurf"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'cc508a.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'cc508a.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model = BDF(debug=True, log=log, mode='msc')
        model._add_disabled_cards()
        model.read_bdf(bdf_filename=bdf_filename, validate=True, xref=False,
                       punch=False, read_includes=True, save_file_structure=False,
                       encoding=None)
        model.safe_cross_reference()
        save_load_deck(model, remove_disabled_cards=False, run_save_load_hdf5=False)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_18(self):
        """checks see101nd.bdf, which tests superelement cards"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'see101nd.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'see101nd.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()
        #save_load_deck(model)

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_19(self):
        """checks see101ta.bdf, which tests superelement cards"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'see101ta.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'see101ta.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()
        save_load_deck(model, run_renumber=False, run_op2_writer=False)

        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_20(self):
        """checks gpst17.bdf, which tests GridPointStressesVolumeDirectArray"""
        log = get_logger(level='error')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'gpst17.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'gpst17.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = BDF(debug=True, log=log, mode='msc')
        model._add_disabled_cards()
        model.read_bdf(bdf_filename=bdf_filename, validate=True, xref=False,
                       punch=False, read_includes=True, save_file_structure=False,
                       encoding=None)
        model.safe_cross_reference()

        # run_op2_reader - super strange PLOAD4 bug
        #save_load_deck(model, xref=False, run_renumber=False,
                       #run_op2_reader=False, run_op2_writer=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_21(self):
        """checks cqra00366.bdf, which tests RealBush1DStressArray"""
        #log = get_logger(level='error')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'cqra00366.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'cqra00366.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_22(self):
        """checks dbxdra7.bdf, which tests RealBush1DStressArray"""
        #log = get_logger(level='error')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'dbxdra7.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'dbxdra7.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_23(self):
        """checks ehbus69.bdf, which tests RealBush1DStressArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'ehbus69.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ehbus69.op2')

        ##  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename,
            #run_skin_solids=False, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        save_load_deck(model, run_test_bdf=False,
                       run_mass_properties=False, run_mirror=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_24(self):
        """checks tst1d3.bdf, which tests RealBar10NodesStrainArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'tst1d3.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'tst1d3.op2')

        ##  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename,
            run_skin_solids=False, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        save_load_deck(model, run_test_bdf=False, run_convert=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_25(self):
        """checks trncomp12.bdf, which tests FailureIndicesArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'trncomp12.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'trncomp12.op2')

        ##  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename,
            run_skin_solids=False, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        save_load_deck(model, run_renumber=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=False, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_26(self):
        """checks tr1091x.bdf, which tests RealBendForceArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'tr1091x.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'tr1091x.op2')

        ##  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename,
            #run_skin_solids=False, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model = BDF(debug=True, log=log, mode='msc')
        model._add_disabled_cards()
        model.read_bdf(bdf_filename=bdf_filename, validate=True, xref=False,
                       punch=False, read_includes=True, save_file_structure=False,
                       encoding=None)
        assert model.card_count['CBEAM'] == 2, model.card_count
        assert model.card_count['CBEND'] == 2, model.card_count
        #model.safe_cross_reference()

        #save_load_deck(model, run_renumber=False)

        log = get_logger(level='warning')
        exclude_results = '*cbend_stress'
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=exclude_results, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_27(self):
        """checks ac10804.bdf, which tests ComplexPlateStressArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'ac10804.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'ac10804.op2')

        ##  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename,
            #run_skin_solids=False, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model = BDF(debug=True, log=log, mode='msc')
        model._add_disabled_cards()
        model.read_bdf(bdf_filename=bdf_filename, validate=True, xref=False,
                       punch=False, read_includes=True, save_file_structure=False,
                       encoding=None)
        #assert model.card_count['CBEND'] == 2, model.card_count
        #model.safe_cross_reference()

        #save_load_deck(model, run_renumber=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_28(self):
        """checks sdr11se_s2dc.bdf, which tests ComplexCBushStressArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'sdr11se_s2dclg.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'sdr11se_s2dclg.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=False, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=False, log=log)

    def test_bdf_op2_other_30(self):
        """checks rot063akd2s_107.bdf, which tests CampbellDiagram"""
        log = get_logger(level='info')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'rot063akd2s_107.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'rot063akd2s_107.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        #model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=True,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_other_31(self):
        """checks htrussx.bdf, which tests getting rid of the Panel"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'htrussx.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'htrussx.op2')

        #  can't parse replication
        #unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            #'', bdf_filename, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_bdf_op2_64_bit(self):
        """
        checks d173.bdf, which tests MSC Nastran 64-bit without the
        op2.is_interlaced flag
        """
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'msc', '64_bit', 'd173.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'msc', '64_bit', 'd173.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_vibroacoustics1(self):
        """
        checks test_vba.bdf, which tests NX Nastran vibroacoustics
        """
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'nx', 'test_vba', 'test_vba.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'nx', 'test_vba', 'test_vba.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        log = get_logger(level='debug')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=True,
                is_mag_phase=False,
                is_sort2=False, is_nx=True, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_op2_superelement_1(self):
        """
        extse04c
        EXTRN
        """
        log = get_logger(level='info')
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'extse04c.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'superelements', 'extse04c.op2')

        #  can't parse replication
        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        unused_model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        #model.safe_cross_reference()

        #save_load_deck(model, run_save_load=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=True,
                write_f06=True, write_op2=True,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    def test_set_results(self):
        """tests setting only a subset of results"""
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.op2')
        f06_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.test_op2.f06')

        op2 = OP2(debug=False, log=log)
        op2.set_results('stress')
        op2.read_op2(op2_filename)
        stress = op2.op2_results.stress
        self.assertEqual(len(stress.cpenta_stress), 0, len(stress.cpenta_stress))
        self.assertEqual(len(stress.chexa_stress), 0, len(stress.chexa_stress))
        self.assertEqual(len(stress.ctetra_stress), 1, len(stress.ctetra_stress))
        self.assertEqual(len(op2.displacements), 0, len(op2.displacements))

        op2 = OP2(debug=False, log=log)
        op2.set_results(['stress', 'displacements'])
        op2.read_op2(op2_filename)
        stress = op2.op2_results.stress
        self.assertEqual(len(stress.cpenta_stress), 0, len(stress.cpenta_stress))
        self.assertEqual(len(stress.chexa_stress), 0, len(stress.chexa_stress))
        self.assertEqual(len(stress.ctetra_stress), 1, len(stress.ctetra_stress))
        self.assertEqual(len(op2.displacements), 1, len(op2.displacements))
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_solid_bending_01(self):
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        f06_filename = os.path.join(folder, 'solid_bending.test_op2.f06')
        debug = False
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = f'{model}.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)

        read_op2(op2_filename, debug=False, log=log)
        run_op2(op2_filename, write_bdf=False,
                write_f06=True,
                debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
                load_as_h5=False, log=log)
        assert os.path.exists(debug_file), os.listdir(folder)

        op2 = run_op2(op2_filename, make_geom=False, write_bdf=False,
                      write_f06=True,
                      debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
                      build_pandas=True, log=log)[0]
        assert os.path.exists(debug_file), os.listdir(folder)

        subcase = op2.case_control_deck.subcases[1]
        value, options = subcase['TITLE']
        assert options == [], subcase
        value, options = subcase['SUBTITLE']
        assert options == [], subcase
        value, options = subcase['LABEL']
        assert options == [], subcase

        op2.save('op2_model.pik')
        op2_load = OP2()
        op2_load.load('op2_model.pik')
        os.remove('op2_model.pik')
        os.remove(debug_file)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    @unittest.skipIf(not IS_H5PY, "No h5py")
    def test_op2_solid_bending_02(self):
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2 = OP2(debug=False, log=log, debug_file=None, mode=None)
        op2.load_as_h5 = True
        op2.read_op2(op2_filename=op2_filename, combine=True,
                     build_dataframe=False, skip_undefined_matrices=False,
                     encoding=None)
        #op2 = read_op2(op2_filename, debug=False)
        del op2

    def test_op2_solid_bending_02_geom(self):
        log = get_logger(level='warning')
        #log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        hdf5_filename = os.path.join(folder, 'solid_bending.test_op2_solid_bending_02_geom.h5')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False, write_hdf5=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, build_pandas=True, log=log,
            name='_02_geom')
        if IS_PANDAS:
            assert op2.displacements[1].data_frame is not None

        op2.print_subcase_key()

        if IS_H5PY:  # TODO: broken on old packages test only...
            op2.export_hdf5_filename(hdf5_filename)  # fails...
            op2b = OP2(debug=False, log=log)
            op2b.load_hdf5_filename(hdf5_filename, combine=True)
            op2b.print_subcase_key()

    def test_op2_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'static_solid_shell_bar.test_op2.f06')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, build_pandas=True, log=log,
            name='_01_geom')
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_mode_solid_shell_bar_01_geom(self):
        """tests reading op2 geometry"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        subcases = [1]
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, log=log)
        op2.get_op2_stats(short=False)
        op2.get_op2_stats(short=True)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_buckling_solid_shell_bar_01_geom(self):
        """single subcase buckling"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling_solid_shell_bar.op2')
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, log=log)

        f06_filename = os.path.join(folder, 'buckling_solid_shell_bar.test_op2_sort2.f06')
        op2.write_f06(f06_filename, is_mag_phase=False, is_sort1=False,
                      #delete_objects=True,
                      end_flag=False, quiet=True,
                      repr_check=False, close=True)

        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_buckling_solid_shell_bar_02_geom(self):
        """multi subcase buckling"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'buckling2_solid_shell_bar.op2')
        unused_op2 = read_op2_geom(op2_filename, debug=False, log=log)
        subcases = 1
        op2 = read_op2_geom(op2_filename, debug=False, subcases=subcases, log=log)
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
            subcases=subcases, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=True, log=log)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = 2
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, log=log)
        assert len(op2.displacements) == 0, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

        subcases = [1, 2]
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=False, write_bdf=False,
            write_f06=True, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=subcases, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, log=log)
        assert len(op2.displacements) == 1, len(op2.displacements)
        assert len(op2.eigenvectors) == 1, len(op2.eigenvectors)

    def test_op2_transient_solid_shell_bar_01_geom(self):
        """transient test"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'transient_solid_shell_bar.test_op2.f06')
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False,
            write_f06=False, write_op2=False,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, build_pandas=True, log=log)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_frequency_solid_shell_bar_01_geom(self):
        """frequency test"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')
        unused_op2 = read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True,
            write_f06=False, write_op2=True,
            is_mag_phase=False, is_sort2=False, delete_f06=False,
            subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, binary_debug=False,
            quiet=True, stop_on_failure=True,
            dev=False, log=log)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

    def test_op2_transfer_function_01(self):
        """tests the transfer function cards work"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'transfer_function')
        #bdf_filename = os.path.join(folder, 'actuator_tf_modeling.bdf')
        op2_filename = os.path.join(folder, 'actuator_tf_modeling.op2')
        f06_filename = os.path.join(folder, 'freq_solid_shell_bar.test_op2.f06')

        unused_op2 = read_op2_geom(op2_filename, debug=False, log=log)

        debug = False
        write_bdf = True
        write_f06 = True
        make_geom = True
        op2, unused_is_passed = run_op2(
            op2_filename, write_bdf=write_bdf, make_geom=make_geom,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)

        #fem1, fem2, diff_cards = self.run_bdf(folder, bdf_filename, log=log)
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
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'aero', 'monpnt3')
        op2_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.op2')
        f06_filename = os.path.join(folder, 'Monitor_Points_data_LINE5000000_10FREQs.test_op2.f06')
        op2 = OP2Geom(log=log, debug=False, debug_file='temp.debug')
        op2.read_op2(op2_filename)
        monitor3 = op2.op2_results.monitor3
        assert len(monitor3.frequencies) == 11, monitor3
        str(monitor3)
        op2.write_f06(f06_filename)
        os.remove(f06_filename)
        os.remove('temp.debug')

    def test_op2_solid_shell_bar_01(self):
        """tests sol_101_elements/static_solid_shell_bar.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
            load_as_h5=False, log=log)

        force = op2.op2_results.force
        stress = op2.op2_results.stress
        #strain = op2.op2_results.strain

        isubcase = 1
        rod_force = force.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (1, 2, 2), rod_force.data.shape

        rod_stress = stress.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (1, 2, 4), rod_stress.data.shape

        cbar_force = force.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (1, 1, 8), cbar_force.data.shape

        cbar_stress = stress.cbar_stress[isubcase]
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (1, 1, 15), cbar_stress.data.shape

        cbeam_force = force.cbeam_force[isubcase]
        assert cbeam_force.nelements == 1, cbeam_force.nelements
        assert cbeam_force.data.shape == (1, 2, 8), cbeam_force.data.shape

        cbeam_stress = stress.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # wrong
        assert cbeam_stress.data.shape == (1, 2, 8), cbeam_stress.data.shape

        cquad4_force = force.cquad4_force[isubcase]
        assert cquad4_force.nelements == 4, cquad4_force.nelements
        assert cquad4_force.data.shape == (1, 20, 8), cquad4_force.data.shape

        cquad4_stress = stress.cquad4_stress[isubcase]
        # TODO: should this be 4; yes by actual count...
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (1, 20, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_force = force.ctria3_force[isubcase]
        assert ctria3_force.nelements == 8, ctria3_force.nelements
        assert ctria3_force.data.shape == (1, 8, 8), ctria3_force.data.shape

        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        ctetra_stress = stress.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (1, 10, 10), ctetra_stress.data.shape
        assert ctetra_stress.is_von_mises, ctetra_stress

        cpenta_stress = stress.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (1, 14, 10), cpenta_stress.data.shape
        assert cpenta_stress.is_von_mises, cpenta_stress

        chexa_stress = stress.chexa_stress[isubcase]
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

    def test_op2_solid_shell_bar_01_straincurvature(self):
        """tests sol_101_elements/static_solid_shell_bar_straincurve.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        isubcase = 1
        stress = op2.op2_results.stress
        strain = op2.op2_results.strain
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = strain.ctria3_strain[isubcase]
        sword = get_scode_word(ctria3_strain.s_code, ctria3_strain.stress_bits)
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance, '%s\n%s' % (ctria3_strain, sword)
        assert ctria3_strain.is_von_mises, '%s\n%s' % (ctria3_strain, sword)

        cquad4_strain = strain.cquad4_strain[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance, cquad4_strain
        assert cquad4_strain.is_von_mises, cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance(self):
        """tests sol_101_elements/static_solid_shell_bar_fiberdist.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        isubcase = 1
        stress = op2.op2_results.stress
        strain = op2.op2_results.strain
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = strain.ctria3_strain[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance, ctria3_strain
        assert ctria3_strain.is_von_mises, ctria3_strain

        cquad4_strain = strain.cquad4_strain[isubcase]
        sword = get_scode_word(cquad4_strain.s_code, cquad4_strain.stress_bits)
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), '%s\n%s' % (cquad4_strain.data.shape, sword)
        assert cquad4_strain.is_fiber_distance, '%s\n%s' % (cquad4_strain, sword)
        assert cquad4_strain.is_von_mises, '%s\n%s' % (cquad4_strain, sword)

    def test_op2_solid_shell_bar_01_straincurvature_shear(self):
        """tests sol_101_elements/static_solid_shell_bar_straincurve_shear.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        isubcase = 1
        stress = op2.op2_results.stress
        strain = op2.op2_results.strain
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert not ctria3_stress.is_von_mises, ctria3_stress

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2?
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert not cquad4_stress.is_von_mises, cquad4_stress

        ctria3_strain = strain.ctria3_strain[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert not ctria3_strain.is_fiber_distance, ctria3_strain
        assert not ctria3_strain.is_von_mises, ctria3_strain

        cquad4_strain = strain.cquad4_strain[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2?
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert not cquad4_strain.is_fiber_distance, cquad4_strain
        assert not cquad4_strain.is_von_mises, cquad4_strain

    def test_op2_solid_shell_bar_01_fiberdistance_shear(self):
        """tests sol_101_elements/static_solid_shell_bar_fiberdist_shear.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        isubcase = 1
        stress = op2.op2_results.stress
        strain = op2.op2_results.strain
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (1, 8, 8), ctria3_stress.data.shape
        assert ctria3_stress.is_fiber_distance, ctria3_stress
        assert ctria3_stress.is_von_mises is False, ctria3_stress

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 4, cquad4_stress.nelements # TODO: this should be 2
        assert cquad4_stress.data.shape == (1, 4, 8), cquad4_stress.data.shape
        assert cquad4_stress.is_fiber_distance, cquad4_stress
        assert cquad4_stress.is_von_mises is False, cquad4_stress

        ctria3_strain = strain.ctria3_strain[isubcase]
        assert ctria3_strain.nelements == 8, ctria3_strain.nelements
        assert ctria3_strain.data.shape == (1, 8, 8), ctria3_strain.data.shape
        assert ctria3_strain.is_fiber_distance, ctria3_strain
        assert ctria3_strain.is_von_mises is False, ctria3_strain

        cquad4_strain = strain.cquad4_strain[isubcase]
        assert cquad4_strain.nelements == 4, cquad4_strain.nelements # TODO: this should be 2
        assert cquad4_strain.data.shape == (1, 4, 8), cquad4_strain.data.shape
        assert cquad4_strain.is_fiber_distance, cquad4_strain
        assert cquad4_strain.is_von_mises is False, cquad4_strain

    def test_op2_solid_shell_bar_mode(self):
        """tests sol_101_elements/mode_solid_shell_bar.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        stress = op2.op2_results.stress
        force = op2.op2_results.force
        isubcase = 1
        rod_force = force.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (3, 2, 2), rod_force.data.shape

        rod_stress = stress.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (3, 2, 4), rod_stress.data.shape

        cbar_force = force.cbar_force[isubcase]
        str(cbar_force.data_frame)
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (3, 1, 8), cbar_force.data.shape

        cbar_stress = stress.cbar_stress[isubcase]
        assert cbar_stress.nelements == 3, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (3, 1, 15), cbar_stress.data.shape

        cbeam_stress = stress.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (3, 2, 8), cbeam_stress.data.shape

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (3, 20, 8), cquad4_stress.data.shape

        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (3, 8, 8), ctria3_stress.data.shape

        ctetra_stress = stress.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (3, 10, 10), ctetra_stress.data.shape

        cpenta_stress = stress.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (3, 14, 10), cpenta_stress.data.shape

        chexa_stress = stress.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (3, 9, 10), chexa_stress.data.shape

        if IS_PANDAS:
            cbar_force.build_dataframe()
        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_buckling(self):
        """tests sol_101_elements/buckling_solid_shell_bar.op2"""
        log = get_logger(level='warning')
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
        read_op2_geom(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)

        isubcases = [(1, 1, 1, 0, 0, '', ''), (1, 8, 1, 0, 0, '', '')]
        isubcase = isubcases[1]

        force = op2.op2_results.force
        stress = op2.op2_results.stress

        try:
            rod_force = force.crod_force[isubcase]
        except KeyError:
            msg = 'isubcase=%s was not found\nkeys=%s' % (isubcase, force.crod_force.keys())
            raise KeyError(msg)
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (4, 2, 2), rod_force.data.shape

        rod_stress = stress.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (4, 2, 4), rod_stress.data.shape

        cbar_force = force.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (4, 1, 8), cbar_force.data.shape

        cbar_stress = stress.cbar_stress[isubcase]
        assert cbar_stress.nelements == 4, cbar_stress.nelements  # TODO: wrong
        assert cbar_stress.data.shape == (4, 1, 15), cbar_stress.data.shape

        cbeam_stress = stress.cbeam_stress[isubcase]
        assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        assert cbeam_stress.data.shape == (4, 2, 8), cbeam_stress.data.shape

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 20, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (4, 20, 8), cquad4_stress.data.shape

        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements
        assert ctria3_stress.data.shape == (4, 8, 8), ctria3_stress.data.shape

        ctetra_stress = stress.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (4, 10, 10), ctetra_stress.data.shape

        cpenta_stress = stress.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (4, 14, 10), cpenta_stress.data.shape

        chexa_stress = stress.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (4, 9, 10), chexa_stress.data.shape

        if IS_PANDAS:
            cbar_force.build_dataframe()
            str(cbar_force.data_frame)

        assert os.path.exists(debug_file), os.listdir(folder)
        os.remove(debug_file)

    def test_op2_solid_shell_bar_freq(self):
        """tests sol_101_elements/freq_solid_shell_bar.op2"""
        log = get_logger(level='warning')
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
        read_op2(op2_filename, debug=False, log=log)
        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True, log=log)
        isubcase = 1

        force = op2.op2_results.force
        stress = op2.op2_results.stress

        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape

        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        rod_force = force.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape

        #if op2.is_nx:
            #return
        rod_stress = stress.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (7, 2, 2), rod_stress.data.shape

        cbar_force = force.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (7, 1, 8), cbar_force.data.shape

        cbar_stress = stress.cbar_stress[isubcase]
        assert cbar_stress.nelements == 1, cbar_stress.nelements
        assert cbar_stress.data.shape == (7, 1, 9), cbar_stress.data.shape

        #print(stress.cbeam_stress.keys())
        #cbeam_stress = op2.cbeam_stress[isubcase]
        #assert cbeam_stress.nelements == 2, cbeam_stress.nelements
        #assert cbeam_stress.data.shape == (7, 2, 8), cbeam_stress.data.shape

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.data.shape == (7, 40, 3), cquad4_stress.data.shape
        assert cquad4_stress.nelements == 40, cquad4_stress.nelements # TODO: wrong?

        #print(stress.ctria3_stress.keys())
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.data.shape == (7, 16, 3), ctria3_stress.data.shape
        assert ctria3_stress.nelements == 8, ctria3_stress.nelements # TODO: wrong

        ctetra_stress = stress.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (7, 10, 6), ctetra_stress.data.shape

        cpenta_stress = stress.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (7, 14, 6), cpenta_stress.data.shape

        chexa_stress = stress.chexa_stress[isubcase]
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

    def test_op2_solid_shell_bar_transient(self):
        """
        MSC 2005r2 Tables : GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPTS, DYNAMICS, DIT
                            OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        NX 10 Tables : PVT0, CASECC, GEOM1S, GEOM2S, GEOM3S, GEOM4S, EPTS, MPTS,
                       DYNAMICS, BGPDTS, EQEXINS, DIT,
                       OQG1, OUGV1, OGPFB1, OEF1X, OES1X1, OSTR1X, OPG1
        """
        log = get_logger(level='warning')
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

        op2 = OP2()
        #result_set = op2._result_set
        assert op2._results.is_not_saved('stress.ctetra_stress') == False

        read_op2(op2_filename, debug=debug, log=log)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom, write_bdf=write_bdf,
            write_f06=write_f06,
            debug=debug, stop_on_failure=True, binary_debug=True, quiet=True,
            build_pandas=True, log=log)
        isubcase = 1

        # rod_force = op2.crod_force[isubcase]
        # assert rod_force.nelements == 2, rod_force.nelements
        # assert rod_force.data.shape == (7, 2, 2), rod_force.data.shape


        # isubcases = [(1, 1, 1, 0, 'DEFAULT'), (1, 8, 1, 0, 'DEFAULT')]
        # isubcase = isubcases[1]

        force = op2.op2_results.force
        stress = op2.op2_results.stress
        #strain = op2.op2_results.strain

        rod_force = force.crod_force[isubcase]
        assert rod_force.nelements == 2, rod_force.nelements
        assert rod_force.data.shape == (21, 2, 2), rod_force.data.shape

        rod_stress = stress.crod_stress[isubcase]
        assert rod_stress.nelements == 2, rod_stress.nelements
        assert rod_stress.data.shape == (21, 2, 4), rod_stress.data.shape

        cbar_force = force.cbar_force[isubcase]
        assert cbar_force.nelements == 1, cbar_force.nelements
        assert cbar_force.data.shape == (21, 1, 8), cbar_force.data.shape

        cbar_stress = stress.cbar_stress[isubcase]
        assert cbar_stress.nelements == 21, cbar_stress.nelements # 1-wrong
        assert cbar_stress.data.shape == (21, 1, 15), cbar_stress.data.shape

        #print(op2.cbeam_stress.keys())
        # cbeam_stress = op2.cbeam_stress[isubcase]
        # assert cbeam_stress.nelements == 11, cbeam_stress.nelements  # TODO: wrong
        # assert cbeam_stress.data.shape == (7, 11, 8), cbeam_stress.data.shape

        cquad4_stress = stress.cquad4_stress[isubcase]
        assert cquad4_stress.nelements == 40, cquad4_stress.nelements
        assert cquad4_stress.data.shape == (21, 40, 8), cquad4_stress.data.shape

        #print(op2.ctria3_stress.keys())
        ctria3_stress = stress.ctria3_stress[isubcase]
        assert ctria3_stress.nelements == 16, ctria3_stress.nelements  # TODO: 8-wrong
        assert ctria3_stress.data.shape == (21, 16, 8), ctria3_stress.data.shape

        assert len(stress.ctetra_stress) == 1, stress.ctetra_stress
        ctetra_stress = stress.ctetra_stress[isubcase]
        assert ctetra_stress.nelements == 2, ctetra_stress.nelements
        assert ctetra_stress.data.shape == (21, 10, 10), ctetra_stress.data.shape

        cpenta_stress = stress.cpenta_stress[isubcase]
        assert cpenta_stress.nelements == 2, cpenta_stress.nelements
        assert cpenta_stress.data.shape == (21, 14, 10), cpenta_stress.data.shape

        chexa_stress = stress.chexa_stress[isubcase]
        assert chexa_stress.nelements == 1, chexa_stress.nelements
        assert chexa_stress.data.shape == (21, 9, 10), chexa_stress.data.shape

        grid_point_forces = op2.grid_point_forces[isubcase]
        #print(grid_point_forces._ntotals)
        assert grid_point_forces.ntotal == 130, grid_point_forces.ntotal
        assert grid_point_forces.data.shape == (21, 130, 6), grid_point_forces.data.shape

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
        test_op2(argv, show_args=False)

    def test_op2_good_sine_01(self):
        """tests freq_sine/good_sine.op2"""
        op2_filename = os.path.join(MODEL_PATH, 'freq_sine', 'good_sine.op2')
        make_geom = False
        write_bdf = False
        write_f06 = False
        log = get_logger(level='warning')
        read_op2(op2_filename, log=log)
        build_pandas = False # IS_TRANSIENT_PANDAS
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

        #if IS_PANDAS:
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
        bdf._add_disabled_cards()
        bdf.read_bdf(bdf_filename, xref=False)

        debug = False
        debug_file = 'debug.out'
        read_op2(op2_filename, debug=False, log=log)
        op2 = OP2(debug=debug, debug_file=debug_file, log=log)
        op2.read_op2(op2_filename)
        assert os.path.exists(debug_file), os.listdir('.')
        os.remove(debug_file)

        _verify_ids(bdf, op2, isubcase=1)

        msg = ''
        #op2.subcase_key = {1: [1], 2: [2], 3: [3], 4: [4], 5: [5], 6: [6], 7: [7], 8: [8], 9: [9], 10: [10], 11: [11], 12: [12], 13: [13], 14: [14], 15: [15], 16: [16], 17: [17], 18: [18], 19: [19], 20: [20], 21: [21], 22: [22], 23: [23], 24: [24], 25: [25], 26: [26], 27: [27], 28: [28], 29: [29], 30: [30], 31: [31], 32: [32], 33: [33], 34: [34], 35: [35], 36: [36], 37: [37], 38: [38], 39: [39], 40: [40], 41: [41], 42: [42], 43: [43], 44: [44], 45: [45], 46: [46], 47: [47], 48: [48], 49: [49], 50: [50], 51: [51], 52: [52], 53: [53], 54: [54], 55: [55], 56: [56], 57: [57], 58: [58], 59: [59], 60: [60], 61: [61], 62: [62], 63: [63], 64: [64], 65: [65], 66: [66], 67: [67], 68: [68], 69: [69], 70: [70], 71: [71], 72: [72], 73: [73], 74: [74], 75: [75], 76: [76], 77: [77], 78: [78], 79: [79], 80: [80], 81: [81], 82: [82], 83: [83], 84: [84], 85: [85], 86: [86], 87: [87], 88: [88], 89: [89], 90: [90], 91: [91], 92: [92], 93: [93], 94: [94], 95: [95], 96: [96], 97: [97], 98: [98], 99: [99], 100: [100], 101: [101], 102: [102], 103: [103], 104: [104], 105: [105], 106: [106], 107: [107], 108: [108], 109: [109], 110: [110], 111: [111], 112: [112], 113: [113], 114: [114], 115: [115], 116: [116], 117: [117], 118: [118], 119: [119], 120: [120]}
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

    def test_cbar_orientation(self):
        """tests cbar_orientation.op2"""
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'cbar_orientation.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'bars', 'cbar_orientation.op2')
        make_geom = True
        write_bdf = True
        write_f06 = True
        log = get_logger(level='warning')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename,
            run_skin_solids=False, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log)
        save_load_deck(model)

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, log=log)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom,
            write_bdf=write_bdf, read_bdf=True,
            write_f06=write_f06, write_op2=True, write_hdf5=True,
            is_mag_phase=False, is_sort2=False, is_nx=None, delete_f06=False,
            build_pandas=False, subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, log=log, binary_debug=True,
            quiet=True, stop_on_failure=True, dev=False, xref_safe=False,
            post=None, load_as_h5=True)
        #compare_elements(fem1, op2)

        assert os.path.exists(debug_file), os.listdir(os.path.dirname(op2_filename))
        os.remove(debug_file)

    def test_op2_cbush_01(self):
        """tests cbush/cbush.op2"""
        op2_filename = os.path.join(MODEL_PATH, 'unit', 'cbush', 'cbush.op2')
        make_geom = True
        write_bdf = True
        write_f06 = True
        log = get_logger(level='warning')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        debug_file = model + '.debug.out'

        if os.path.exists(debug_file):
            os.remove(debug_file)
        read_op2(op2_filename, log=log)

        op2, unused_is_passed = run_op2(
            op2_filename, make_geom=make_geom,
            write_bdf=write_bdf, read_bdf=None,
            write_f06=write_f06, write_op2=True, write_hdf5=True,
            is_mag_phase=False, is_sort2=False, is_nx=None, delete_f06=True,
            build_pandas=True, subcases=None, exclude_results=None, short_stats=False,
            compare=True, debug=False, log=log, binary_debug=True,
            quiet=True, stop_on_failure=True, dev=False, xref_safe=False,
            post=None, load_as_h5=True)
        isubcase = 1

        cbush_stress = op2.op2_results.stress.cbush_stress[isubcase]
        assert cbush_stress.nelements == 1, cbush_stress.nelements
        assert cbush_stress.data.shape == (1, 1, 6), cbush_stress.data.shape

        cbush_strain = op2.op2_results.strain.cbush_strain[isubcase]
        assert cbush_strain.nelements == 1, cbush_strain.nelements
        assert cbush_strain.data.shape == (1, 1, 6), cbush_strain.data.shape

        cbush_force = op2.cbush_force[isubcase]
        assert cbush_force.nelements == 1, cbush_force.nelements
        assert cbush_force.data.shape == (1, 1, 6), cbush_force.data.shape
        if IS_PANDAS:
            op2.build_dataframe()

        assert os.path.exists(debug_file), os.listdir(os.path.dirname(op2_filename))
        os.remove(debug_file)

    def test_cgap_01(self):
        """checks cc188b.bdf"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'cc188b.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'cc188b.op2')
        unused_fem1, unused_fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=False,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=True, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

        model = read_bdf(bdf_filename, debug=False, log=log)
        save_load_deck(model, run_op2_reader=False)

    def test_bdf_op2_cbush_nonlinear(self):
        """checks tr1091x.bdf, which tests RealBendForceArray"""
        log = get_logger(level='info')
        bdf_filename = os.path.join(OP2_TEST_PATH, 'cbush_nonlinear', 'cbush_106_large_disp.bdf')
        op2_filename = os.path.join(OP2_TEST_PATH, 'cbush_nonlinear', 'cbush_106_large_disp.op2')

        unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            '', bdf_filename, log=log)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        model.safe_cross_reference()

        save_load_deck(model, run_test_bdf=False)

        log = get_logger(level='warning')
        run_op2(op2_filename, make_geom=True, write_bdf=True, read_bdf=True,
                write_f06=True, write_op2=True,
                is_mag_phase=False,
                is_sort2=False, is_nx=None, delete_f06=True,
                subcases=None, exclude_results=None, short_stats=False,
                compare=False, debug=False, binary_debug=True,
                quiet=True,
                stop_on_failure=True, dev=False,
                build_pandas=True, log=log)

    @unittest.expectedFailure
    def test_set_times_01(self):
        """specify the modes to extract"""
        log = get_logger(level='warning')
        model = OP2(debug=False, log=log, debug_file=None, mode='msc')

        # subcase : times
        times = {1 : [2, 3]}
        model.set_transient_times(times)
        op2_filename = os.path.abspath(os.path.join(MODEL_PATH, 'sol_101_elements',
                                                    'mode_solid_shell_bar.op2'))
        model.read_op2(op2_filename, combine=True, build_dataframe=False,
                       skip_undefined_matrices=False,
                       encoding=None)
        isubcase = 1
        #print(model.get_op2_stats(short=False))
        eigenvector = model.eigenvectors[isubcase]
        #print(eigenvector)
        assert len(eigenvector.modes) == 2, eigenvector.modes

    def test_random_ctria3(self):
        """runs a random test"""
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'random')
        op2_filename = os.path.join(folder, 'random_test_bar_plus_tri.op2')
        f06_filename = os.path.join(folder, 'random_test_bar_plus_tri.test_op2.f06')
        op2 = read_op2_geom(op2_filename, debug=False, log=log)

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
        log = get_logger(level='warning')
        folder = os.path.join(MODEL_PATH, 'random')
        op2_filename = os.path.join(folder, 'rms_tri_oesrmx1.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, debug=False, log=log)
        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=None, write_f06=True,
            write_op2=False, write_hdf5=False, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=False, build_pandas=False, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

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
            write_op2=True, write_hdf5=False, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=False, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
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
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_ogstr(self):
        """test grid_point_strains"""
        log = get_logger(level='warning')
        op2_filename = os.path.join(MODEL_PATH, 'ogs', 'nx_ogstr.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=False, write_hdf5=True, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=False, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_cbeam3_cbend(self):
        """test CBEAM3/CBEND"""
        log = get_logger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'b3bend.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'b3bend.op2')
        model = read_bdf(bdf_filename, debug=False, log=log)
        save_load_deck(model)

        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=True, read_bdf=None, write_f06=True,
            write_op2=False, write_hdf5=True, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_ougv1pat(self):
        """test OUGV1PAT table"""
        log = get_logger(level='warning')
        op2_filename1 = os.path.join(OP2_TEST_PATH, 'ougv1pat', 'winkel_2005r3_ougcord_basic.op2')
        op2_filename2 = os.path.join(OP2_TEST_PATH, 'ougv1pat', 'winkel_2013.1_ougcord_basic.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        unused_op2, unused_is_passed = run_op2(
            op2_filename1, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

        unused_op2, unused_is_passed = run_op2(
            op2_filename2, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    def test_xsop2dir(self):
        """test NX 2019 XSOP2DIR"""
        log = get_logger(level='warning')
        #bdf_filename = os.path.join(MODEL_PATH, 'other', 'extse04c_cnv2_0.bdf')
        op2_filename = os.path.join(MODEL_PATH, 'other', 'extse04c_cnv2_0.op2')
        #model = read_bdf(bdf_filename, encoding='ascii', debug=False, log=log)
        #bdf_filename_out = os.path.join(MODEL_PATH, 'other', 'extse04c_cnv2_0.bdf')
        #model.write_bdf(bdf_filename_out)
        #os.remove(bdf_filename_out)

        unused_op2, unused_is_passed = run_op2(
            op2_filename, make_geom=True, write_bdf=False, read_bdf=False, write_f06=True,
            write_op2=True, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=True)

    def test_sol_106(self):
        """tests SOL 106 pandas bug"""
        log = get_logger(level='warning')
        op2_filename1 = os.path.join(MODEL_PATH, 'bugs', 'sol_106_pandas', 'test.op2')
        #bdf_filename = os.path.join(folder, 'rms_tri_oesrmx1.bdf')
        #unused_op2 = read_op2_geom(op2_filename, xref=False, log=log)

        WRITE_OP2 = False
        unused_op2, unused_is_passed = run_op2(
            op2_filename1, make_geom=True, write_bdf=False, read_bdf=None, write_f06=True,
            write_op2=WRITE_OP2, write_hdf5=IS_H5PY, is_mag_phase=False, is_sort2=False,
            is_nx=None, delete_f06=True, build_pandas=True, subcases=None,
            exclude_results=None, short_stats=False, compare=True, debug=False, log=log,
            binary_debug=True, quiet=True, stop_on_failure=True,
            dev=False, xref_safe=False, post=None, load_as_h5=False)

    #def test_bdf_op2_random_small_plate(self):
        #"""checks small_plate.op2, which tests PSD tables"""
        #log = get_logger(level='info')
        ##bdf_filename = os.path.join(MODEL_PATH, 'small_plate', 'small_plate.dat')
        #op2_filename = os.path.join(MODEL_PATH, 'small_plate', 'small_plate.op2')

        ##  can't parse replication
        ##unused_fem1, unused_fem2, diff_cards = self.run_bdf(
            ##'', bdf_filename, log=log)
        ##diff_cards2 = list(set(diff_cards))
        ##diff_cards2.sort()
        ##assert len(diff_cards2) == 0, diff_cards2

        ##model = read_bdf(bdf_filename, debug=False, log=log, xref=False)
        ##model.safe_cross_reference()

        #save_load_deck(model, run_save_load=True, run_renumber=False)

        #log = get_logger(level='warning')
        #run_op2(op2_filename, make_geom=True, write_bdf=False, read_bdf=False,
                #write_f06=True, write_op2=False,
                #is_mag_phase=False,
                #is_sort2=False, is_nx=None, delete_f06=True,
                #subcases=None, exclude_results=None, short_stats=False,
                #compare=False, debug=False, binary_debug=True,
                #quiet=True,
                #stop_on_failure=True, dev=False,
                #build_pandas=True, log=log)


def _verify_ids(bdf, op2: OP2, isubcase=1):
    """helper function for tests"""
    types = ['CQUAD4', 'CTRIA3', 'CHEXA', 'CPENTA', 'CTETRA', 'CROD', 'CONROD', 'CTUBE']
    out = bdf.get_card_ids_by_card_types(types)

    stress = op2.op2_results.stress
    strain = op2.op2_results.strain
    force = op2.op2_results.force

    card_type = 'CQUAD4'
    if stress.cquad4_stress:
        try:
            case = stress.cquad4_stress[isubcase]
        except KeyError:
            raise KeyError('getting cquad4_stress; isubcase=%s; keys=%s' % (
                isubcase, stress.cquad4_stress.keys()))
        eids = np.unique(case.element_node[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if strain.cquad4_strain:
        try:
            case = strain.cquad4_strain[isubcase]
        except KeyError:
            raise KeyError('getting cquad4_strain; isubcase=%s; keys=%s' % (
                isubcase, strain.cquad4_strain.keys()))
        eids = np.unique(case.element_node[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if strain.cquad4_composite_strain:
        try:
            case = strain.cquad4_composite_strain[isubcase]
        except KeyError:
            raise KeyError('getting cquad4_composite_strain; isubcase=%s; keys=%s' % (
                isubcase, strain.cquad4_composite_strain.keys()))
        eids = np.unique(case.element_layer[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if stress.cquad4_composite_stress:
        try:
            case = stress.cquad4_composite_stress[isubcase]
        except KeyError:
            raise KeyError('getting cquad4_composite_stress; isubcase=%s; keys=%s' % (
                isubcase, stress.cquad4_composite_stress.keys()))

        eids = np.unique(case.element_layer[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if force.cquad4_force:
        try:
            case = force.cquad4_force[isubcase]
        except KeyError:
            raise KeyError('getting cquad4_force; isubcase=%s; keys=%s' % (
                isubcase, op2.cquad4_force.keys()))

        eids = np.unique(case.element)
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)


    card_type = 'CTRIA3'
    if stress.ctria3_stress:
        case = stress.ctria3_stress[isubcase]
        eids = np.unique(case.element_node[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if strain.ctria3_strain:
        case = strain.ctria3_strain[isubcase]
        eids = np.unique(case.element_node[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if strain.ctria3_composite_strain:
        case = strain.ctria3_composite_strain[isubcase]
        eids = np.unique(case.element_layer[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if stress.ctria3_composite_stress:
        case = stress.ctria3_composite_stress[isubcase]
        eids = np.unique(case.element_layer[:, 0])
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)
    if force.ctria3_force:
        case = force.ctria3_force[isubcase]
        eids = np.unique(case.element)
        for eid in eids:
            assert eid in out[card_type], 'eid=%s eids=%s card_type=%s'  % (eid, out[card_type], card_type)


if __name__ == '__main__':  # pragma: no cover
    ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
    if not ON_RTD:
        unittest.main()
