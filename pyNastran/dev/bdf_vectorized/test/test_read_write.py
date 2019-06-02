"""tests the pyNastran solver"""
import os
import unittest
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.test.test_bdf import run_and_compare_fems
from pyNastran.dev.bdf_vectorized.bdf import read_bdf as read_bdfv
from pyNastran.bdf.bdf import read_bdf
#from pyNastran.dev.bdf_vectorized.solver.solver import Solver

PKG_PATH = pyNastran.__path__[0]
test_path = os.path.join(PKG_PATH, '..', 'models')

def read_write_compare(bdf_filename, bdf_filename_out):
    """runs checks on the two bdfs"""
    vmodel = read_bdfv(bdf_filename)
    vmodel.write_bdf(bdf_filename_out)
    run_and_compare_fems(
        bdf_filename, bdf_filename_out, debug=False, xref=True, check=True,
        punch=False, mesh_form=None,
        print_stats=False, encoding=None,
        sum_load=True, size=8, is_double=False,
        stop=False, nastran='', post=-1, dynamic_vars=None,
        quiet=False, dumplines=False, dictsort=False,
        nerrors=0, dev=False, crash_cards=None,
        run_extract_bodies=False,
    )

class TestReadWriteVectorized(unittest.TestCase):
    """tests the vectorized bdf read/write against the non-vectorized version"""

    @staticmethod
    def test_solid_bending():
        """vectorized vs. standard test on solid_bending.bdf"""
        bdf_filename = os.path.join(test_path, 'solid_bending', 'solid_bending.bdf')
        bdf_filename_out = os.path.join(test_path, 'solid_bending', 'solid_bending2.bdf')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_static_solid_shell_bar():
        """vectorized vs. standard test on static_solid_shell_bar.bdf"""
        bdf_filename = os.path.join(test_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        bdf_filename_out = os.path.join(test_path, 'sol_101_elements', 'static_solid_shell_bar2.bdf')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_mode_solid_shell_bar():
        """vectorized vs. standard test on mode_solid_shell_bar.bdf"""
        bdf_filename = os.path.join(test_path, 'sol_101_elements', 'mode_solid_shell_bar.bdf')
        bdf_filename_out = os.path.join(test_path, 'sol_101_elements', 'mode_solid_shell_bar2.bdf')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_freq_solid_shell_bar():
        """vectorized vs. standard test on freq_solid_shell_bar.bdf"""
        bdf_filename = os.path.join(test_path, 'sol_101_elements', 'freq_solid_shell_bar.bdf')
        bdf_filename_out = os.path.join(test_path, 'sol_101_elements', 'freq_solid_shell_bar2.bdf')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_bwb():
        """vectorized vs. standard test on bwb_saero.bdf"""
        bdf_filename = os.path.join(test_path, 'bwb', 'bwb_saero.bdf')
        bdf_filename_out = os.path.join(test_path, 'bwb', 'bwb_saero2.bdf')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_isat_01():
        """vectorized vs. standard test on ISat_Dploy_Sm.dat"""
        bdf_filename = os.path.join(test_path, 'iSat', 'ISat_Dploy_Sm.dat')
        bdf_filename_out = os.path.join(test_path, 'iSat', 'ISat_Dploy_Sm2.dat')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    @staticmethod
    def test_isat_02():
        """vectorized vs. standard test on ISat_Launch_Sm_4pt.dat"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_4pt.dat')
        bdf_filename_outv = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_4ptv.dat')
        bdf_filename_out = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_4pt2.dat')

        vmodel = read_bdfv(bdf_filename)
        vmodel.write_bdf(bdf_filename_outv)
        model = read_bdf(bdf_filename, log=log)
        model.write_bdf(bdf_filename_out)

        run_and_compare_fems(
            bdf_filename, bdf_filename_outv, debug=False, xref=True, check=True,
            punch=False, mesh_form=None,
            print_stats=False, encoding=None,
            sum_load=False, size=8, is_double=False,
            stop=False, nastran='', post=-1, dynamic_vars=None,
            quiet=False, dumplines=False, dictsort=False,
            nerrors=0, dev=False, crash_cards=None,
        )
        os.remove(bdf_filename_out)

    @staticmethod
    def test_isat_03():
        """vectorized vs. standard test on ISat_Launch_Sm_Rgd.dat"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_Rgd.dat')
        #bdf_filename_outv = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_Rgdv.dat')
        bdf_filename_out = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_Rgd2.dat')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)

    def test_isat_04(self):
        """vectorized vs. standard test on iSat_launch_100Hz.dat"""
        bdf_filename = os.path.join(test_path, 'iSat', 'iSat_launch_100Hz.dat')
        #bdf_filename_outv = os.path.join(test_path, 'iSat', 'ISat_Launch_Sm_Rgdv.dat')
        bdf_filename_out = os.path.join(test_path, 'iSat', 'iSat_launch_100Hz2.dat')

        read_write_compare(bdf_filename, bdf_filename_out)
        os.remove(bdf_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
