import os
import unittest
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.nastranIOv import NastranIO
import pyNastran
#from pyNastran.utils.log import get_logger2

class NastranGUI(NastranIO, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)

pkg_path = pyNastran.__path__[0]
model_path = os.path.join(pkg_path, '..', 'models')



class TestNastranGUI(unittest.TestCase):

    def test_solid_shell_bar_01(self):
        bdf_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        op2_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_solid_shell_bar_02(self):
        bdf_filename = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.bdf')
        op2_filename = os.path.join(model_path, 'sol_101_elements', 'mode_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_solid_shell_bar_03(self):
        bdf_filename = os.path.join(model_path, 'sol_101_elements', 'buckling_solid_shell_bar.bdf')
        op2_filename = os.path.join(model_path, 'sol_101_elements', 'buckling_solid_shell_bar.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_solid_bending(self):
        bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
        op2_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_beam_modes_01(self):
        bdf_filename = os.path.join(model_path, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(model_path, 'beam_modes', 'beam_modes_m1.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_beam_modes_02(self):
        bdf_filename = os.path.join(model_path, 'beam_modes', 'beam_modes.dat')
        op2_filename = os.path.join(model_path, 'beam_modes', 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_beam_modes_03(self):
        dirname = os.path.join(model_path, 'beam_modes')
        bdf_filename = os.path.join(dirname, 'beam_modes.dat')
        op2_filename = os.path.join(dirname, 'beam_modes_m1.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        #test.load_nastran_results(op2_filename, None)

        test.load_nastran_geometry(bdf_filename, dirname)
        #test.load_nastran_results(op2_filename, dirname)

        test.load_nastran_geometry(bdf_filename, '')
        test.load_nastran_results(op2_filename, dirname)

    def test_beam_modes_04(self):
        dirname = os.path.join(model_path, 'beam_modes')
        bdf_filename = os.path.join(dirname, 'beam_modes.dat')
        op2_filename = os.path.join(dirname, 'beam_modes_m2.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

        test.load_nastran_geometry(bdf_filename, dirname)
        test.load_nastran_results(op2_filename, dirname)

        test.load_nastran_geometry(bdf_filename, '')


    @unittest.expectedFailure
    def test_contact(self):
        """this test fails because of a misparsed card"""
        bdf_filename = os.path.join(model_path, 'contact', 'contact.bdf')
        op2_filename = os.path.join(model_path, 'contact', 'contact.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_fsi(self):
        bdf_filename = os.path.join(model_path, 'fsi', 'fsi.bdf')
        op2_filename = os.path.join(model_path, 'fsi', 'fsi.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_thermal_01(self):
        dirname = os.path.join(model_path, 'thermal')
        bdf_filename = os.path.join(dirname, 'thermal_test_153.bdf')
        op2_filename = os.path.join(dirname, 'thermal_test_153.op2')

        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_bwb_gui(self):
        bdf_filename = os.path.join(model_path, 'bwb', 'BWB_saero.bdf')
        test = NastranGUI()
        #test.log = get_logger2()
        test.load_nastran_geometry(bdf_filename, None)

    def test_femap_rougv1_01(self):
        """tests the exhaust manifold and it's funny eigenvectors"""
        dirname = os.path.join(model_path, 'femap_exhaust')
        #bdf_filename = os.path.join(dirname, 'modal_example.bdf')
        op2_filename = os.path.join(dirname, 'modal_example.op2')

        test = NastranGUI()
        test.load_nastran_geometry(op2_filename, None)
        test.load_nastran_results(op2_filename, None)

    def test_aero(self):
        """tests the bah_plane"""
        bdf_filename = os.path.join(model_path, 'aero', 'bah_plane', 'bah_plane.bdf')
        test = NastranGUI()
        test.load_nastran_geometry(bdf_filename, None)


def test_bottle():  # pragma: no cover
    """
    Tests Nastran GUI loading
    """
    test = NastranGUI()
    test.load_nastran_geometry('bottle_shell_w_holes_pmc.bdf', '')
    test.load_nastran_results('bottle_shell_w_holes_pmc.op2', '')

    keys = test.result_cases.keys()
    assert (1, 'Stress1', 1, 'centroid', '%.3f') in keys, keys

if __name__ == '__main__':  # pragma: no cover
    #test_solid_shell_bar_01()
    unittest.main()

