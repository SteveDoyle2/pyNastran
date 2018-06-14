from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
from numpy import allclose, array
from six import StringIO

import pyNastran
from pyNastran.utils import object_attributes, object_methods
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.write_path import write_include, _split_path
from pyNastran.bdf.test.test_bdf import run_bdf, run_all_files_in_folder

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'test')
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

class Tester(unittest.TestCase):

    def run_bdf(self, folder, bdf_filename, xref=False, cid=None, size=8,
                mesh_form='combined', dynamic_vars=None, debug=False, quiet=True,
                run_extract_bodies=True):
        cid = 0
        #xref = False
        return run_bdf(folder, bdf_filename, xref=xref, cid=cid, size=size,
                       is_folder=True,
                       mesh_form=mesh_form, dynamic_vars=dynamic_vars,
                       debug=debug, quiet=quiet,
                       sum_load=True, run_extract_bodies=run_extract_bodies)

    def run_all_files_in_folder(self, folder, xref=False, cid=None, debug=False):
        run_all_files_in_folder(folder, xref=xref, cid=cid, debug=debug)

is_windows = 'nt' in os.name

class TestBDF(Tester):

    @unittest.skipIf(not is_windows, 'write_include doesnt writing INCLUDEs on mac/linux')
    def test_write_path(self):
        if is_windows:
            include_name = r'C:\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'
            msg1 = write_include(include_name, is_windows=True)
            sline1 = _split_path(include_name)

            include_name = r'/opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
            msg2 = write_include(include_name, is_windows=False)
            sline2 = _split_path(include_name)

            include_name = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
            msg3 = write_include(include_name, is_windows=False)
            sline3 = _split_path(include_name)

            include_name = r'opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
            msg4 = write_include(include_name, is_windows=True)
            sline4 = _split_path(include_name)

            msg1_expected = r'INCLUDE C:\\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py' '\n'
            msg2_expected = 'INCLUDE /opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py\n'
            msg3_expected = ('INCLUDE /opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/\n'
                             '        pyNastran/bdf/writePath.py\n')
            msg4_expected = (r'INCLUDE opt\NASA\test1\test2\test3\test4\formats\pynastran_v0.6' '\\\n'
                             r'        pyNastran\bdf\writePath.py' '\n')
            assert msg1 == msg1_expected, 'test1 actual:\n%r\nexpected:\n%r\n%s' % (msg1, msg1_expected, str(sline1))
            assert msg2 == msg2_expected, 'test2 actual:\n%r\nexpected:\n%r\n%s' % (msg2, msg2_expected, str(sline2))
            assert msg3 == msg3_expected, 'test3 actual:\n%r\nexpected:\n%r\n%s' % (msg3, msg3_expected, str(sline3))
            assert msg4 == msg4_expected, 'test4 actual:\n%r\nexpected:\n%r\n%s' % (msg4, msg4_expected, str(sline4))

    def test_object_attributes_01(self):
        """tests getting object attributes"""
        model = BDF(debug=False)
        model.object_attributes(mode='public', keys_to_skip=None)

    def test_object_attributes_02(self):
        """tests getting object attributes with key skipping"""
        model = BDF(debug=False)
        keys = []
        object_attributes(model, mode='public', keys_to_skip=keys)

    def test_object_attributes_03(self):
        """tests getting object attributes with a card"""
        model = BDF(debug=False)
        model.add_card(['GRID', 1], 'GRID')
        grid = model.nodes[1]
        grid.object_attributes(mode='public', keys_to_skip=None)

    def test_object_methods_01(self):
        """tests getting object methods using the builtin method"""
        model = BDF(debug=False)
        keys = []
        model.object_methods(mode="public", keys_to_skip=keys)

    def test_object_methods_02(self):
        """tests getting object methods from the method"""
        model = BDF(debug=False)
        keys = []
        object_methods(model, mode="public", keys_to_skip=keys)

    def test_object_methods_03(self):
        """tests getting object attributes with a card"""
        model = BDF(debug=False)
        model.add_card(['GRID', 1], 'GRID')
        grid = model.nodes[1]
        grid.object_methods(mode='public', keys_to_skip=None)

    def test_bdf_01(self):
        """checks solid_bending.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        self.run_bdf('', bdf_filename)
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, xref=True)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 2, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 72, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 186, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        mass, cg, I = fem1.mass_properties()

        assert allclose(mass, 6.0), 'mass = %s' % mass
        cg_exact = array([0.5, 1., 1.5])
        for i, (cgi, cgie) in enumerate(zip(cg, cg_exact)):
            assert allclose(cgi, cgie), 'i=%s cg=%s' % (i, str(cg))

        compare_mass_cg_inertia(fem1)
        compare_mass_cg_inertia(fem1, reference_point=None)

    def test_bdf_01_hdf5(self):
        """checks solid_bending.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        model = BDF(debug=False)
        cards = model._read_bdf_cards(bdf_filename, validate=True, xref=False,
                                      punch=False, read_includes=True, encoding=None)
        assert len(cards) == 9, len(cards)

    def test_bdf_02(self):
        """checks plate_py.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.dat')
        self.run_bdf('', bdf_filename)
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, xref=True)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.coords) == 3, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.nodes) == 231, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 200, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        compare_mass_cg_inertia(fem1)

    def test_bdf_03(self):
        """checks cbush.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'cbush', 'cbush.dat')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, debug=False)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 5, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 2, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 0, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 1, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue

        compare_mass_cg_inertia(fem1)

        self.run_bdf('', bdf_filename, xref=True, debug=False)

    def test_bdf_04(self):
        """checks beam_modes.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 12, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 10, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.masses) == 1, 'len(masses) = %i' % len(fem.masses)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 3, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue
            assert len(fem.properties_mass) == 0, 'len(properties_mass) = %i' % len(fem.properties_mass)
        fem1.cross_reference()
        compare_mass_cg_inertia(fem1)

        #self.run_bdf(folder, bdf_filename, xref=True, debug=False) # PBEAML is not supported

    def _test_bdf_slash(self):
        """tests a / in a deck"""
        lines = [
            '$ DEC/CMS REPLACEMENT HISTORY, Element D10912R.DAT',
            '$ *1    15-JUN-1990 17:41:35 CMSMGR "66B PLUS/G 66B/ Initial installation of TPL test problems"',
            '$ DEC/CMS REPLACEMENT HISTORY, Element D10912R.DAT',
            'RESTART VERSION=LAST,KEEP $ RESTART FROM D10911R DBS=D10911D',
            'ID EDS, D10912R  $',
            '$ID EDS, D2712R   $',
            '$ID EDS,D27D12R',
            'SOL 109 $',
            '$SOL 27,0',
            '$DIAG 8,14',
            'TIME 5',
            '$READ 10 $ D27D11',
            'CEND',
            'TITLE=NEW RIGID FORMATS - CANTILEVER BEAM                      D10912R',
            'SUBTITLE=DIRECT TRANSIENT',
            'SET 1000=10,30,40',
            'SET 2000=111,200',
            'METHOD=1',
            'DISP(SORT2)=2000',
            'SPC=200',
            'MPC=100',
            'TSTEP=100',
            'FORCE=1000',
            'SUBCASE 1',
            'DLOAD=10',
            'BEGIN BULK',
            '/       47',
            'TLOAD1,10,2,0,0,10',
            'ENDDATA',
        ]
        bdf_file = StringIO()
        bdf_file.writelines(lines)
        bdf_file.seek(0)
        #with self.assertRaises(NotImplementedError):
        model = read_bdf(bdf_file, validate=True, xref=True,
                         punch=False, skip_cards=None,
                         read_cards=None,
                         encoding=None, log=None,
                         debug=True, mode='msc')


    def test_bdf_05(self):
        """checks testA.dat"""
        bdf_filename = os.path.join(PKG_PATH, 'bdf', 'test', 'unit', 'testA.bdf')
        (fem1, fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=False, run_extract_bodies=False,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2
        #os.remove(bdf_filename + '_out')
        #self.run_bdf(folder, bdf_filename, xref=True) # PBEAML is not supported

    def test_bdf_06(self):
        """checks bar3truss/vared_bar3.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'bar3truss', 'vared_bar3.bdf')

        dynamic_vars = {
            'bar1_a': 1.0,
            'bar2_a': 1.0,
            'bar3_a': 1.0,
            'loadx': 1.0,
            'loady': 1.0,
            'loadmag': 10000.,
            'youngs' : 1e7,
            'rho': 0.01,
        }
        fem1, fem2, diff_cards = self.run_bdf(
            '', bdf_filename, dynamic_vars=dynamic_vars)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 4, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 4, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 3, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 3, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue

        fem1.cross_reference()
        fem1.pop_xref_errors()
        compare_mass_cg_inertia(fem1)

    def test_bdf_thermal_01(self):
        """checks time_thermal_elements.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert len(fem.params) == 1, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 9, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 0, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 7, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.masses) == 0, 'len(masses) = %i' % len(fem.masses)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
            assert len(fem.properties_mass) == 0, 'len(properties_mass) = %i' % len(fem.properties_mass)
        fem1.cross_reference()
        compare_mass_cg_inertia(fem1)

    def test_bdf_transfer_function_01(self):
        """checks transfer_function/actuator_tf_modeling.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'transfer_function', 'actuator_tf_modeling.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert fem.card_count['CONM2'] == 3, fem.card_count
            assert fem.card_count['SPOINT'] == 1, fem.card_count
            assert fem.card_count['EPOINT'] == 1, fem.card_count
            assert fem.card_count['PARAM'] == 1, fem.card_count
            assert fem.card_count['CELAS2'] == 2, fem.card_count
            assert fem.card_count['GRID'] == 3, fem.card_count
            assert fem.card_count['EIGR'] == 1, fem.card_count
            assert fem.card_count['EIGC'] == 1, fem.card_count
            assert fem.card_count['MPC'] == 1, fem.card_count
            assert fem.card_count['TF'] == 2, fem.card_count

    def test_bdf_aero_01(self):
        """checks aero/aerobeam.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'aerobeam.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename)
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        for fem in [fem1, fem2]:
            assert fem.card_count['MAT1'] == 3, fem.card_count
            assert fem.card_count['DCONADD'] == 2, fem.card_count
            assert fem.card_count['MAT1'] == 3, fem.card_count
            assert fem.card_count['DMI'] == 6, fem.card_count
            assert fem.card_count['PAERO1'] == 1, fem.card_count
            assert fem.card_count['EIGRL'] == 2, fem.card_count
            assert fem.card_count['PBAR'] == 1, fem.card_count
            assert fem.card_count['DESVAR'] == 3, fem.card_count
            assert fem.card_count['DRESP1'] == 11, fem.card_count
            assert fem.card_count['DRESP2'] == 6, fem.card_count

            assert fem.card_count['SPC1'] == 4, fem.card_count
            assert fem.card_count['AESTAT'] == 10, fem.card_count
            assert fem.card_count['TRIM'] == 4, fem.card_count
            assert fem.card_count['SPLINE2'] == 3, fem.card_count
            assert fem.card_count['DVPREL1'] == 6, fem.card_count
            assert fem.card_count['SUPORT1'] == 2, fem.card_count
            assert fem.card_count['DCONSTR'] == 10, fem.card_count
            assert fem.card_count['AELIST'] == 3, fem.card_count
            assert fem.card_count['CORD2R'] == 6, fem.card_count
            assert fem.card_count['CONM2'] == 10, fem.card_count

            assert fem.card_count['ENDDATA'] == 1, fem.card_count
            assert fem.card_count['AERO'] == 1, fem.card_count
            assert fem.card_count['PARAM'] == 4, fem.card_count
            assert fem.card_count['CBEAM'] == 3, fem.card_count
            assert fem.card_count['GRID'] == 14, fem.card_count
            assert fem.card_count['SET1'] == 5, fem.card_count
            assert fem.card_count['MKAERO1'] == 1, fem.card_count
            assert fem.card_count['PBEAML'] == 3, fem.card_count
            assert fem.card_count['FLFACT'] == 5, fem.card_count
            assert fem.card_count['AESURF'] == 3, fem.card_count
            assert fem.card_count['DEQATN'] == 3, fem.card_count
            assert fem.card_count['CBAR'] == 4, fem.card_count
            assert fem.card_count['CAERO1'] == 3, fem.card_count
            assert fem.card_count['AEROS'] == 1, fem.card_count
            assert fem.card_count['FLUTTER'] == 4, fem.card_count
            assert fem.card_count['DOPTPRM'] == 1, fem.card_count

def compare_mass_cg_inertia(fem1, reference_point=None, sym_axis=None):
    mass1, cg1, I1 = fem1.mass_properties(reference_point=reference_point, sym_axis=sym_axis)
    #mass1, cg1, I1 = fem1.mass_properties_no_xref(reference_point=reference_point, sym_axis=sym_axis)


class TestBaseCard(Tester):
    """Tests methods used by ``BaseCard``"""
    def test_base_card_01_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        msg = 'expected=%s actual=%s' % (expected, collapse_thru_by(data))
        self.assertEqual(collapse_thru_by(data), expected, msg)

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 7, 9, 11, 12, 14, 16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2]
        expected = [1, 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 3, 5, 7, 9, 11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4]
        expected = [1, 'THRU', 4]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3]
        expected = [1, 2, 3]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEqual(collapse_thru_by(data), expected, collapse_thru_by(data))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
