import os
import unittest
from io import StringIO
from numpy import allclose, array
from cpylog import get_logger

import pyNastran
from pyNastran.utils import object_attributes, object_methods
#from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf import BDF, read_bdf, CrossReferenceError
from pyNastran.bdf.write_path import write_include, _split_path
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
from pyNastran.bdf.test.test_bdf import run_bdf, compare, run_lots_of_files, main as test_bdf

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'test')
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

class Tester(unittest.TestCase):

    def run_bdf(self, folder, bdf_filename, xref=False, size=8,
                mesh_form='combined', dynamic_vars=None, debug=False, quiet=True,
                run_extract_bodies=True,
                run_skin_solids=True, save_file_structure=False, log=None):
        #xref = False
        if quiet:
            debug = None
        return run_bdf(folder, bdf_filename, xref=xref, size=size,
                       is_folder=True,
                       mesh_form=mesh_form, dynamic_vars=dynamic_vars,
                       debug=debug, quiet=quiet,
                       sum_load=True, run_extract_bodies=run_extract_bodies,
                       run_skin_solids=run_skin_solids,
                       save_file_structure=save_file_structure, log=log)


class TestBDF(Tester):

    def test_bdf_test(self):
        filenames = [
            #os.path.join(MODEL_PATH, 'superelements', 'resvec23.bdf'),
            os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf'),
        ]
        run_lots_of_files(filenames, folder='', debug=False, xref=True, check=True, punch=False,
                          nastran='', encoding=None, size=None, is_double=None, post=None,
                          sum_load=True, dev=True, crash_cards=None, pickle_obj=True,
                          write_hdf5=True, quiet=False)


    def test_write_path(self):
        include_name = r'C:\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'
        msg1 = write_include(include_name, is_windows=True)
        sline1 = _split_path(include_name, is_windows=True)

        include_name = r'/opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
        msg2 = write_include(include_name, is_windows=False)
        sline2 = _split_path(include_name, is_windows=False)

        include_name = r'/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
        msg3 = write_include(include_name, is_windows=False)
        sline3 = _split_path(include_name, is_windows=False)

        include_name = r'opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'
        msg4 = write_include(include_name, is_windows=True)
        sline4 = _split_path(include_name, is_windows=True)

        msg1_expected = r"INCLUDE 'C:\\NASA\formats\pynastran_v0.6\pyNastran\bdf\writePath.py'" + '\n'
        msg2_expected = "INCLUDE '/opt/NASA/formats/pynastran_v0.6/pyNastran/bdf/writePath.py'\n"
        msg3_expected = ("INCLUDE '/opt/NASA/test1/test2/test3/test4/formats/pynastran_v0.6/\n"
                         "        pyNastran/bdf/writePath.py'\n")
        msg4_expected = (r"INCLUDE 'opt\NASA\test1\test2\test3\test4\formats\pynastran_v0.6" + '\\\n' +
                         r"        pyNastran\bdf\writePath.py'" + '\n')
        assert msg1 == msg1_expected, 'test1 actual:\n%r\nexpected:\n%r\n%s' % (msg1, msg1_expected, str(sline1))
        assert msg2 == msg2_expected, 'test2 actual:\n%r\nexpected:\n%r\n%s' % (msg2, msg2_expected, str(sline2))
        assert msg3 == msg3_expected, 'test3 actual:\n%r\nexpected:\n%r\n%s' % (msg3, msg3_expected, str(sline3))
        assert msg4 == msg4_expected, 'test4 actual:\n%s\nexpected:\n%s\n%s' % (msg4, msg4_expected, str(sline4))

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
        log = get_logger(log=None, level='error', encoding='utf-8')
        self.run_bdf('', bdf_filename, log=log)
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, xref=True, log=log)
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
        mass, cg, unused_I = mass_properties(fem1)

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
        cards = model._read_bdf_cards(bdf_filename,
                                      punch=False, read_includes=True, encoding=None)
        assert len(cards) == 9, len(cards)


        etype_to_eids_pids_nids = model.get_elements_properties_nodes_by_element_type()
        assert len(etype_to_eids_pids_nids) == 1, list(etype_to_eids_pids_nids.keys())
        #etype_to_eids_pids_nids[etype] : [eids, pids, nids]
        unused_eids, pids, unused_node_ids = etype_to_eids_pids_nids['CTETRA4']
        assert pids.min() == 1, pids.min()
        assert pids.max() == 1, pids.max()

        out = model.get_elements_nodes_by_property_type(
            dtype='int32', save_element_types=False)
        etype_pid_to_eids_nids, _etype_to_eids_pids_nids = out
        #etype_pid_to_eids_nids[(etype, pid)] : [eids, nids]
        assert _etype_to_eids_pids_nids is None, _etype_to_eids_pids_nids
        assert len(etype_pid_to_eids_nids) == 1, list(etype_pid_to_eids_nids.keys())

    def test_bdf_02(self):
        """checks plate_py.dat"""
        bdf_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.dat')
        self.run_bdf('', bdf_filename)
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, xref=True)

        etype_to_eids_pids_nids = fem1.get_elements_properties_nodes_by_element_type()
        assert len(etype_to_eids_pids_nids) == 1, list(etype_to_eids_pids_nids.keys())
        #etype_to_eids_pids_nids[etype] : [eids, pids, nids]
        unused_eids, pids, unused_node_ids = etype_to_eids_pids_nids['CQUAD4']
        assert pids.min() == 1, pids.min()
        assert pids.max() == 1, pids.max()

        etype_pid_to_eids_nids, _etype_to_eids_pids_nids = fem1.get_elements_nodes_by_property_type(
            dtype='int32', save_element_types=False)
        #etype_pid_to_eids_nids[(etype, pid)] : [eids, nids]
        assert _etype_to_eids_pids_nids is None, _etype_to_eids_pids_nids
        assert len(etype_pid_to_eids_nids) == 1, list(etype_pid_to_eids_nids.keys())

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
        bdf_filename = os.path.join(MODEL_PATH, 'unit', 'cbush', 'cbush.dat')
        log = get_logger(log=None, level='error', encoding='utf-8')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)

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

        self.run_bdf('', bdf_filename, xref=True, debug=False, log=log)

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
        unused_model = read_bdf(bdf_file, validate=True, xref=True,
                                punch=False, skip_cards=None,
                                read_cards=None,
                                encoding=None, log=None,
                                debug=True, mode='msc')

    def test_bdf_xref_fail(self):
        """tests various xref's failing"""
        model = BDF(debug=False, log=None, mode='msc')
        def _run(model, bdf_filename):
            """helper for ``test_bdf_xref_fail``"""
            bdf_filename.seek(0)
            with self.assertRaises(CrossReferenceError):
                model.read_bdf(bdf_filename)
            model.clear_attributes()
            model.safe_cross_reference()
            model.uncross_reference()
            model.safe_cross_reference(xref=False)
            model.clear_attributes()
        #-------------------------------------------------
        # missing node_id on element
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1\n"
            "CONROD,10,1,2,1000,\n"
            "MAT1,1000,3.0e7,,0.3\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing coord_id on grid
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1,2\n"
            "GRID,2\n"
            "CONROD,10,1,2,1000,\n"
            "MAT1,1000,3.0e7,,0.3\n"
        )
        #_run(model, bdf_filename)
        bdf_filename.seek(0)
        with self.assertRaises(KeyError):
            model.read_bdf(bdf_filename)
        model.clear_attributes()
        model.safe_cross_reference()
        model.clear_attributes()

        #-------------------------------------------------
        # missing material on property
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1\n"
            "GRID,2\n"
            "CROD,10,100,1,2\n"
            "PROD,100,1000,1.0\n"
            #"MAT1,1000,3.0e7,,0.3\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing node_id on conm2
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            #"GRID,1\n"
            "CONM2,10,1\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing node_id on plotel
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1\n"
            #"GRID,2\n"
            "PLOTEL,10,1,2\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing load on LOAD
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1\n"
            "GRID,2\n"
            "CONROD,10,1,2,1000,\n"
            "MAT1,1000,3.0e7,,0.3\n"
            "LOAD,1,1.0,1.0,10000\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing load on LOAD
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "GRID,1\n"
            "GRID,2\n"
            "CONROD,10,1,2,1000,\n"
            "MAT1,1000,3.0e7,,0.3\n"
            "LOAD,1,1.0,1.0,10000\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing element_id on PLOAD4
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            #"GRID,1\n"
            #"GRID,2\n"
            #"GRID,3\n"
            #"GRID,4\n"
            #"CQUAD4,10,1000,1,2,3,4\n"
            "MAT1,1000,3.0e7,,0.3\n"
            "LOAD,1,1.0,1.0,10000\n"
            "PLOAD4,10000,10\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing node_id on DPHASE
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "DPHASE,200, 1,4\n"
            #"GRID,1\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing TLOAD1/RLOAD1 on DLOAD
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "DLOAD,1,1.0,1.0,10000\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing node_id on RBE2
        bdf_filename = StringIO()
        bdf_filename.write(
            "CEND\n"
            "BEGIN BULK\n"
            "RBE2,10,1,123456\n"
        )
        _run(model, bdf_filename)

        #-------------------------------------------------
        # missing node_id on TIC
        bdf_filename = StringIO()
        bdf_filename.write(
            'CEND\n'
            'BEGIN BULK\n'
            'TIC,14,1,-1\n'
            'TIC,10,1,0\n'
            'TIC,11,1,1\n'
            'TIC,12,1,2\n'
            'TIC,13,1,3\n'
            'TIC,14,1,4\n'
            'TIC,15,1,5\n'
            'TIC,16,1,6\n'
        )
        #_run(model, bdf_filename)
        bdf_filename.seek(0)
        with self.assertRaises(CrossReferenceError):
            model.read_bdf(bdf_filename)
        model.pop_parse_errors()
        model.clear_attributes()
        model.safe_cross_reference()
        model.uncross_reference()
        model.safe_cross_reference(xref=False)
        model.clear_attributes()

        with self.assertRaises(SyntaxError):
            model.add_card(['TIC', 15, -2, -5], 'TIC', comment='tic', is_list=True, has_none=True)
        #with self.assertRaises(SyntaxError):
            #model.pop_parse_errors()

    def test_bdf_xref_safe(self):
        """testing various safe_xref methods"""
        model = BDF(debug=False, log=None, mode='msc')
        aefact_id = 2
        caero_id = 1
        ref_id = caero_id
        xref_errors = {'aefact' : []}
        assert model.safe_aefact(aefact_id, ref_id, xref_errors, msg='') is None
        assert xref_errors == {'aefact' : [(1, 2)]}, xref_errors

        #--------------------------------------------------------
        cid = 2
        eid = 1
        ref_id = eid
        xref_errors = {'cid' : []}
        assert model.safe_coord(cid, ref_id, xref_errors, msg='') is None
        assert xref_errors == {'cid' : [(1, 2)]}, xref_errors

        #--------------------------------------------------------
        mid = 2
        pid = 1
        ref_id = eid
        xref_errors = {'mid' : []}
        assert model.safe_material(mid, ref_id, xref_errors, msg='') is None
        assert xref_errors == {'mid' : [(1, 2)]}, xref_errors

        #--------------------------------------------------------
        pid = 2
        eid = 1
        ref_id = eid
        xref_errors = {'pid' : []}
        assert model.safe_property(pid, ref_id, xref_errors, msg='') is None
        assert xref_errors == {'pid' : [(1, 2)]}, xref_errors

        #--------------------------------------------------------
        pid = 2
        eid = 1
        ref_id = eid
        xref_errors = {'pid' : []}
        assert model.safe_property_mass(pid, ref_id, xref_errors, msg='') is None
        assert xref_errors == {'pid' : [(1, 2)]}, xref_errors

        #--------------------------------------------------------
        point_ids = [1, 2, 3]
        eid = 1
        ref_id = eid
        #xref_errors = {'pid' : []}
        assert model.safe_get_points(point_ids, msg='')[0] == point_ids, model.safe_get_points(point_ids, msg='')
        #assert xref_errors == {'pid' : [(1, 2)]}, xref_errors

    def test_bdf_05(self):
        """checks testA.dat"""
        bdf_filename = os.path.join(PKG_PATH, 'bdf', 'test', 'unit', 'testA.bdf')
        (unused_fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=False, run_extract_bodies=False,
            #save_file_structure=True,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        if len(diff_cards2) != 0:  # pragma: no cover
            msg = 'check testA.test_bdf.out\ndiff_cards2=%s\n' % (diff_cards2)
            raise AssertionError(msg)

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

    def test_bdf_superelement_1(self):
        """checks resvec23.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'resvec23.bdf')
        log = get_logger(log=None, level='error', encoding='utf-8')
        (unused_fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=True, run_extract_bodies=False,
            save_file_structure=False, log=log,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2
        #bdf_filenames = {bdf_filename : 'cat.bdf',}
        #fem1.write_bdfs(bdf_filenames)
        #os.remove('cat.bdf')

    def test_bdf_superelement_2(self):
        """checks superelement.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'superelement.bdf')
        log = get_logger(log=None, level='error', encoding='utf-8')
        (fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=True, run_extract_bodies=False,
            save_file_structure=True, log=log,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2
        bdf_filenames = {bdf_filename : 'cat.bdf',}
        fem1.write_bdfs(bdf_filenames)
        os.remove('cat.bdf')

        #os.remove(bdf_filename + '_out')
        #self.run_bdf(folder, bdf_filename, xref=True) # PBEAML is not supported

    def test_bdf_superelement_3(self):
        """checks cqrsee101b2.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'cqrsee101b2.bdf')
        log = get_logger(log=None, level='error', encoding='utf-8')
        (fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=True, run_extract_bodies=False,
            save_file_structure=True, log=log,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2
        bdf_filenames = {bdf_filename : 'cat.bdf',}
        fem1.write_bdfs(bdf_filenames)
        os.remove('cat.bdf')

    def test_bdf_superelement_4(self):
        """checks see101l8.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'see101l8.bdf')
        log = get_logger(log=None, level='error', encoding='utf-8')
        (unused_fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=True, run_extract_bodies=False,
            save_file_structure=True, log=log,
        )
        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2
        #bdf_filenames = {bdf_filename : 'cat.bdf',}
        #fem1.write_bdfs(bdf_filenames)
        #os.remove('cat.bdf')

    def test_bdf_superelement_5(self):
        """checks flyswatter.bdf"""
        from pyNastran.bdf.mesh_utils.bdf_renumber import superelement_renumber
        model_path = os.path.join(MODEL_PATH, 'superelements', 'flyswatter')
        bdf_filename = os.path.join(model_path, 'flyswatter.bdf')
        bdf_filename_out = os.path.join(model_path, 'flyswatter.re.bdf')
        #log = get_logger(log=None, level='error', encoding='utf-8')

        fem1 = read_bdf(bdf_filename, validate=True, xref=True, punch=False,
                        save_file_structure=False, skip_cards=None, read_cards=None,
                        encoding=None, log=None, debug=True, mode='msc')

        superelement_renumber(
            fem1, bdf_filename_out=bdf_filename_out,
            starting_id_dict=None)

    def test_bdf_other_1(self):
        """checks axisymmetric model"""
        bdf_filename = os.path.join(MODEL_PATH, 'other', 'd07d2.bdf')
        bdf_filename_test = os.path.join(MODEL_PATH, 'other', 'd07d2.test_bdf.bdf')
        fem1 = read_bdf(bdf_filename, validate=True, xref=True, punch=False,
                        skip_cards=None, read_cards=None,
                        encoding=None, log=None, debug=False, mode='msc')
        fem1.write_bdf(bdf_filename_test)
        fem2 = read_bdf(bdf_filename_test, debug=None)

        diff_cards = compare(fem1, fem2, xref=True, check=False,
                             print_stats=True, quiet=True)

        diff_cards2 = list(set(diff_cards))
        diff_cards2.sort()
        assert len(diff_cards2) == 0, diff_cards2

        #os.remove(bdf_filename + '_out')
        #self.run_bdf(folder, bdf_filename, xref=True) # PBEAML is not supported

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

    def test_aero_02(self):
        """checks 0012_flutter.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.bdf')
        #log = get_logger(log=None, level='error', encoding='utf-8')
        argv = ['test_bdf', bdf_filename, '-q']
        test_bdf(argv=argv)
        #self.run_bdf('', bdf_filename, log=log)
        #fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, xref=True, log=log)
        #diff_cards2 = list(set(diff_cards))
        #diff_cards2.sort()
        #assert len(diff_cards2) == 0, diff_cards2

def compare_mass_cg_inertia(fem1, reference_point=None, sym_axis=None):
    unused_mass1, unused_cg1, unused_I1 = mass_properties(
        fem1, reference_point=reference_point, sym_axis=sym_axis)
    #mass1, cg1, I1 = fem1.mass_properties_no_xref(reference_point=reference_point, sym_axis=sym_axis)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
