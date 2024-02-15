import os
import unittest
from io import StringIO
import numpy as np
from numpy import allclose, array
from cpylog import SimpleLogger

import pyNastran
from pyNastran.utils import object_attributes, object_methods
#from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf import BDF, read_bdf, CrossReferenceError
from pyNastran.bdf.write_path import write_include, _split_path
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
from pyNastran.bdf.mesh_utils.forces_moments import get_forces_moments_array
from pyNastran.bdf.test.test_bdf import run_bdf, compare, run_lots_of_files, main as test_bdf

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'test')
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
GUI_MODEL_DIRNAME = os.path.join(PKG_PATH, 'converters', 'nastran', 'models')

class Tester(unittest.TestCase):

    def run_bdf(self, folder: str, bdf_filename: str,
                xref=False, size=8,
                mesh_form='combined', dynamic_vars=None, debug=False, quiet=True,
                run_extract_bodies=True,
                run_skin_solids=True, save_file_structure=False,
                validate_case_control=True, log=None):
        #xref = False
        if quiet:
            debug = None
        return run_bdf(folder, bdf_filename, xref=xref, size=size,
                       is_folder=True,
                       mesh_form=mesh_form, dynamic_vars=dynamic_vars,
                       debug=debug, quiet=quiet,
                       sum_load=True, run_extract_bodies=run_extract_bodies,
                       run_skin_solids=run_skin_solids,
                       save_file_structure=save_file_structure,
                       validate_case_control=validate_case_control,
                       log=log)


class TestBDFUnit(Tester):

    def test_bdf_pickle_copy(self):
        """verify we get 5 include files if they are one after the other"""
        model = BDF(debug=False)
        #from copy import deepcopy
        #from pyNastran.bdf.bdf import read_bdf
        #from pyNastran.bdf.mesh_utils.mass_properties import mass_properties

        dirname = os.path.join(MODEL_PATH, 'bugs', 'euler_column_linear_buckling')

        bdf_filename = os.path.join(dirname, 'euler_column_linear_buckling.bdf')
        obj_filename = os.path.join(dirname, 'model.obj')
        model = read_bdf(bdf_filename)
        model.log.debug('model')
        model2 = model.__deepcopy__({})
        model2.log.debug('model2')

        #model2.cross_reference()
        mass = mass_properties(model2)[0]  # fails because of missing cross-reference, comment to reach next line
        model3 = model2.__deepcopy__({})  # fails while deep-copying
        model3.log.debug('model3')

        #import pickle
        model3.save(obj_filename=obj_filename, unxref=True)

    def test_bdf_include5(self):
        """verify we get 5 include files if they are one after the other"""
        model = BDF(debug=False)
        unit_dir = bdf_filename = os.path.join(TEST_PATH, 'unit')
        bdf_filename = os.path.join(unit_dir, 'include5.bdf')
        #bdf_filename_out = os.path.join(unit_dir, 'include5_out.bdf')
        model.read_bdf(bdf_filename, save_file_structure=False, read_includes=False)

        #print(model.include_filenames)
        file0_include_filenames = model.include_filenames[0]
        assert len(model.include_filenames) == 1, model.include_filenames
        assert len(file0_include_filenames) == 5, len(file0_include_filenames)
        #out_filenames = {}
        #for ifile, include_filenames in model.include_filenames.items():
            #for include_filename in include_filenames:
                #out_filenames[include_filename] = include_filename

        #model.write_bdfs(out_filenames, relative_dirname=unit_dir, encoding=None,
                         #size=8, is_double=False, enddata=None, close=True, is_windows=None)
        #model.write_bdf(bdf_filename_out, encoding=None, size=8,
                        #nodes_size=None, elements_size=None, loads_size=None,
                        #is_double=False, interspersed=False, enddata=None,
                        #write_header=True, close=True)
        x = 1

    def test_bdf_test(self):
        #log = SimpleLogger(level='warning', encoding='utf-8')
        filenames = [
            #os.path.join(MODEL_PATH, 'superelements', 'resvec23.bdf'),
            os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf'),
        ]
        run_lots_of_files(filenames, folder='', debug=None, xref=True, check=True, punch=False,
                          nastran='', encoding=None, size=None, is_double=None, post=None,
                          sum_load=True, dev=True, crash_cards=None, pickle_obj=True,
                          write_hdf5=True, quiet=True)

    def test_forces_moments_ctria3(self):
        bdf_filename = os.path.join(GUI_MODEL_DIRNAME, 'ctria3_pload_pload2_pload4.bdf')
        model = read_bdf(bdf_filename, debug=False)

        nnodes = 8
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3,
                   5: 4, 6: 5, 7: 6, 8: 7, }
        normals = np.array([
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
        ])
        eid_map = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5,}
        p0 = [0., 0., 0.]

        #-------------------------------
        # 3 load cases that are all the same (other than the card type)
        expected_pressures = [2., 2., 1., 1., 0., 0.]
        #expected_moments = []
        #expected_spcd = []

        # way harder to prove than the cquad4 version...we're kinda relying on the
        # PLOAD vs. PLOAD2 vs. PLOAD4 comparison
        expected_forces = array([
            [0., 0., 0.6666667 ],
            [0., 0., 0.6666667 ],
            [0., 0., 0.16666667],
            [0., 0., 0.        ],
            [0., 0., 0.33333334],
            [0., 0., 0.8333334 ],
            [0., 0., 0.33333334],
            [0., 0., 0.        ]], dtype='float32')
        #-------------------------------
        #  PLOAD
        load_case_id = 10000
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        assert is_loads
        assert len(out) == 4, out
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(centroidal_pressures, [0., 0., 0., 0., 0., 0.])
        assert np.allclose(expected_forces, forces)

        #-------------------------------
        #  PLOAD2
        load_case_id = 10002
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)

        #-------------------------------
        #  PLOAD2
        load_case_id = 10002
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)
        #-------------------------------
        #  PLOAD4
        load_case_id = 10004
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)
        x = 1

    def test_forces_moments_cquad4(self):
        bdf_filename = os.path.join(GUI_MODEL_DIRNAME, 'cquad4_pload_pload2_pload4.bdf')
        model = read_bdf(bdf_filename, debug=False)

        nnodes = 8
        nid_map = {1: 0, 2: 1, 3: 2, 4: 3,
                   5: 4, 6: 5, 7: 6, 8: 7, }
        normals = np.array([
            [0., 0., 1.],
            [0., 0., 1.],
            [0., 0., 1.],
        ])
        eid_map = {1: 0, 2: 1, 3: 2}
        p0 = [0., 0., 0.]

        #-------------------------------
        # 3 load cases that are all the same (other than the card type)
        expected_pressures = [2., 1., 0.]
        #expected_moments = []
        #expected_spcd = []
        expected_forces = array([
            [0.  , 0.  , 0.5 ],
            [0.  , 0.  , 0.75],
            [0.  , 0.  , 0.25],
            [0.  , 0.  , 0.  ],
            [0.  , 0.  , 0.5 ],
            [0.  , 0.  , 0.75],
            [0.  , 0.  , 0.25],
            [0.  , 0.  , 0.  ]], dtype='float32')
        #-------------------------------
        #  PLOAD
        load_case_id = 10000
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        assert is_loads
        assert len(out) == 4, out
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(centroidal_pressures, [0., 0., 0.])
        assert np.allclose(expected_forces, forces)

        #-------------------------------
        #  PLOAD2
        load_case_id = 10002
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)

        #-------------------------------
        #  PLOAD2
        load_case_id = 10002
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)
        #-------------------------------
        #  PLOAD4
        load_case_id = 10004
        dependents_nodes = set()
        is_loads, out = get_forces_moments_array(
            model, p0, load_case_id,
            eid_map, nnodes, normals, dependents_nodes, nid_map,
            include_grav=False, fdtype='float32')
        centroidal_pressures, forces, moments, spcd = out
        assert np.allclose(expected_pressures, centroidal_pressures)
        assert np.allclose(expected_forces, forces)
        x = 1

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
        log = SimpleLogger(level='error', encoding='utf-8')
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
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'plate_py', 'plate_py.dat')
        self.run_bdf('', bdf_filename, log=log)
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
        log = SimpleLogger(level='error', encoding='utf-8')
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
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'beam_modes', 'beam_modes.dat')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
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
        #self.run_bdf(folder, bdf_filename, xref=True, log=log) # PBEAML is not supported

    def test_bdf_06(self):
        """checks bar3truss/vared_bar3.bdf"""
        log = SimpleLogger(level='warning', encoding='utf-8')
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
            '', bdf_filename, dynamic_vars=dynamic_vars, log=log)
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
        log = SimpleLogger(level='error', encoding='utf-8')
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
        log = SimpleLogger(level='error', encoding='utf-8')
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
        #self.run_bdf(folder, bdf_filename, xref=True, log=log) # PBEAML is not supported

    def test_bdf_superelement_3(self):
        """checks cqrsee101b2.bdf"""
        bdf_filename = os.path.join(MODEL_PATH, 'superelements', 'cqrsee101b2.bdf')
        log = SimpleLogger(level='error', encoding='utf-8')
        (fem1, unused_fem2, diff_cards) = self.run_bdf(
            '', bdf_filename, xref=True, run_extract_bodies=False,
            save_file_structure=True, log=log, validate_case_control=False,
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
        log = SimpleLogger(level='error', encoding='utf-8')
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
        #log = SimpleLogger(level='error', encoding='utf-8')

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
        #self.run_bdf(folder, bdf_filename, xref=True, log=log) # PBEAML is not supported

    def test_bdf_thermal_01(self):
        """checks time_thermal_elements.bdf"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'elements', 'time_thermal_elements.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
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
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'transfer_function', 'actuator_tf_modeling.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
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
        log = SimpleLogger(level='warning', encoding='utf-8')
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'aerobeam.bdf')
        fem1, fem2, diff_cards = self.run_bdf('', bdf_filename, log=log)
        #fa2j_1, _rows, _cols = fem1.dmi['FA2J'].get_matrix()
        #fa2j_2, _rows, _cols = fem2.dmi['FA2J'].get_matrix()
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
        #log = SimpleLogger(level='error', encoding='utf-8')
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
    #mass1, cg1, I1 = mass_properties_no_xref(fem1, reference_point=reference_point, sym_axis=sym_axis)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
