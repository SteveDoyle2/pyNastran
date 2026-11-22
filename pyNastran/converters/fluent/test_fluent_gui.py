import os
import unittest
from pathlib import Path
import numpy as np
from cpylog import get_logger

import pyNastran
from pyNastran.gui.testing_methods import FakeGUIMethods
#from pyNastran.bdf.bdf import BDF
from pyNastran.converters.fluent.fluent import Fluent
from pyNastran.converters.fluent.fluent_io import FluentIO
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent
from pyNastran.converters.fluent.ugrid_to_fluent import ugrid_to_fluent_filename
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.qt_files.mark_actions import label_values


PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH /  '..' / 'models'
BWB_PATH = MODEL_PATH / 'bwb'
TEST_PATH = PKG_PATH / 'converters' / 'fluent'
UGRID_PATH = PKG_PATH / 'converters' / 'aflr' / 'ugrid' / 'models'


class FluentGui(FakeGUIMethods):
    def __init__(self):
        FakeGUIMethods.__init__(self)
        self.model = FluentIO(self)
        self.build_fmts(['fluent'], stop_on_failure=True)


class TestFluentGui(unittest.TestCase):

    def test_fluent_geometry_01(self):
        """tests the bwb model"""
        log = get_logger(level='warning', encoding='utf-8')
        #geometry_filename = MODEL_PATH / 'threePlugs.a.tri'

        nastran_filename = BWB_PATH / 'bwb_saero.bdf'
        #vrt_filename2 = BWB_PATH / 'bwb_saero2.vrt'
        vrt_filename = BWB_PATH / 'bwb_saero.vrt'
        #cel_filename = BWB_PATH / 'bwb_saero.cel'
        #daten_filename = BWB_PATH / 'bwb_saero.daten'
        tecplot_filename = BWB_PATH / 'bwb_saero.plt'
        nastran_to_fluent(nastran_filename, vrt_filename, log=log)

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.on_load_geometry(
            vrt_filename, geometry_format='fluent',
            stop_on_failure=True)

    def test_fluent_gui_ugrid3d_gui_box(self):
        """simple UGRID3D box model"""
        ugrid_filename = UGRID_PATH / 'box.b8.ugrid'
        fluent_filename = UGRID_PATH / 'box.vrt'
        h5_filename = UGRID_PATH / 'box.h5'
        if h5_filename.exists():
            os.remove(h5_filename)
        fluent_model = ugrid_to_fluent_filename(ugrid_filename, fluent_filename)

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.on_load_geometry(
            fluent_filename, geometry_format='fluent', stop_on_failure=True)

    def test_fluent_gui_flat_plate_h5(self):
        """
        Flat plate h5 with extra result columns, loaded like ``pyNastranGUI -f fluent plate.h5``:

        - a 3x2 quad plate (6 elements) with results Cp (float), Float1 (float),
          Int1 (whole numbers, e.g. a per-element id/flag) and String1 (a string
          column: integer codes + Fluent.result_labels)
        - every column is a GUI result, in file order, after the standard ones
        - Cp/Float1 round trip exactly as floats ('%.2f'); Int1 is shown as an
          int result ('%i') with the same values and no labels
        - String1: the labels survive the h5 round trip; the legend annotates
          each code with its label and hides the numeric ticks; switching back
          to a numeric result removes the annotations; the probe shows the label
        - a column with labels must use whole-number codes, and every code
          needs a label (GuiResult raises otherwise)
        """
        log = get_logger(level='warning', encoding='utf-8')
        nx, ny = 3, 2
        x, y = np.meshgrid(np.linspace(0., 1., nx + 1), np.linspace(-0.5, 0.5, ny + 1),
                           indexing='ij')
        model = Fluent(log=log)
        model.node_id = np.arange(1, x.size + 1)
        model.xyz = np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])
        nid = model.node_id.reshape(nx + 1, ny + 1)
        quads = []
        for i in range(nx):
            for j in range(ny):
                quads.append([nid[i, j], nid[i + 1, j], nid[i + 1, j + 1], nid[i, j + 1]])
        nquad = len(quads)
        eids = np.arange(1, nquad + 1)
        model.quads = np.column_stack([eids, np.ones(nquad, dtype='int32'), quads])
        model.tris = np.zeros((0, 5), dtype='int32')
        model.element_ids = eids
        model.result_element_id = eids
        model.titles = ['ElementID', 'Cp', 'Float1', 'Int1', 'String1']
        cp = np.linspace(0.05, 0.3, nquad)
        float1 = np.linspace(-5.0, 5.0, nquad)
        int1 = np.array([4, 4, 4, 23, 23, 20])
        string1 = np.array([1, 1, 2, 2, 3, 1])
        labels = {1: 'alpha', 2: 'beta', 3: 'gamma'}
        model.results = np.column_stack([cp, float1, int1, string1])
        model.result_labels = {'String1': labels}

        h5_filename = TEST_PATH / 'tmp_flat_plate.h5'
        model.write_h5(h5_filename)
        test = FluentGui()
        test.log = log
        test.on_load_geometry(h5_filename, geometry_format='fluent', stop_on_failure=True)

        cases = {name: obj for obj, (_itime, name) in test.result_cases.values()}
        names = [name for _obj, (_itime, name) in test.result_cases.values()]
        self.assertEqual(names[-4:], ['Cp', 'Float1', 'Int1', 'String1'])
        np.testing.assert_array_equal(cases['Cp'].scalar, cp)
        np.testing.assert_array_equal(cases['Float1'].scalar, float1)
        self.assertEqual(cases['Cp'].data_format, '%.2f')
        self.assertEqual(cases['Int1'].scalar.dtype.kind, 'i')
        self.assertEqual(cases['Int1'].data_format, '%i')
        np.testing.assert_array_equal(cases['Int1'].scalar, int1)
        self.assertIsNone(cases['Int1'].get_labels(0, 'Int1'))

        string_res = cases['String1']
        self.assertEqual(string_res.get_labels(0, 'String1'), labels)
        np.testing.assert_array_equal(string_res.scalar, string1)

        test.setup_fake_text_actors()
        icase = {name: icase for icase, (_obj, (_i, name)) in test.result_cases.items()}
        lut = test.scalar_bar.color_function
        test._set_case('dummy', icase['String1'], explicit=False, cycle=False,
                       skip_click_check=False, min_value=None, max_value=None,
                       is_legend_shown=None, show_msg=False)
        annotations = {lut.GetAnnotatedValue(j).ToInt(): lut.GetAnnotation(j)
                       for j in range(lut.GetNumberOfAnnotatedValues())}
        self.assertEqual(annotations, labels)
        self.assertFalse(test.scalar_bar.scalar_bar.GetDrawTickLabels())
        self.assertEqual(label_values(string_res, 0, 'String1', string1[4]), 'gamma')
        self.assertEqual(list(label_values(string_res, 0, 'String1', string1[:3])),
                         ['alpha', 'alpha', 'beta'])
        self.assertEqual(label_values(cases['Int1'], 0, 'Int1', 23), 23)

        test._set_case('dummy', icase['Cp'], explicit=False, cycle=False,
                       skip_click_check=False, min_value=None, max_value=None,
                       is_legend_shown=None, show_msg=False)
        self.assertEqual(lut.GetNumberOfAnnotatedValues(), 0)
        self.assertTrue(test.scalar_bar.scalar_bar.GetDrawTickLabels())

        with self.assertRaises(AssertionError):   # float codes
            GuiResult(0, 'bad', 'bad', 'centroid', np.array([1.5, 2.0]), labels=labels)
        with self.assertRaises(AssertionError):   # code 4 has no label
            GuiResult(0, 'bad', 'bad', 'centroid', np.array([1, 4]), labels=labels)
        os.remove(h5_filename)

    def test_fluent_gui_missing_nodes(self):
        model = Fluent()
        model.node_id = np.array([1, 2, 3, 4])
        model.xyz = np.array([
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
        ])
        model.tris = np.array([
            [1, 10, 1, 2, 3],
        ])
        model.quads = np.array([
            [2, 12, 1, 2, 3, 4],
        ])
        model.result_element_id = np.array([1, 2])
        model.element_ids = np.array([1, 2])
        model.titles = ['ShellID', 'Pi']
        model.results = np.ones((len(model.result_element_id), 1)) * 3.14

        log = get_logger(level='warning', encoding='utf-8')
        test = FluentGui()
        test.log = log
        test.model.load_fluent_geometry(model)

    def test_fluent_gui_cp_force_moment_unitless(self):
        """
        Cp integration with the old defaults (sref=cref=bref=1, xyz_ref=0)
        and no unit transform

        1x1 quad at z=0 (+z normal); Cp=1 -> CFz=1
        centroid=(0.5, 0.5, 0) -> CM = r x F = (0.5, -0.5, 0)
        """
        model = _two_element_cp_model()
        test = FluentGui()
        test.log = get_logger(level='warning', encoding='utf-8')
        test.model.load_fluent_geometry(model)

        cf, cm, region_dict = test.model.force_moment_coefficients
        # quad: A=1,   n=+z, Cp=1.0,  centroid=(0.5, 0.5, 0)
        # tri:  A=0.5, n=+z, Cp=-0.5, centroid=(4/3, 1/3, 0)
        cf_quad = np.array([0., 0., 1.0])
        cf_tri = np.array([0., 0., -0.25])
        cm_quad = np.cross([0.5, 0.5, 0.], cf_quad)
        cm_tri = np.cross([4/3, 1/3, 0.], cf_tri)
        assert np.allclose(cf, cf_quad + cf_tri), cf
        assert np.allclose(cm, cm_quad + cm_tri), cm
        assert np.allclose(region_dict[10][0], cf_quad)
        assert np.allclose(region_dict[10][1], cm_quad)
        assert np.allclose(region_dict[20][0], cf_tri)
        assert np.allclose(region_dict[20][1], cm_tri)

    def test_fluent_gui_cp_force_moment_units(self):
        """
        Model is in meters; sref/cref/bref/xyz_ref are in inches.

        Checks:
         - area is converted m^2 -> in^2 before dividing by sref (in^2)
         - centroid is converted m -> in before subtracting xyz_ref (in)
         - lref = [bref, cref, bref] for [roll, pitch, yaw]
         - the Pa->psi pressure conversion hits 'Pressure', not the Cp column
        """
        model = _two_element_cp_model()
        test = FluentGui()
        test.log = get_logger(level='warning', encoding='utf-8')
        other_settings = test.settings.other_settings
        other_settings.units_model_in = ('m', 'N', 's', 'Pa')
        other_settings.units_length = 'in'
        other_settings.units_pressure = 'psi'
        other_settings.sref = 1.0         # in^2
        other_settings.cref = 2.0         # in
        other_settings.bref = 4.0         # in
        other_settings.xyz_ref = np.array([10., 20., 0.])  # in
        test.model.load_fluent_geometry(model)

        m_to_in = 1. / 0.0254
        sref = 1.0

        lref = np.array([4.0, 2.0, 4.0])
        xyz_ref = np.array([10., 20., 0.])
        # quad: A=1 m^2,   Cp=1.0,  centroid=(0.5, 0.5, 0) m
        # tri:  A=0.5 m^2, Cp=-0.5, centroid=(4/3, 1/3, 0) m
        cf_quad = np.array([0., 0., 1.0 * 1.0 * m_to_in**2 / sref])
        cf_tri = np.array([0., 0., -0.5 * 0.5 * m_to_in**2 / sref])
        r_quad = np.array([0.5, 0.5, 0.]) * m_to_in - xyz_ref
        r_tri = np.array([4/3, 1/3, 0.]) * m_to_in - xyz_ref
        cm_quad = np.cross(r_quad, cf_quad) / lref
        cm_tri = np.cross(r_tri, cf_tri) / lref

        cf, cm, region_dict = test.model.force_moment_coefficients
        assert np.allclose(cf[2], 1550.0031 * 0.75), cf  # 1 m^2 = 1550.0031 in^2
        assert np.allclose(cf, cf_quad + cf_tri), cf
        assert np.allclose(cm, cm_quad + cm_tri), cm
        assert np.allclose(region_dict[10][0], cf_quad)
        assert np.allclose(region_dict[10][1], cm_quad)
        assert np.allclose(region_dict[20][0], cf_tri)
        assert np.allclose(region_dict[20][1], cm_tri)

        # the GUI results: Pressure converted to psi; Cp untouched
        cases = test.result_cases
        titles = {case[1][1]: case[0] for case in cases.values()}
        pressure = titles['Pressure'].scalar
        cp = titles['Pressure Coefficient'].scalar
        assert np.allclose(pressure, np.array([101325., 50000.]) / 6894.757), pressure
        assert np.allclose(cp, [1.0, -0.5]), cp

    def test_fluent_cp_force_moment_function(self):
        """same as the GUI case, but calls the function directly; m -> ft"""
        from pyNastran.converters.fluent.fluent_io import compute_cp_force_moment_coefficients
        cp = np.array([2.0])
        area = np.array([0.3048**2])          # 1 ft^2 in m^2
        centroid = np.array([[0.3048, 0., 0.]])  # x=1 ft
        normal = np.array([[0., 0., 1.]])
        region = np.array([1])
        cf, cm, region_dict = compute_cp_force_moment_coefficients(
            cp, area, centroid, normal, region,
            sref=0.5, cref=0.25, bref=1.0, xyz_ref=[0.5, 0., 0.],
            units_length_in='m', units_length_out='ft')
        # CFz = Cp*A/sref = 2*1/0.5 = 4
        # CMy = -(x-xref)*CFz/cref = -(0.5)*4/0.25 = -8
        assert np.allclose(cf, [0., 0., 4.]), cf
        assert np.allclose(cm, [0., -8., 0.]), cm


def _two_element_cp_model() -> Fluent:
    """
    Flat plate at z=0 with a +z normal

    quad (region 10): (0,0)-(1,0)-(1,1)-(0,1); area=1.0
    tri  (region 20): (1,0)-(2,0)-(1,1);       area=0.5
    """
    model = Fluent()
    model.node_id = np.array([1, 2, 3, 4, 5])
    model.xyz = np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [2., 0., 0.],
    ])
    # [eid, region, n1, n2, ...]
    model.quads = np.array([
        [1, 10, 1, 2, 3, 4],
    ])
    model.tris = np.array([
        [2, 20, 2, 5, 3],
    ])
    model.result_element_id = np.array([1, 2])
    model.element_ids = np.array([1, 2])
    # Cp is intentionally not the last column
    model.titles = ['ShellID', 'Pressure', 'Pressure Coefficient', 'Pi']
    model.results = np.array([
        [101325., 1.0, 3.14],
        [50000., -0.5, 3.14],
    ])
    return model


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
