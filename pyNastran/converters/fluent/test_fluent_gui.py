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
