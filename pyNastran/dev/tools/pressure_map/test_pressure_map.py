from pathlib import Path
import tempfile
import unittest
import numpy as np

from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.aero.export_caero_mesh import export_caero_mesh
from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.fluent.fluent import Fluent
from pyNastran.converters.fluent.nastran_to_fluent import nastran_to_fluent

from pyNastran.dev.tools.pressure_map.pressure_map import (
    pressure_map, pressure_filename_to_fa2j, pressure_filename_to_wkk_diag,
    map_panel_force_moment_centroid, pressure_map_to_panel_model,
    pressure_map_to_structure_model, pressure_map_to_lcp_model,
    _check_lcp_group_location)
from pyNastran.dev.tools.pressure_map.setup_aero import (
    get_aero_model, get_aero_pressure_centroid, log_aero_force_moment)

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = (PKG_PATH / '..' / 'models').resolve()
DIRNAME = Path(__file__).parent


class TestPressureMap(unittest.TestCase):
    def test_pressure_map_cart3d(self):
        aero_format = 'cart3d'
        cart3d_filename = PKG_PATH / 'converters' / 'cart3d' / 'models' / 'threePlugs.bin.tri'
        bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.bdf'
        caero_bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.caero.bdf'

        skip_cards = ['CBAR']
        log = SimpleLogger(level='warning')
        bdf_model = read_bdf(bdf_filename, skip_cards=skip_cards, log=log)
        # if not caero_bdf_filename.exists():  # pragma: no cover
        export_caero_mesh(
            bdf_filename, caero_bdf_filename,
            is_aerobox_model=True,
            write_panel_xyz=False,
            write_header=False,
            write_end_data=False)

        aero_model, variables = get_aero_model(
            cart3d_filename, aero_format,
                   aero_xyz_scale=1.0,
                   xyz_units_out='in',
                   stop_on_failure=True, log=log)

        neids = len(aero_model.elements)
        eids = np.arange(neids)
        cp = np.sin(eids/neids)
        aero_model.regions = np.ones(neids)
        aero_model.loads['Cp'] = cp

        pressure_map(
            aero_model,
            bdf_model,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='pressure',
            force_to_pressure_tol=0.02,
            method='full_model',
            xyz_units_out='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'cart3d_pressure_fullmodel_1.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None,
            log=log, is_obj=False)

        pressure_map(
            aero_model,
            bdf_model,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force',
            method='full_model',
            xyz_units_out='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'cart3d_force_fullmodel_2.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None,
            log=log)

        with self.assertRaises(RuntimeError):
            pressure_map(
                aero_model,
                bdf_model,
                # eids_structure=np.array([]),
                # eid_csv_filename='',
                eid_load_id=-1,
                aero_format=aero_format,
                map_type='force_moment',
                method='full_model',
                xyz_units_out='in',
                pressure_units='psi',
                pressure_sid=1,
                force_sid=2,
                moment_sid=3,
                idtype='int32', fdtype='float64',
                pressure_filename=DIRNAME/'cart3d_forcemoment_fullmodel_3.bdf',
                aero_xyz_scale=1.0, qinf=1.0,
                sref=1.0, cref=1.0, bref=1.0,
                reference_point=None,
                regions_to_include=None,
                regions_to_remove=None,
                log=log)

        # log = SimpleLogger(level='debug')
        pressure_filename = DIRNAME / 'cart3d_forcemoment_panelmodel_4.bdf'
        main_pressure_filename = DIRNAME / 'main_cart3d_forcemoment_panelmodel_4.bdf'
        caero_model = read_bdf(caero_bdf_filename,
                               punch=True, xref=True)
        pressure_map(
            aero_model, #cart3d_filename,
            caero_model, #caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force_moment',
            method='panel_model',
            xyz_units_out='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=pressure_filename,
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None,
            log=log)

        fa2j_filename = DIRNAME / 'main_cart3d_fa2j_5.bdf'
        log.debug('working on fa2j writer')
        pressure_filename_to_fa2j(
            main_pressure_filename, fa2j_filename, sid=1, log=log)

        wkk_filename = DIRNAME / 'cart3d_wkk_6.bdf'
        pressure_filename1 = main_pressure_filename
        pressure_filename2 = main_pressure_filename
        log.debug('working on wkk writer')
        pressure_filename_to_wkk_diag(
            pressure_filename1,
            pressure_filename2,
            np.zeros(2),
            np.ones(2),
            wkk_filename,
            force_sid1=2, moment_sid1=3,
            force_sid2=1, moment_sid2=3,
            log=log)

    def test_pressure_map_fluent(self):
        aero_format = 'fluent'
        # map_type = 'pressure'
        bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.bdf'
        vrt_filename = MODEL_DIR / 'bwb' / 'bwb_saero.vrt'
        caero_bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.caero.bdf'
        # if not caero_bdf_filename.exists():  # pragma: no cover
        export_caero_mesh(
            bdf_filename, caero_bdf_filename,
            is_aerobox_model=True,
            write_panel_xyz=False,
            write_header=True,
            write_end_data=True)

        log = SimpleLogger(level='warning')
        # log = SimpleLogger(level='info')
        # if not vrt_filename.exists():  # pragma: no cover
        nastran_to_fluent(bdf_filename, vrt_filename, log=log)

        aero_model, variables = get_aero_model(
            vrt_filename, aero_format,
                   aero_xyz_scale=1.0,
                   xyz_units_out='in',
                   stop_on_failure=True, log=log)
        aero_model.titles = ['ElementID', 'Pressure Coefficient']
        # get_aero_pressure_centroid(
        #     aero_model, aero_format,
        #     map_type, variable='Cp',
        #     regions_to_include=None,
        #     regions_to_remove=None)
        log.info('running 1st pressure map')
        # log.level = 'debug'
        pressure_map(
            aero_model, #cart3d_filename,
            caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force_moment',
            method='panel_model',
            xyz_units_out='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'fluent_forcemoment_panelmodel_1.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None,
            log=log)

        log.info('running 2nd pressure map')
        pressure_map(
            aero_model, #cart3d_filename,
            caero_bdf_filename,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force',
            method='full_model',
            xyz_units_out='in',
            pressure_units='psi',
            pressure_sid=1,
            force_sid=2,
            moment_sid=3,
            idtype='int32', fdtype='float64',
            pressure_filename=DIRNAME/'fluent_force_fullmodel_2.bdf',
            aero_xyz_scale=1.0, qinf=1.0,
            sref=1.0, cref=1.0, bref=1.0,
            reference_point=None,
            regions_to_include=None,
            regions_to_remove=None,
            log=log)

        log.info('running 3rd pressure map (force -> pressure)')
        pressure_filename = DIRNAME / 'fluent_force_pressure_fullmodel_3.bdf'
        pressure_model = pressure_map(
            aero_model,
            caero_bdf_filename,
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='force',
            method='full_model',
            xyz_units_out='in',
            pressure_units='psi',
            force_sid=2,
            pressure_filename=pressure_filename,
            aero_xyz_scale=1.0, qinf=1.0,
            force_to_pressure_tol=0.05,
            log=log)
        load_types = {load.type for load in pressure_model.loads[2]}
        self.assertTrue(load_types <= {'FORCE', 'PLOAD2'}, load_types)
        if 'PLOAD2' in load_types:
            self.assertTrue(pressure_filename.with_suffix('.pressure.csv').exists())
        if 'FORCE' in load_types:
            self.assertTrue(pressure_filename.with_suffix('.force.csv').exists())

    def test_fluent_lcp(self):
        bwb_dirname = MODEL_DIR / 'bwb'
        aero_format = 'fluent'
        log = SimpleLogger(level='warning')
        # log = SimpleLogger(level='debug')
        #bdf_filename = MODEL_DIR / 'bwb' / 'bwb_saero.bdf'

        bdf_filename = bwb_dirname / 'bwb_saero.bdf'
        vrt_filename = bwb_dirname / 'bwb_saero.vrt'
        nastran_to_fluent(bdf_filename, vrt_filename, log=log)

        reference_point = np.array([807, 0, 0], dtype='float64')
        aero_model, variables = get_aero_model(
            vrt_filename, aero_format,
                   aero_xyz_scale=1.0,
                   xyz_units_in='in',
                   xyz_units_out='in',
                   reference_point=reference_point,
                   stop_on_failure=True, log=log)
        aero_model.titles = aero_model.titles.tolist()

        #from pyNastran.utils import object_stats
        #print(object_stats(aero_model))
        neid = len(aero_model.results)
        eids = np.arange(neid)
        cp = np.sin(eids/neid)
        aero_model.results[:, -1] = cp
        aero_model.titles[-1] = 'Pressure Coefficient'

        #model = read_bdf(bdf_filename)
        #nid0 = max(model.nodes) + 1
        model = BDF(log=log)
        model.bdf_filename = bwb_dirname / 'bwb_lcp_nodes.bdf'
        nid0 = 1
        x = np.linspace(0., 1614., num=21)
        for inid, xi in enumerate(x):
            nid = nid0 + inid
            model.add_grid(nid, [xi, 0., 0.])

        log.info('running LCP / force_moment map')
        pressure_filename = bwb_dirname / 'fluent_LCP.bdf'
        pressure_model = pressure_map(
            aero_model,
            model,
            eid_load_id=-1,
            reference_point=reference_point,
            aero_format='fluent',
            map_type='force_moment',
            method='LCP',
            xyz_units_in='in',
            xyz_units_out='in',
            pressure_units='psi',
            force_sid=2,
            moment_sid=3,
            pressure_filename=pressure_filename,
            aero_xyz_scale=1.0, qinf=1.0,
            log=log)

def _empty_quad_arrays():
    """Return empty aero quad arrays (area, Cp, normal, centroid, force/q)."""
    return (
        np.zeros(0, dtype='float64'),       # area
        np.zeros(0, dtype='float64'),       # Cp
        np.zeros((0, 3), dtype='float64'),  # normal
        np.zeros((0, 3), dtype='float64'),  # centroid
        np.zeros((0, 3), dtype='float64'),  # force/q
    )


def _call_panel_force_moment(
        panel_dim,
        structure_eids, structure_xyz, structure_area, structure_z_sign,
        aero_tri_centroid, aero_tri_cp, aero_tri_area, aero_tri_normal,
        iaero_tri,
        qinf=1.0,
        reference_point=None,
        sref=1.0, bref=1.0, cref=1.0,
    ):
    """Thin wrapper around map_panel_force_moment_centroid for testing.

    Builds the derived arrays (force/q, empty quads, BDF output model)
    so each test only specifies the physically meaningful inputs.

    Returns
    -------
    struct_force : (3,) float ndarray
        panel-normal structural force  [0,0,ΣFz] or [0,ΣFy,0]
    struct_moment : (3,) float ndarray
        structural moment about reference_point (normal-force only)
    """
    log = SimpleLogger(level='warning')
    bdf_out = BDF(log=log)

    # F/q per aero tri = Cp * A * n
    aero_tri_force_per_q = (
        (aero_tri_cp * aero_tri_area)[:, np.newaxis] * aero_tri_normal
    )
    quad = _empty_quad_arrays()
    iaero_quad = np.array([], dtype='int32')

    struct_force, struct_moment = map_panel_force_moment_centroid(
        bdf_out, panel_dim,
        aero_tri_area, aero_tri_cp, aero_tri_normal,
        aero_tri_centroid, aero_tri_force_per_q,
        *quad,
        structure_eids, structure_xyz, structure_area, structure_z_sign,
        iaero_tri, iaero_quad,
        flip_Cp=False, qinf=qinf,
        reference_point=reference_point,
        sref=sref, bref=bref, cref=cref,
    )
    return struct_force, struct_moment


class TestPanelModelBalance(unittest.TestCase):
    """Unit tests for the force/moment balance in map_panel_force_moment_centroid.

    Each test constructs synthetic aero and structural data with known
    analytical answers, then verifies the returned structural force and
    moment vectors.
    """

    # ------------------------------------------------------------------
    # Test 1: force conservation on z-panels
    # ------------------------------------------------------------------
    def test_z_panel_force_conservation(self):
        """Total structural Fz must equal total aero Fz; Fx=Fy=0."""
        # Two z-panels at y = ±25
        structure_eids = np.array([1, 2], dtype='int32')
        structure_xyz = np.array([[5., -25., 0.],
                                  [5.,  25., 0.]])
        structure_area = np.array([500., 500.])
        structure_z_sign = np.array([1., 1.])

        # 4 aero tris, all with normal pointing +z
        aero_tri_centroid = np.array([
            [4., -30., 0.5],
            [6., -20., 0.5],
            [4.,  20., 0.5],
            [6.,  30., 0.5],
        ])
        aero_tri_cp = np.array([1.0, 1.0, 2.0, 2.0])
        aero_tri_area = np.array([100., 100., 100., 100.])
        aero_tri_normal = np.tile([0., 0., 1.], (4, 1))

        # mapping: first two tris → panel 0, last two → panel 1
        iaero_tri = np.array([0, 0, 1, 1])

        sf, sm = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri)

        # F/q = Cp*A*n, qinf=1 → F = F/q
        # Panel 0: Fz = 1*100 + 1*100 = 200
        # Panel 1: Fz = 2*100 + 2*100 = 400
        # Total Fz = 600
        expected_fz = 600.0
        self.assertAlmostEqual(sf[0], 0.0, places=10, msg='Fx must be 0')
        self.assertAlmostEqual(sf[1], 0.0, places=10, msg='Fy must be 0')
        self.assertAlmostEqual(sf[2], expected_fz, places=10,
                               msg='Fz not conserved')

    # ------------------------------------------------------------------
    # Test 2: force conservation on y-panels
    # ------------------------------------------------------------------
    def test_y_panel_force_conservation(self):
        """Total structural Fy must equal total aero Fy; Fx=Fz=0."""
        structure_eids = np.array([10, 20], dtype='int32')
        structure_xyz = np.array([[5., 0., -25.],
                                  [5., 0.,  25.]])
        structure_area = np.array([500., 500.])
        structure_z_sign = np.array([1., 1.])  # sign(normal_y)

        # aero tris with normal in +y
        aero_tri_centroid = np.array([
            [5., 0.5, -30.],
            [5., 0.5, -20.],
            [5., 0.5,  20.],
            [5., 0.5,  30.],
        ])
        aero_tri_cp = np.array([1.5, 1.5, 0.5, 0.5])
        aero_tri_area = np.array([80., 80., 80., 80.])
        aero_tri_normal = np.tile([0., 1., 0.], (4, 1))
        iaero_tri = np.array([0, 0, 1, 1])

        sf, sm = _call_panel_force_moment(
            'y', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri)

        # Panel 0: Fy = 1.5*80 + 1.5*80 = 240
        # Panel 1: Fy = 0.5*80 + 0.5*80 = 80
        expected_fy = 320.0
        self.assertAlmostEqual(sf[0], 0.0, places=10, msg='Fx must be 0')
        self.assertAlmostEqual(sf[1], expected_fy, places=10,
                               msg='Fy not conserved')
        self.assertAlmostEqual(sf[2], 0.0, places=10, msg='Fz must be 0')

    # ------------------------------------------------------------------
    # Test 3: z-panel moment balance (Mx, My, Mz)
    # ------------------------------------------------------------------
    def test_z_panel_moment_balance(self):
        """Structural moment about origin must match aero moment for
        the Fz-only terms: Mx = Σ(y·Fz), My = -Σ(x·Fz), Mz = 0."""
        structure_eids = np.array([1, 2], dtype='int32')
        structure_xyz = np.array([[5., -25., 0.],
                                  [5.,  25., 0.]])
        structure_area = np.array([500., 500.])
        structure_z_sign = np.array([1., 1.])

        aero_tri_centroid = np.array([
            [4., -30., 0.5],
            [6., -20., 0.5],
            [4.,  20., 0.5],
            [6.,  30., 0.5],
        ])
        aero_tri_cp = np.array([1.0, 1.0, 2.0, 2.0])
        aero_tri_area = np.array([100., 100., 100., 100.])
        aero_tri_normal = np.tile([0., 0., 1.], (4, 1))
        iaero_tri = np.array([0, 0, 1, 1])

        ref = np.zeros(3)
        sf, sm = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            reference_point=ref)

        # Force per aero tri (qinf=1, F = Cp*A*n):
        #   tri 0: [0, 0, 100] at (4, -30, 0.5)
        #   tri 1: [0, 0, 100] at (6, -20, 0.5)
        #   tri 2: [0, 0, 200] at (4,  20, 0.5)
        #   tri 3: [0, 0, 200] at (6,  30, 0.5)
        #
        # For z-normal aero (Fx=Fy=0):
        #   Mx = Σ(y·Fz) = -30*100 + -20*100 + 20*200 + 30*200 = 5000
        #   My = -Σ(x·Fz) = -(4*100 + 6*100 + 4*200 + 6*200) = -3000
        #   Mz = 0  (no Fy, no Fx)
        expected_mx = 5000.0
        expected_my = -3000.0
        expected_mz = 0.0
        self.assertAlmostEqual(sm[0], expected_mx, places=6,
                               msg=f'Mx wrong: {sm[0]} != {expected_mx}')
        self.assertAlmostEqual(sm[1], expected_my, places=6,
                               msg=f'My wrong: {sm[1]} != {expected_my}')
        self.assertAlmostEqual(sm[2], expected_mz, places=6,
                               msg=f'Mz wrong: {sm[2]} != {expected_mz}')

    # ------------------------------------------------------------------
    # Test 4: y-panel moment balance
    # ------------------------------------------------------------------
    def test_y_panel_moment_balance(self):
        """Structural moment about origin must match aero moment for
        the Fy-only terms: Mx = -Σ(z·Fy), My = 0, Mz = Σ(x·Fy)."""
        structure_eids = np.array([10, 20], dtype='int32')
        structure_xyz = np.array([[5., 0., -25.],
                                  [5., 0.,  25.]])
        structure_area = np.array([500., 500.])
        structure_z_sign = np.array([1., 1.])

        aero_tri_centroid = np.array([
            [4., 0.5, -30.],
            [6., 0.5, -20.],
            [4., 0.5,  20.],
            [6., 0.5,  30.],
        ])
        aero_tri_cp = np.array([1.5, 1.5, 0.5, 0.5])
        aero_tri_area = np.array([80., 80., 80., 80.])
        aero_tri_normal = np.tile([0., 1., 0.], (4, 1))
        iaero_tri = np.array([0, 0, 1, 1])

        ref = np.zeros(3)
        sf, sm = _call_panel_force_moment(
            'y', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            reference_point=ref)

        # Force per tri (Fy = Cp*A):
        #   tri 0: Fy=120 at z=-30   tri 1: Fy=120 at z=-20
        #   tri 2: Fy=40  at z=20    tri 3: Fy=40  at z=30
        #
        # Mx = -Σ(z·Fy) = -(-30*120 + -20*120 + 20*40 + 30*40) = -(−6000) + 0 ... let me compute
        #   = -((-30)*120 + (-20)*120 + 20*40 + 30*40)
        #   = -(−3600 − 2400 + 800 + 1200) = -(−4000) = 4000
        # My = 0  (no Fz)
        # Mz = Σ(x·Fy) = 4*120 + 6*120 + 4*40 + 6*40
        #     = 480 + 720 + 160 + 240 = 1600
        expected_mx = 4000.0
        expected_my = 0.0
        expected_mz = 1600.0
        self.assertAlmostEqual(sm[0], expected_mx, places=6,
                               msg=f'Mx wrong: {sm[0]} != {expected_mx}')
        self.assertAlmostEqual(sm[1], expected_my, places=6,
                               msg=f'My wrong: {sm[1]} != {expected_my}')
        self.assertAlmostEqual(sm[2], expected_mz, places=6,
                               msg=f'Mz wrong: {sm[2]} != {expected_mz}')

    # ------------------------------------------------------------------
    # Test 5: z + y combined — no double-counting
    # ------------------------------------------------------------------
    def test_z_plus_y_no_double_counting(self):
        """Summing z-call and y-call results must give correct Fy, Fz, and
        moments without double-counting.

        Uses aero elements with BOTH Fy and Fz components (tilted normal)
        so that each call sees the full force vector but only returns its
        normal component.
        """
        # -- Structural panels --
        # One z-panel at y=50 and one y-panel at z=30
        z_eids = np.array([1], dtype='int32')
        z_xyz = np.array([[10., 50., 0.]])
        z_area = np.array([500.])
        z_sign = np.array([1.])

        y_eids = np.array([2], dtype='int32')
        y_xyz = np.array([[10., 0., 30.]])
        y_area = np.array([500.])
        y_sign = np.array([1.])

        # -- Aero elements --
        # 2 tris with tilted normals (both Fy and Fz components)
        aero_centroid = np.array([
            [10., 40., 5.],
            [10., 60., 5.],
        ])
        # normal = [0, sin(30°), cos(30°)] ≈ [0, 0.5, 0.866]
        ny, nz = 0.5, np.sqrt(3)/2
        aero_normal = np.array([[0., ny, nz],
                                [0., ny, nz]])
        aero_cp = np.array([1.0, 1.0])
        aero_area = np.array([200., 200.])
        # F/q per tri = Cp*A*n = 1*200*[0, 0.5, 0.866] = [0, 100, 173.2]
        fz_per = aero_area[0] * nz
        fy_per = aero_area[0] * ny

        # Both tris map to the single panel in each call
        iaero_1 = np.array([0, 0])

        ref = np.zeros(3)
        sf_z, sm_z = _call_panel_force_moment(
            'z', z_eids, z_xyz, z_area, z_sign,
            aero_centroid, aero_cp, aero_area, aero_normal,
            iaero_1, reference_point=ref)

        sf_y, sm_y = _call_panel_force_moment(
            'y', y_eids, y_xyz, y_area, y_sign,
            aero_centroid, aero_cp, aero_area, aero_normal,
            iaero_1, reference_point=ref)

        total_force = sf_z + sf_y
        total_moment = sm_z + sm_y

        # Expected forces (qinf=1):
        expected_fz = 2 * fz_per  # both tris Fz
        expected_fy = 2 * fy_per  # both tris Fy
        self.assertAlmostEqual(total_force[0], 0.0, places=10)
        self.assertAlmostEqual(total_force[1], expected_fy, places=6)
        self.assertAlmostEqual(total_force[2], expected_fz, places=6)

        # No double-counting: Fy comes only from y-call, Fz only from z-call
        self.assertAlmostEqual(sf_z[1], 0.0, places=10,
                               msg='z-call must not carry Fy')
        self.assertAlmostEqual(sf_y[2], 0.0, places=10,
                               msg='y-call must not carry Fz')

        # Expected moments about origin
        # M = Σ(r_i × F_normal_i), decomposed:
        #   z-call: each tri has F_normal = [0, 0, fz_per]
        #     Mx_z = Σ(y·Fz) = 40*fz_per + 60*fz_per = 100*fz_per
        #     My_z = -Σ(x·Fz) = -(10+10)*fz_per = -20*fz_per
        #     Mz_z = 0
        #   y-call: each tri has F_normal = [0, fy_per, 0]
        #     Mx_y = -Σ(z·Fy) = -(5+5)*fy_per = -10*fy_per
        #     My_y = 0
        #     Mz_y = Σ(x·Fy) = (10+10)*fy_per = 20*fy_per
        expected_mx = 100 * fz_per - 10 * fy_per
        expected_my = -20 * fz_per
        expected_mz = 20 * fy_per
        self.assertAlmostEqual(total_moment[0], expected_mx, places=4,
                               msg=f'Mx wrong: {total_moment[0]} != {expected_mx}')
        self.assertAlmostEqual(total_moment[1], expected_my, places=4,
                               msg=f'My wrong: {total_moment[1]} != {expected_my}')
        self.assertAlmostEqual(total_moment[2], expected_mz, places=4,
                               msg=f'Mz wrong: {total_moment[2]} != {expected_mz}')

    # ------------------------------------------------------------------
    # Test 6: z_sign does not affect the global balance
    # ------------------------------------------------------------------
    def test_z_sign_does_not_affect_balance(self):
        """Flipping z_sign (panel normal direction) must not change the
        global force/moment totals — it only affects the PLOAD2 sign."""
        structure_eids = np.array([1, 2], dtype='int32')
        structure_xyz = np.array([[5., -25., 0.],
                                  [5.,  25., 0.]])
        structure_area = np.array([500., 500.])

        aero_tri_centroid = np.array([
            [5., -25., 1.],
            [5.,  25., 1.],
        ])
        aero_tri_cp = np.array([1.0, 2.0])
        aero_tri_area = np.array([100., 100.])
        aero_tri_normal = np.tile([0., 0., 1.], (2, 1))
        iaero_tri = np.array([0, 1])

        ref = np.zeros(3)
        # z_sign = +1 for both
        sf_pos, sm_pos = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            np.array([1., 1.]),
            aero_tri_centroid, aero_tri_cp, aero_tri_area,
            aero_tri_normal, iaero_tri, reference_point=ref)

        # z_sign = -1 for first panel (as if its normal points -z)
        sf_neg, sm_neg = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            np.array([-1., 1.]),
            aero_tri_centroid, aero_tri_cp, aero_tri_area,
            aero_tri_normal, iaero_tri, reference_point=ref)

        np.testing.assert_allclose(sf_pos, sf_neg, atol=1e-12,
                                   err_msg='z_sign changed global force')
        np.testing.assert_allclose(sm_pos, sm_neg, atol=1e-12,
                                   err_msg='z_sign changed global moment')

    # ------------------------------------------------------------------
    # Test 7: within-panel moment_normal captures Mx from span-wise
    #         Fz distribution
    # ------------------------------------------------------------------
    def test_moment_normal_captures_roll(self):
        """When aero elements within a single panel have different y offsets,
        moment_normal must capture the resulting Mx (roll).

        One wide z-panel centered at y=0.  Aero tris at y=+50 and y=-50
        with unequal Fz → net Mx ≠ 0.  Without moment_normal this would
        be lost (panel center is at y=0).
        """
        structure_eids = np.array([1], dtype='int32')
        structure_xyz = np.array([[5., 0., 0.]])  # panel center at y=0
        structure_area = np.array([1000.])
        structure_z_sign = np.array([1.])

        # Two aero tris, one on each side of the panel center
        aero_tri_centroid = np.array([
            [5., -50., 0.],
            [5.,  50., 0.],
        ])
        # Asymmetric Cp: left side has more lift
        aero_tri_cp = np.array([3.0, 1.0])
        aero_tri_area = np.array([100., 100.])
        aero_tri_normal = np.tile([0., 0., 1.], (2, 1))
        iaero_tri = np.array([0, 0])  # both → same panel

        ref = np.zeros(3)
        sf, sm = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            reference_point=ref)

        # Fz left = 3*100 = 300, Fz right = 1*100 = 100
        # Total Fz = 400
        self.assertAlmostEqual(sf[2], 400.0, places=10)

        # Mx = y_left * Fz_left + y_right * Fz_right
        #    = (-50)*300 + 50*100 = -15000 + 5000 = -10000
        #
        # Without moment_normal, Mx would be:
        #   panel_y * total_Fz = 0 * 400 = 0  (WRONG)
        #
        # With moment_normal, we get the correction:
        #   moment_normal[0] = (-50-0)*300 + (50-0)*100 = -10000
        #   cross(dr_ref=[5,0,0], [0,0,400])[0] = 0
        #   Mx = -10000 + 0 = -10000  (CORRECT)
        expected_mx = -10000.0
        self.assertAlmostEqual(sm[0], expected_mx, places=6,
                               msg=f'Mx wrong: {sm[0]} != {expected_mx}; '
                                   f'moment_normal not capturing roll')

        # My = -Σ(x*Fz) = -(5*300 + 5*100) = -2000
        expected_my = -2000.0
        self.assertAlmostEqual(sm[1], expected_my, places=6)

    # ------------------------------------------------------------------
    # Test 8: non-origin reference point
    # ------------------------------------------------------------------
    def test_nonzero_reference_point(self):
        """Moments about a non-origin reference point must be correct."""
        structure_eids = np.array([1], dtype='int32')
        structure_xyz = np.array([[10., 0., 0.]])
        structure_area = np.array([500.])
        structure_z_sign = np.array([1.])

        aero_tri_centroid = np.array([[10., 0., 0.]])
        aero_tri_cp = np.array([1.0])
        aero_tri_area = np.array([200.])
        aero_tri_normal = np.array([[0., 0., 1.]])
        iaero_tri = np.array([0])

        # ref at [5, 0, 0] → r = centroid - ref = [5, 0, 0]
        ref = np.array([5., 0., 0.])
        sf, sm = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            reference_point=ref)

        # F = [0, 0, 200]
        # r = panel_center - ref = [5, 0, 0]
        # M = r × F = [0*200-0*0, 0*0-5*200, 5*0-0*0] = [0, -1000, 0]
        self.assertAlmostEqual(sf[2], 200.0, places=10)
        self.assertAlmostEqual(sm[0], 0.0, places=10)
        self.assertAlmostEqual(sm[1], -1000.0, places=6)
        self.assertAlmostEqual(sm[2], 0.0, places=10)

    # ------------------------------------------------------------------
    # Test 9: qinf scaling
    # ------------------------------------------------------------------
    def test_qinf_scaling(self):
        """Forces and moments must scale linearly with qinf."""
        structure_eids = np.array([1], dtype='int32')
        structure_xyz = np.array([[5., 20., 0.]])
        structure_area = np.array([500.])
        structure_z_sign = np.array([1.])

        aero_tri_centroid = np.array([[5., 20., 0.]])
        aero_tri_cp = np.array([1.0])
        aero_tri_area = np.array([100.])
        aero_tri_normal = np.array([[0., 0., 1.]])
        iaero_tri = np.array([0])

        ref = np.zeros(3)
        sf1, sm1 = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            qinf=1.0, reference_point=ref)

        sf3, sm3 = _call_panel_force_moment(
            'z', structure_eids, structure_xyz, structure_area,
            structure_z_sign, aero_tri_centroid, aero_tri_cp,
            aero_tri_area, aero_tri_normal, iaero_tri,
            qinf=3.0, reference_point=ref)

        np.testing.assert_allclose(sf3, 3.0 * sf1, atol=1e-12,
                                   err_msg='force does not scale with qinf')
        np.testing.assert_allclose(sm3, 3.0 * sm1, atol=1e-12,
                                   err_msg='moment does not scale with qinf')


#---------------------------------------------------------------------------
# PLOAD2 values written by the panel model
#
# Structure: one 4x2 CQUAD4 in z=0 (area=8), nodes counter-clockwise
# from (0,0), normal +z -> z_sign=+1.
#   leading-edge midpoint p14 = (0, 1, 0), chord = 4
#   quarter-chord moment center = (1, 1, 0)
#   (the old code added chord/4 to x, y AND z -> (1, 2, 1))
#
# Aero: one flat element (n=+z, Cp=1) and one element tilted 45 deg in x
# (n=[-1,0,1]/sqrt(2), Cp=2), both at z>0, so Fx contributes to My.
# Everything maps to the one panel.
#
# PLOAD2 sids written per panel (z-panels):
#   1 pressure = Fz/A        2 force = Fz       3 moment = My about the
#   4 Cp = p/qinf            5 signed Cp = z_sign*Cp      moment center
#---------------------------------------------------------------------------
PANEL_QINF = 2.0


def _one_quad_structure(log: SimpleLogger) -> BDF:
    """single 4x2 CQUAD4 in the z=0 plane; eid=1"""
    model = BDF(log=log)
    xyzs = [(0., 0., 0.), (4., 0., 0.), (4., 2., 0.), (0., 2., 0.)]
    for nid, xyz in enumerate(xyzs, start=1):
        model.add_grid(nid, list(xyz))
    model.add_cquad4(1, 10, [1, 2, 3, 4])
    model.add_pshell(10, mid1=100, t=0.1)
    model.add_mat1(100, 1.0e7, None, 0.3)
    model.cross_reference()
    return model


def _two_tri_cart3d(log: SimpleLogger, xyz_scale: float=1.0) -> Cart3D:
    """
    flat tri:   area=1,         n=+z,            Cp=1, centroid=(7/6, 5/6, 1/2)
    tilted tri: area=sqrt(2)/2, n=[-1,0,1]/sqrt2, Cp=2, centroid=(17/6, 5/6, 5/6)
    """
    model = Cart3D(log=log)
    model.points = xyz_scale * np.array([
        [0.5, 0.5, 0.5],
        [2.5, 0.5, 0.5],
        [0.5, 1.5, 0.5],
        [3.5, 0.5, 1.5],
        [2.5, 1.5, 0.5],
    ])
    model.regions = np.array([1, 1], dtype='int32')
    model.elements = np.array([[0, 1, 2],
                               [1, 3, 4]], dtype='int32')
    model.loads = {'Cp': np.array([1.0, 2.0])}
    return model


def _two_quad_fluent(log: SimpleLogger) -> Fluent:
    """
    no tris
    flat quad:   area=2,       n=+z,             Cp=1, centroid=(1.5, 1, 0.5)
    tilted quad: area=sqrt(2), n=[-1,0,1]/sqrt2, Cp=2, centroid=(3, 1, 1)
    """
    node_id = np.arange(1, 7, dtype='int32')
    xyz = np.array([
        [0.5, 0.5, 0.5],
        [2.5, 0.5, 0.5],
        [2.5, 1.5, 0.5],
        [0.5, 1.5, 0.5],
        [3.5, 0.5, 1.5],
        [3.5, 1.5, 1.5],
    ])
    # [eid, region, n1, n2, n3, n4]
    quads = np.array([[1, 3, 1, 2, 3, 4],
                      [2, 3, 2, 5, 6, 3]], dtype='int32')
    tris = np.zeros((0, 5), dtype='int32')
    titles = np.array(['ElementID', 'Pressure Coefficient'])
    quad_results = np.array([[1.0], [2.0]])
    tri_results = np.zeros((0, 1))
    return Fluent.from_data(
        node_id, xyz, tris, quads,
        np.array([1, 2], dtype='int32'), titles,
        quad_results, tri_results,
        auto_read_write_h5=False, log=log, debug=False)


def _pload2_values(model: BDF) -> dict[int, float]:
    """{sid: pressure} for a model with one PLOAD2 per sid"""
    values = {}
    for sid, loads in model.loads.items():
        assert len(loads) == 1, (sid, loads)
        values[sid] = loads[0].pressure
    return values


class TestPanelModelLoads(unittest.TestCase):
    """End-to-end checks of the PLOAD2 card values from the panel model."""

    def _map(self, aero_model):
        log = aero_model.log
        structure_model = _one_quad_structure(log)
        out_model = pressure_map_to_panel_model(
            aero_model, structure_model,
            np.array([1], dtype='int32'),
            reference_point=None,
            map_type='force_moment',
            qinf=PANEL_QINF)
        return _pload2_values(out_model)

    def test_pload2_values_tri_aero(self):
        """pressure/force/moment/Cp cards for a known 2-tri aero model"""
        log = SimpleLogger(level='warning')
        values = self._map(_two_tri_cart3d(log))

        # F/q = sum(Cp*A*n) = [0,0,1] + [-1,0,1] = [-1, 0, 2]
        # F = qinf*F/q = [-2, 0, 4]
        # My about (1,1,0) = sum(dz*Fx - dx*Fz):
        #   flat:   d=(1/6, -1/6, 1/2),   F=(0,0,2)  -> -1/3
        #   tilted: d=(11/6, -1/6, 5/6),  F=(-2,0,2) -> -16/3
        #   total = -17/3
        #   (with the old (1,2,1) center this was -11/3)
        expected = {
            1: 4.0 / 8.0,                 # pressure = Fz/A
            2: 4.0,                       # Fz
            3: -17.0 / 3.0,               # My about the quarter chord
            4: 4.0 / 8.0 / PANEL_QINF,    # Cp
            5: 4.0 / 8.0 / PANEL_QINF,    # signed Cp (z_sign=+1)
        }
        self.assertEqual(sorted(values), sorted(expected))
        for sid, value in expected.items():
            self.assertAlmostEqual(values[sid], value, places=10, msg=f'sid={sid}')

    def test_pload2_values_quad_only_aero(self):
        """a quad-only aero model must still get the quarter-chord moment center"""
        log = SimpleLogger(level='warning')
        values = self._map(_two_quad_fluent(log))

        # F/q = [0,0,2] + [-2,0,2] = [-2, 0, 4];  F = [-4, 0, 8]
        # My about (1,1,0):
        #   flat:   d=(0.5, 0, 0.5), F=(0,0,4)  -> -2
        #   tilted: d=(2, 0, 1),     F=(-4,0,4) -> -12
        #   total = -14
        #   (old (1,2,1) center -> -10; old loop never ran without
        #    tris -> center=(0,0,0) -> -22)
        expected = {
            1: 8.0 / 8.0,
            2: 8.0,
            3: -14.0,
            4: 8.0 / 8.0 / PANEL_QINF,
            5: 8.0 / 8.0 / PANEL_QINF,
        }
        self.assertEqual(sorted(values), sorted(expected))
        for sid, value in expected.items():
            self.assertAlmostEqual(values[sid], value, places=10, msg=f'sid={sid}')

    def test_no_aero_elements_raises(self):
        """no tris and no quads (e.g., region filter removed everything) must fail"""
        log = SimpleLogger(level='error')
        with self.assertRaises(RuntimeError):
            get_aero_pressure_centroid(
                _two_quad_fluent(log), 'fluent', 'force_moment',
                regions_to_include=[99])


def _tri_quad_fluent(log: SimpleLogger) -> Fluent:
    """
    mixed mesh; the quad is listed first in the results (eid 1 < eid 2)
    quad: eid=1, z=0,  2x1 -> area=2,   n=+z, Cp=1
    tri:  eid=2, z=5,  1x1 -> area=0.5, n=+z, Cp=3
    """
    node_id = np.arange(1, 8, dtype='int32')
    xyz = np.array([
        [0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.],
        [0., 0., 5.], [1., 0., 5.], [0., 1., 5.],
    ])
    quads = np.array([[1, 3, 1, 2, 3, 4]], dtype='int32')
    tris = np.array([[2, 3, 5, 6, 7]], dtype='int32')
    titles = np.array(['ElementID', 'Pressure Coefficient'])
    return Fluent.from_data(
        node_id, xyz, tris, quads,
        np.array([1, 2], dtype='int32'), titles,
        np.array([[1.0]]), np.array([[3.0]]),
        auto_read_write_h5=False, log=log, debug=False)


class TestAeroCombinedArrays(unittest.TestCase):
    def test_fluent_tri_quad_rows_align(self):
        """area/Cp/centroid/normal must all be stacked [tri, quad]"""
        log = SimpleLogger(level='warning')
        aero_dict = get_aero_pressure_centroid(
            _tri_quad_fluent(log), 'fluent', 'pressure')
        # row 0 = tri (z=5), row 1 = quad (z=0)
        np.testing.assert_allclose(aero_dict['centroid'][:, 2], [5., 0.])
        np.testing.assert_allclose(aero_dict['area'], [0.5, 2.0])
        np.testing.assert_allclose(aero_dict['Cp_centroid'], [3.0, 1.0])
        np.testing.assert_allclose(aero_dict['normal'], [[0., 0., 1.], [0., 0., 1.]])


def _capture_log() -> tuple[SimpleLogger, list[str]]:
    """info-level logger whose messages are appended to the returned list"""
    msgs = []
    log = SimpleLogger(level='info',
                       log_func=lambda typ, filename, lineno, msg: msgs.append(msg))
    return log, msgs


def _parse_balance(msgs: list[str], title: str) -> dict[str, np.ndarray]:
    """
    Reads the dimensional rows of the logged force/moment balance table.

    Returns {'aero_f', 'aero_m', 'struct_f', 'struct_m'}; values are
    printed with 6 significant figures.
    """
    lines = '\n'.join(msgs).splitlines()
    istart = [i for i, line in enumerate(lines) if line.strip() == title]
    assert len(istart) == 1, f'{title!r} found {len(istart)} times'
    labels = {
        'Aero  Force  (dim)': 'aero_f',
        'Struct Force (dim)': 'struct_f',
        'Aero  Moment (dim)': 'aero_m',
        'Struct Moment(dim)': 'struct_m',
    }
    out = {}
    for line in lines[istart[0]:]:
        sline = line.strip()
        for label, key in labels.items():
            if sline.startswith(label) and key not in out:
                out[key] = np.array(sline[len(label):].split(), dtype='float64')
        if len(out) == len(labels):
            break
    assert len(out) == len(labels), out
    return out


BALANCE_PRESSURE = 'Force/Moment Balance: Aero vs Mapped Structure (full_model)'
BALANCE_FORCE = 'Force/Moment Balance: Aero vs Mapped Structure (full_model, force)'
BALANCE_FORCE_PRESSURE = 'Force/Moment Balance: Aero vs Mapped Structure (full_model, force+pressure)'


def _plate_fin_structure(log: SimpleLogger, flip_plate: bool=False) -> BDF:
    """
    eid=1: 4x2 plate at z=0, n=+z (n=-z if flip_plate); nodes 1-4
    eid=2: 2x2 fin in the y=1 plane, x=10..12, z=0..2; nodes 5-8
    """
    model = BDF(log=log)
    xyzs = [(0., 0., 0.), (4., 0., 0.), (4., 2., 0.), (0., 2., 0.),
            (10., 1., 0.), (12., 1., 0.), (12., 1., 2.), (10., 1., 2.)]
    for nid, xyz in enumerate(xyzs, start=1):
        model.add_grid(nid, list(xyz))
    model.add_cquad4(1, 10, [1, 4, 3, 2] if flip_plate else [1, 2, 3, 4])
    model.add_cquad4(2, 10, [5, 6, 7, 8])
    model.add_pshell(10, mid1=100, t=0.1)
    model.add_mat1(100, 1.0e7, None, 0.3)
    model.cross_reference()
    return model


def _skin_fin_cart3d(log: SimpleLogger) -> Cart3D:
    """
    skin tri:         area=1,   n=+z, Cp=1, centroid=(5/3, 5/6, 0.1)
                      -> F = qinf*[0,0,1] = [0,0,2]
    fin leading edge: area=0.5, n=-x, Cp=2, centroid=(10, 5/6, 5/6)
                      -> F = qinf*[-1,0,0] = [-2,0,0]
    """
    model = Cart3D(log=log)
    model.points = np.array([
        [1., 0.5, 0.1], [3., 0.5, 0.1], [1., 1.5, 0.1],
        [10., 0.5, 0.5], [10., 0.5, 1.5], [10., 1.5, 0.5],
    ])
    model.regions = np.ones(2, dtype='int32')
    model.elements = np.array([[0, 1, 2],
                               [3, 4, 5]], dtype='int32')
    model.loads = {'Cp': np.array([1.0, 2.0])}
    return model
LOG_TOL = dict(rtol=1e-5, atol=1e-9)  # table prints 6 significant figures


class TestFullModelBalance(unittest.TestCase):
    """full_model logs aero vs structure force/moment about reference_point"""

    def _balance(self, aero_model, reference_point=None):
        log, msgs = _capture_log()
        aero_model.log = log
        structure_model = _one_quad_structure(log)
        out_model = pressure_map_to_structure_model(
            aero_model, structure_model, np.array([1], dtype='int32'),
            reference_point, map_type='pressure', qinf=PANEL_QINF,
            sref=4.0, bref=2.0, cref=0.5)
        text = '\n'.join(msgs)
        self.assertIn('sref=4.0; bref=2.0; cref=0.5', text)
        bal = _parse_balance(msgs, BALANCE_PRESSURE)
        return out_model, bal['aero_f'], bal['aero_m'], bal['struct_f'], bal['struct_m']

    def test_full_model_balance_tri_aero(self):
        log = SimpleLogger(level='warning')
        out_model, aero_f, aero_m, struct_f, struct_m = self._balance(
            _two_tri_cart3d(log))

        # aero: F_flat = 2*1*1*[0,0,1] = [0,0,2] at (7/6, 5/6, 1/2)
        #       F_tilt = 2*2*(sqrt2/2)*[-1,0,1]/sqrt2 = [-2,0,2] at (17/6, 5/6, 5/6)
        #       M = r x F = [5/3, -7/3, 0] + [5/3, -22/3, 5/3]
        np.testing.assert_allclose(aero_f, [-2., 0., 4.], **LOG_TOL)
        np.testing.assert_allclose(aero_m, [10/3, -29/3, 5/3], **LOG_TOL)

        # structure: centroid (2,1,0) is nearest the flat tri -> p = 2*1 = 2
        #   F = 2*8*[0,0,1] = [0,0,16];  M = (2,1,0) x F = [16,-32,0]
        #   (same as pyNastran's sum_forces_moments on the PLOAD2)
        self.assertAlmostEqual(_pload2_values(out_model)[1], 2.0)
        np.testing.assert_allclose(struct_f, [0., 0., 16.], **LOG_TOL)
        np.testing.assert_allclose(struct_m, [16., -32., 0.], **LOG_TOL)

    def test_full_model_balance_reference_point(self):
        """M(ref) = M(0) - ref x F on both sides"""
        log = SimpleLogger(level='warning')
        ref = np.array([1., -2., 3.])
        unused, aero_f0, aero_m0, struct_f0, struct_m0 = self._balance(
            _two_tri_cart3d(log))
        unused, aero_f, aero_m, struct_f, struct_m = self._balance(
            _two_tri_cart3d(log), reference_point=ref)
        np.testing.assert_allclose(aero_f, aero_f0)
        np.testing.assert_allclose(struct_f, struct_f0)
        np.testing.assert_allclose(aero_m, aero_m0 - np.cross(ref, aero_f0), **LOG_TOL)
        np.testing.assert_allclose(struct_m, struct_m0 - np.cross(ref, struct_f0), **LOG_TOL)

    def test_full_model_balance_logs(self):
        """the table actually reaches the log"""
        log, msgs = _capture_log()
        pressure_map_to_structure_model(
            _two_tri_cart3d(log), _one_quad_structure(log),
            np.array([1], dtype='int32'), None, qinf=PANEL_QINF)
        text = '\n'.join(msgs)
        self.assertIn(BALANCE_PRESSURE, text)
        self.assertIn('Struct Moment(CM=M/qSL)', text)


class TestFullModelForce(unittest.TestCase):
    """full_model map_type='force': aero force vectors -> nearest shell node"""

    def _map(self, structure_model, aero_model, eids=(1,)):
        log, msgs = _capture_log()
        structure_model.log = log
        out_model = pressure_map_to_structure_model(
            aero_model, structure_model, np.array(eids, dtype='int32'),
            None, map_type='force', qinf=PANEL_QINF, force_sid=2)
        text = '\n'.join(msgs)
        self.assertEqual(text.count(BALANCE_FORCE), 1)
        return out_model, _parse_balance(msgs, BALANCE_FORCE), text

    def _forces(self, out_model, sid=2) -> dict[int, np.ndarray]:
        self.assertEqual(list(out_model.loads), [sid])
        return {load.node: load.scaled_vector for load in out_model.loads[sid]}

    def test_force_vectors_and_balance(self):
        """drag (Fx) is kept and the total force is conserved exactly"""
        log = SimpleLogger(level='warning')
        out_model, bal, text = self._map(
            _one_quad_structure(log), _two_tri_cart3d(log))

        # flat tri   (7/6, 5/6, 1/2) -> nearest node 1 (0,0,0): F=[0,0,2]
        # tilted tri (17/6,5/6, 5/6) -> nearest node 2 (4,0,0): F=[-2,0,2]
        forces = self._forces(out_model)
        self.assertEqual(sorted(forces), [1, 2])
        np.testing.assert_allclose(forces[1], [0., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(forces[2], [-2., 0., 2.], atol=1e-12)
        for load in out_model.loads[2]:
            self.assertEqual(load.type, 'FORCE')
            self.assertAlmostEqual(np.linalg.norm(load.xyz), 1.0)

        # exact conservation from the FORCE cards themselves
        np.testing.assert_allclose(sum(forces.values()), [-2., 0., 4.], atol=1e-12)

        # logged balance table (6 significant figures)
        np.testing.assert_allclose(bal['aero_f'], [-2., 0., 4.], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        np.testing.assert_allclose(bal['aero_m'], [10/3, -29/3, 5/3], **LOG_TOL)
        # moment at the nodes: (4,0,0) x (-2,0,2) = (0,-8,0)
        np.testing.assert_allclose(bal['struct_m'], [0., -8., 0.], **LOG_TOL)

    def test_force_only_on_mapped_element_nodes(self):
        """a closer node on an unmapped element (e.g. a bulkhead) gets nothing"""
        log = SimpleLogger(level='warning')
        structure_model = _one_quad_structure(log)
        # small quad (eid=2) right at the flat tri centroid; not mapped
        xyzs = [(1.1, 0.8, 0.5), (1.3, 0.8, 0.5), (1.3, 0.9, 0.5), (1.1, 0.9, 0.5)]
        for nid, xyz in enumerate(xyzs, start=5):
            structure_model.add_grid(nid, list(xyz))
        structure_model.add_cquad4(2, 10, [5, 6, 7, 8])
        structure_model.cross_reference()

        out_model, unused_args, unused_text = self._map(
            structure_model, _two_tri_cart3d(log), eids=(1,))
        self.assertEqual(sorted(self._forces(out_model)), [1, 2])

        out_model, unused_args, unused_text = self._map(
            structure_model, _two_tri_cart3d(log), eids=(1, 2))
        forces = self._forces(out_model)
        self.assertIn(5, forces)  # now it's mapped, the flat tri goes there
        np.testing.assert_allclose(sum(forces.values()), [-2., 0., 4.], atol=1e-12)

    def _map_tol(self, structure_model, aero_model, tol, eids=(1,)):
        log, msgs = _capture_log()
        structure_model.log = log
        out_model = pressure_map_to_structure_model(
            aero_model, structure_model, np.array(eids, dtype='int32'),
            None, map_type='force', qinf=PANEL_QINF, force_sid=2,
            force_to_pressure_tol=tol)
        # 2: force (FORCE + PLOAD2); 5: Cp_signed (only if there are PLOAD2s)
        self.assertIn(sorted(out_model.loads), ([2], [2, 5]))
        loads = out_model.loads[2]
        pressures = {load.eids[0]: load.pressure for load in loads if load.type == 'PLOAD2'}
        forces = {load.node: load.scaled_vector for load in loads if load.type == 'FORCE'}
        bal = _parse_balance(msgs, BALANCE_FORCE_PRESSURE)
        return pressures, forces, bal, '\n'.join(msgs)

    def test_force_to_pressure_skin_and_fin(self):
        """skin force (normal) -> PLOAD2; fin leading-edge drag (in-plane) -> FORCE"""
        log = SimpleLogger(level='warning')
        pressures, forces, bal, text = self._map_tol(
            _plate_fin_structure(log), _skin_fin_cart3d(log), tol=0.01, eids=(1, 2))

        # skin: F=[0,0,2] along +z on the 8 in^2 plate -> p = 2/8
        self.assertEqual(sorted(pressures), [1])
        self.assertAlmostEqual(pressures[1], 0.25)
        # fin: F=[-2,0,0] is in the fin plane -> FORCE at node 5 (10,1,0)
        self.assertEqual(sorted(forces), [5])
        np.testing.assert_allclose(forces[5], [-2., 0., 0.], atol=1e-12)

        # nothing dropped -> force conserved
        np.testing.assert_allclose(bal['aero_f'], [-2., 0., 2.], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        self.assertIn('1 of 2 elements -> PLOAD2', text)

    def test_force_to_pressure_flipped_normal(self):
        """the PLOAD2 sign follows F.n, so a reversed shell still pushes the right way"""
        log = SimpleLogger(level='warning')
        pressures, unused_forces, bal, unused_text = self._map_tol(
            _plate_fin_structure(log, flip_plate=True), _skin_fin_cart3d(log),
            tol=0.01, eids=(1, 2))
        self.assertAlmostEqual(pressures[1], -0.25)
        np.testing.assert_allclose(bal['struct_f'], [-2., 0., 2.], **LOG_TOL)

    def _cp_signed(self, structure_model, aero_model, tol, eids=(1, 2)):
        structure_model.log = SimpleLogger(level='warning')
        out_model = pressure_map_to_structure_model(
            aero_model, structure_model, np.array(eids, dtype='int32'),
            None, map_type='force', qinf=PANEL_QINF, force_sid=2,
            cp_signed_sid=5, force_to_pressure_tol=tol)
        self.assertEqual(sorted(out_model.loads), [2, 5])
        return {load.eids[0]: load.pressure for load in out_model.loads[5]}

    def test_cp_signed_ignores_structure_normal(self):
        """skin tri Cp=1: Cp_signed=+Cp/area ratio either way the plate points"""
        log = SimpleLogger(level='warning')
        # F=[0,0,2] over 8 in^2 -> p=+-0.25 -> Cp_signed = 0.25/qinf
        for flip_plate in [False, True]:
            cp_signed = self._cp_signed(
                _plate_fin_structure(log, flip_plate=flip_plate),
                _skin_fin_cart3d(log), tol=0.01)
            self.assertEqual(sorted(cp_signed), [1], msg=f'flip_plate={flip_plate}')
            self.assertAlmostEqual(cp_signed[1], 0.25 / PANEL_QINF,
                                   msg=f'flip_plate={flip_plate}')

    def test_cp_signed_negative_cp(self):
        """a suction Cp stays negative"""
        log = SimpleLogger(level='warning')
        aero_model = _skin_fin_cart3d(log)
        aero_model.loads['Cp'] = np.array([-1.0, 2.0])
        for flip_plate in [False, True]:
            cp_signed = self._cp_signed(
                _plate_fin_structure(log, flip_plate=flip_plate), aero_model, tol=0.01)
            self.assertAlmostEqual(cp_signed[1], -0.25 / PANEL_QINF,
                                   msg=f'flip_plate={flip_plate}')

    def test_force_to_pressure_mixed_element_sizes(self):
        """
        fine aero mesh, uniform Cp=1, structure = one 4x1 quad + two 0.5x1 quads

        nearest *centroid* would give the small quad at x=4..4.5 every aero
        element with x in (3.125, 4.5) -> p = 2.8*qinf (and the big quad
        0.775*qinf).  The closest *surface* gives each quad its footprint,
        so p = qinf*Cp everywhere and the force is conserved.
        """
        log = SimpleLogger(level='warning')
        structure_model = BDF(log=log)
        xs = [0., 4., 4.5, 5.]
        nid = 1
        for x in xs:
            structure_model.add_grid(nid, [x, 0., 0.])
            structure_model.add_grid(nid + 1, [x, 1., 0.])
            nid += 2
        for eid in range(1, 4):
            n1 = 2 * eid - 1
            structure_model.add_cquad4(eid, 10, [n1, n1 + 2, n1 + 3, n1 + 1])
        structure_model.add_pshell(10, mid1=100, t=0.1)
        structure_model.add_mat1(100, 1.0e7, None, 0.3)
        structure_model.cross_reference()

        # 0.05 aero grid on z=0.01; no tri centroid lies on x=4 or x=4.5
        nx, ny = 100, 20
        x = np.linspace(0., 5., nx + 1)
        y = np.linspace(0., 1., ny + 1)
        xg, yg = np.meshgrid(x, y, indexing='ij')
        points = np.column_stack([xg.ravel(), yg.ravel(), np.full(xg.size, 0.01)])
        tris = []
        for i in range(nx):
            for j in range(ny):
                n1 = i * (ny + 1) + j
                n2 = n1 + (ny + 1)
                tris.extend([(n1, n2, n2 + 1), (n1, n2 + 1, n1 + 1)])
        ntri = len(tris)
        aero_model = Cart3D(log=log)
        aero_model.points = points
        aero_model.regions = np.ones(ntri, dtype='int32')
        aero_model.elements = np.array(tris, dtype='int32')
        aero_model.loads = {'Cp': np.ones(len(tris))}

        pressures, forces, bal, unused_text = self._map_tol(
            structure_model, aero_model, tol=0.01, eids=(1, 2, 3))
        self.assertEqual(forces, {})
        self.assertEqual(sorted(pressures), [1, 2, 3])
        for eid, p in pressures.items():
            self.assertAlmostEqual(p, PANEL_QINF, places=10, msg=f'eid={eid}')
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)

    def test_closest_point_triangles(self):
        """vectorized point-triangle distance vs brute-force sampling"""
        from pyNastran.dev.tools.pressure_map.pressure_map import _closest_point_triangles
        rng = np.random.default_rng(0)
        n = 200
        tri = rng.normal(size=(n, 3, 3))
        pts = 2. * rng.normal(size=(n, 3))
        dist = _closest_point_triangles(pts, tri)
        # barycentric grid on each triangle (includes edges/vertices)
        m = 60
        u, v = np.meshgrid(np.linspace(0, 1, m + 1), np.linspace(0, 1, m + 1))
        keep = (u + v) <= 1.0
        u = u[keep]
        v = v[keep]
        for i in range(n):
            a, b, c = tri[i]
            samples = a + np.outer(u, b - a) + np.outer(v, c - a)
            brute = np.linalg.norm(samples - pts[i], axis=1).min()
            self.assertLessEqual(dist[i], brute + 1e-12)
            self.assertGreater(dist[i], brute - 0.05 * np.linalg.norm(b - a) - 0.05 * np.linalg.norm(c - a))

    def test_force_to_pressure_tolerance(self):
        """
        both tris -> the one plate: F=[-2,0,4], p=4/8, in-plane r=[-2,0,0]
        |r|/|F| = 2/sqrt(20) = 0.447
        """
        log = SimpleLogger(level='warning')
        # within tol: PLOAD2; the in-plane part is dropped and logged
        pressures, forces, bal, text = self._map_tol(
            _one_quad_structure(log), _two_tri_cart3d(log), tol=0.5)
        self.assertEqual(pressures, {1: 0.5})
        self.assertEqual(forces, {})
        np.testing.assert_allclose(bal['struct_f'], [0., 0., 4.], **LOG_TOL)
        self.assertIn('dropped in-plane force = [-2.', text)

        # outside tol: same FORCE cards as force_to_pressure_tol=None
        pressures, forces, bal, unused_text = self._map_tol(
            _one_quad_structure(log), _two_tri_cart3d(log), tol=0.4)
        self.assertEqual(pressures, {})
        self.assertEqual(sorted(forces), [1, 2])
        np.testing.assert_allclose(forces[1], [0., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(forces[2], [-2., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(bal['struct_f'], [-2., 0., 4.], **LOG_TOL)

    def test_force_to_pressure_tol_ignored_for_pressure(self):
        """map_type='pressure' warns and ignores force_to_pressure_tol"""
        log, msgs = _capture_log()
        out_model = pressure_map_to_structure_model(
            _two_tri_cart3d(log), _one_quad_structure(log),
            np.array([1], dtype='int32'), None, map_type='pressure',
            qinf=PANEL_QINF, force_to_pressure_tol=0.05, cp_sid=0)
        self.assertIn("force_to_pressure_tol=0.05 requires map_type='force'", '\n'.join(msgs))
        self.assertEqual(_pload2_values(out_model), {1: 2.0})

    def test_force_to_pressure_tol_invalid_raises(self):
        """unused or not, a bad force_to_pressure_tol is an error"""
        log = SimpleLogger(level='error')
        for map_type in ['pressure', 'force']:
            for tol, error in [('0.05', TypeError), (True, TypeError),
                               (-0.1, ValueError), (np.nan, ValueError)]:
                with self.assertRaises(error, msg=f'map_type={map_type} tol={tol!r}'):
                    pressure_map_to_structure_model(
                        _two_tri_cart3d(log), _one_quad_structure(log),
                        np.array([1], dtype='int32'), None, map_type=map_type,
                        force_to_pressure_tol=tol)

    def test_force_moment_map_type_raises(self):
        log = SimpleLogger(level='warning')
        with self.assertRaises(NotImplementedError):
            pressure_map_to_structure_model(
                _two_tri_cart3d(log), _one_quad_structure(log),
                np.array([1], dtype='int32'), None, map_type='force_moment')


class TestAeroXyzScale(unittest.TestCase):
    """aero_xyz_scale is the only unit conversion; check it end to end"""

    def test_get_aero_model_scales_in_place(self):
        """a passed-in model is scaled in place (scale exactly once)"""
        log = SimpleLogger(level='warning')
        aero_model = _two_tri_cart3d(log)
        points0 = aero_model.points.copy()
        model2, unused_variables = get_aero_model(
            aero_model, 'cart3d', aero_xyz_scale=2.0, log=log)
        self.assertIs(model2, aero_model)
        np.testing.assert_allclose(aero_model.points, 2.0 * points0)

    def test_get_aero_model_force_moment(self):
        """get_aero_model logs F/M about reference_point before returning"""
        log, msgs = _capture_log()
        get_aero_model(_two_tri_cart3d(log), 'cart3d',
                       qinf=PANEL_QINF, xyz_units_out='in', pressure_units='psi')
        text = '\n'.join(msgs)
        self.assertIn('aero force/moment sum (all 2 aero elements)', text)
        # same totals as TestFullModelBalance
        self.assertIn('F = [-2, 0, 4]', text)
        self.assertIn('M = [3.33333, -9.66667, 1.66667]', text)

    def test_aero_force_moment_scaled_and_reference_point(self):
        """sums use the scaled xyz; F ~ s^2, M ~ s^3; M(ref) = M(0) - ref x F"""
        log = SimpleLogger(level='warning')
        f0, m0 = log_aero_force_moment(_two_tri_cart3d(log), 'cart3d', qinf=PANEL_QINF)

        s = 2.0
        aero_model = _two_tri_cart3d(log)
        get_aero_model(aero_model, 'cart3d', aero_xyz_scale=s, qinf=PANEL_QINF, log=log)
        f, m = log_aero_force_moment(aero_model, 'cart3d', qinf=PANEL_QINF)
        np.testing.assert_allclose(f, s**2 * f0, atol=1e-12)
        np.testing.assert_allclose(m, s**3 * m0, atol=1e-12)

        ref = np.array([1., -2., 3.])
        f_ref, m_ref = log_aero_force_moment(
            _two_tri_cart3d(log), 'cart3d', qinf=PANEL_QINF, reference_point=ref)
        np.testing.assert_allclose(f_ref, f0, atol=1e-12)
        np.testing.assert_allclose(m_ref, m0 - np.cross(ref, f0), atol=1e-12)

    def test_aero_force_moment_fluent(self):
        """Fluent (quads), including the already-scaled early-return path"""
        log, msgs = _capture_log()
        f, m = log_aero_force_moment(_two_quad_fluent(log), 'fluent', qinf=PANEL_QINF)
        # F/q = [0,0,2] + [-2,0,2]
        np.testing.assert_allclose(f, [-4., 0., 8.], atol=1e-12)

        msgs.clear()
        get_aero_model(_two_quad_fluent(log), 'fluent', qinf=PANEL_QINF)
        self.assertIn('F = [-4, 0, 8]', '\n'.join(msgs))

    def test_aero_force_moment_no_cp_skips(self):
        log, msgs = _capture_log()
        aero_model = _two_tri_cart3d(log)
        aero_model.loads = {'rho': np.ones(2)}
        get_aero_model(aero_model, 'cart3d')
        self.assertIn("aero force/moment sum: no 'Cp' result; skipping", '\n'.join(msgs))

    def test_pressure_map_aero_xyz_scale(self):
        """aero in m, structure in in -> same PLOAD2s as an all-inch model"""
        log = SimpleLogger(level='warning')
        m_to_in = 39.3701

        def _run(aero_model, aero_xyz_scale: float) -> dict[int, float]:
            out_model = pressure_map(
                aero_model, _one_quad_structure(log),
                eids_structure=np.array([1], dtype='int32'),
                eid_load_id=0,
                aero_format='cart3d',
                map_type='force_moment',
                method='panel_model',
                xyz_units_out='in', pressure_units='psi',
                pressure_filename='',  # don't write files
                aero_xyz_scale=aero_xyz_scale,
                qinf=PANEL_QINF,
                log=log)
            return _pload2_values(out_model)

        values_in = _run(_two_tri_cart3d(log), 1.0)
        values_m = _run(_two_tri_cart3d(log, xyz_scale=1.0 / m_to_in), m_to_in)
        self.assertEqual(sorted(values_in), sorted(values_m))
        for sid, value in values_in.items():
            self.assertAlmostEqual(values_m[sid], value, places=8, msg=f'sid={sid}')


BALANCE_LCP_GROUP = 'Force/Moment Balance: Aero vs Mapped Structure (LCP; group={})'
BALANCE_LCP = BALANCE_LCP_GROUP.format('total')


def _lcp_structure(log: SimpleLogger,
                   xyzs: list[tuple[float, float, float]]) -> BDF:
    """GRIDs only (no elements), nid = 1, 2, ..."""
    model = BDF(log=log)
    model.bdf_filename = 'lcp_nodes.bdf'
    for nid, xyz in enumerate(xyzs, start=1):
        model.add_grid(nid, list(xyz))
    model.cross_reference()
    return model


class TestLCP(unittest.TestCase):
    """
    method='LCP': each aero tri/quad force qinf*Cp*A*n goes to the nearest
    load collection point (node).  The node also gets a MOMENT card with the
    transfer moment sum((c_i - x_node) x F_i), so FORCE + MOMENT reproduce
    the aero force and moment exactly.

    _two_tri_cart3d (qinf=PANEL_QINF=2):
      flat tri:   c=(7/6, 5/6, 1/2),  F=(0, 0, 2)
      tilted tri: c=(17/6, 5/6, 5/6), F=(-2, 0, 2)
    """

    def _map(self, structure_model, aero_model, qinf=PANEL_QINF,
             reference_point=None):
        log, msgs = _capture_log()
        structure_model.log = log
        if reference_point is None:
            reference_point = np.zeros(3)
        out_model = pressure_map_to_lcp_model(
            aero_model, structure_model, np.array([], dtype='int32'),
            reference_point, {}, qinf=qinf, force_sid=2, moment_sid=3)
        text = '\n'.join(msgs)
        self.assertEqual(text.count(BALANCE_LCP), 1)
        return out_model, _parse_balance(msgs, BALANCE_LCP)

    def _map_groups(self, structure_model, aero_model, load_control_nids):
        """
        Returns
        -------
        out_model : BDF
        bals : {group_name: balance}
        text : str
            the full log (for warnings)
        """
        log, msgs = _capture_log()
        structure_model.log = log
        out_model = pressure_map_to_lcp_model(
            aero_model, structure_model, np.array([], dtype='int32'),
            np.zeros(3), load_control_nids,
            qinf=PANEL_QINF, force_sid=2, moment_sid=3)
        text = '\n'.join(msgs)
        bals = {}
        for group_name in load_control_nids:
            title = BALANCE_LCP_GROUP.format(group_name)
            if title in text:
                bals[group_name] = _parse_balance(msgs, title)
        return out_model, bals, text

    def _cards(self, out_model) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
        """({nid: force}, {nid: moment}); one FORCE and one MOMENT per node"""
        self.assertEqual(sorted(out_model.loads), [2, 3])
        forces = {}
        moments = {}
        for load in out_model.loads[2]:
            self.assertEqual(load.type, 'FORCE')
            self.assertNotIn(load.node, forces)
            forces[load.node] = load.scaled_vector
        for load in out_model.loads[3]:
            self.assertEqual(load.type, 'MOMENT')
            self.assertNotIn(load.node, moments)
            moments[load.node] = load.scaled_vector
        self.assertEqual(sorted(forces), sorted(moments))
        return forces, moments

    def _total_about(self, structure_model, forces, moments, point):
        """total F and M about point from the FORCE/MOMENT cards"""
        total_f = np.zeros(3)
        total_m = np.zeros(3)
        for nid, force in forces.items():
            xyz = structure_model.nodes[nid].get_position()
            total_f += force
            total_m += moments[nid] + np.cross(xyz - point, force)
        return total_f, total_m

    def test_lcp_nearest_node_force_moment(self):
        """each panel goes to its own nearest node; the far node gets nothing"""
        log = SimpleLogger(level='warning')
        # nid=1 is far away and listed first, so a node/force pairing
        # mix-up would put a load on it
        structure_model = _lcp_structure(
            log, [(100., 100., 100.), (1., 1., 0.), (3., 1., 1.)])
        out_model, bal = self._map(structure_model, _two_tri_cart3d(log))
        forces, moments = self._cards(out_model)

        self.assertEqual(sorted(forces), [2, 3])
        # flat tri -> nid 2 (1,1,0); d=(1/6, -1/6, 1/2); d x (0,0,2)
        np.testing.assert_allclose(forces[2], [0., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(moments[2], [-1/3, -1/3, 0.], atol=1e-12)
        # tilted tri -> nid 3 (3,1,1); d=(-1/6, -1/6, -1/6); d x (-2,0,2)
        np.testing.assert_allclose(forces[3], [-2., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(moments[3], [-1/3, 2/3, -1/3], atol=1e-12)

        # aero about the origin: F=(-2,0,4), M=(10/3, -29/3, 5/3)
        np.testing.assert_allclose(bal['aero_f'], [-2., 0., 4.], **LOG_TOL)
        np.testing.assert_allclose(bal['aero_m'], [10/3, -29/3, 5/3], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)

    def test_lcp_all_panels_to_one_node(self):
        """every panel is summed once (no duplicates) onto a single node"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(10., 0., 0.)])
        out_model, bal = self._map(structure_model, _two_tri_cart3d(log))
        forces, moments = self._cards(out_model)

        self.assertEqual(sorted(forces), [1])
        np.testing.assert_allclose(forces[1], [-2., 0., 4.], atol=1e-12)
        # M about (10,0,0) = M_origin - (10,0,0) x F
        #                  = (10/3, -29/3, 5/3) - (0, -40, 0)
        np.testing.assert_allclose(moments[1], [10/3, 91/3, 5/3], atol=1e-12)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)

    def test_lcp_quads(self):
        """quad-only fluent model: F=(-4,0,8) at qinf=2"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(10., 0., 0.)])
        out_model, bal = self._map(structure_model, _two_quad_fluent(log))
        forces, moments = self._cards(out_model)
        np.testing.assert_allclose(forces[1], [-4., 0., 8.], atol=1e-12)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)

    def test_lcp_tri_and_quad(self):
        """mixed fluent mesh: tri and quad both mapped, each once"""
        log = SimpleLogger(level='warning')
        # quad: c=(1, 0.5, 0), A=2, Cp=1;  tri: c=(1/3, 1/3, 5), A=0.5, Cp=3
        structure_model = _lcp_structure(log, [(1., 0.5, 0.), (0., 0., 5.)])
        out_model, bal = self._map(structure_model, _tri_quad_fluent(log))
        forces, moments = self._cards(out_model)
        self.assertEqual(sorted(forces), [1, 2])
        np.testing.assert_allclose(forces[1], [0., 0., PANEL_QINF * 2.0], atol=1e-12)
        np.testing.assert_allclose(forces[2], [0., 0., PANEL_QINF * 1.5], atol=1e-12)
        np.testing.assert_allclose(moments[1], [0., 0., 0.], atol=1e-12)
        np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
        np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)

    def test_lcp_qinf_scaling(self):
        """FORCE/MOMENT cards are linear in qinf"""
        log = SimpleLogger(level='warning')
        xyzs = [(1., 1., 0.), (3., 1., 1.)]
        out1, unused_bal1 = self._map(_lcp_structure(log, xyzs), _two_tri_cart3d(log), qinf=1.0)
        out5, unused_bal5 = self._map(_lcp_structure(log, xyzs), _two_tri_cart3d(log), qinf=5.0)
        forces1, moments1 = self._cards(out1)
        forces5, moments5 = self._cards(out5)
        self.assertEqual(sorted(forces1), sorted(forces5))
        for nid in forces1:
            np.testing.assert_allclose(forces5[nid], 5.0 * forces1[nid], atol=1e-12)
            np.testing.assert_allclose(moments5[nid], 5.0 * moments1[nid], atol=1e-12)

    def test_lcp_reference_point(self):
        """the cards don't depend on the reference point; the totals do"""
        log = SimpleLogger(level='warning')
        xyzs = [(0., 0., 0.), (2., 1., 0.), (4., 1., 1.)]
        ref = np.array([1., -2., 3.])
        structure_model = _lcp_structure(log, xyzs)
        out0, unused_bal0 = self._map(_lcp_structure(log, xyzs), _two_tri_cart3d(log))
        out_ref, bal_ref = self._map(structure_model, _two_tri_cart3d(log),
                                     reference_point=ref)
        forces0, moments0 = self._cards(out0)
        forces_ref, moments_ref = self._cards(out_ref)
        self.assertEqual(sorted(forces0), sorted(forces_ref))
        for nid in forces0:
            np.testing.assert_allclose(forces_ref[nid], forces0[nid], atol=1e-12)
            np.testing.assert_allclose(moments_ref[nid], moments0[nid], atol=1e-12)

        # aero about ref = M_origin - ref x F
        aero_f = np.array([-2., 0., 4.])
        aero_m = np.array([10/3, -29/3, 5/3]) - np.cross(ref, aero_f)
        total_f, total_m = self._total_about(structure_model, forces_ref, moments_ref, ref)
        np.testing.assert_allclose(total_f, aero_f, atol=1e-12)
        np.testing.assert_allclose(total_m, aero_m, atol=1e-12)
        np.testing.assert_allclose(bal_ref['aero_m'], aero_m, **LOG_TOL)
        np.testing.assert_allclose(bal_ref['struct_m'], aero_m, **LOG_TOL)

    def test_lcp_pressure_map_writes_files(self):
        """pressure_map(method='LCP') writes the main deck, loads, and csv"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(
            log, [(100., 100., 100.), (1., 1., 0.), (3., 1., 1.)])
        with tempfile.TemporaryDirectory() as dirname:
            dirname = Path(dirname)
            structure_model.bdf_filename = str(dirname / 'lcp_nodes.bdf')
            pressure_filename = dirname / 'lcp.bdf'
            pressure_map(
                _two_tri_cart3d(log), structure_model,
                aero_format='cart3d',
                map_type='force_moment',
                method='LCP',
                xyz_units_in='in', xyz_units_out='in', pressure_units='psi',
                force_sid=2, moment_sid=3,
                pressure_filename=str(pressure_filename),
                qinf=PANEL_QINF,
                log=log)
            self.assertTrue((dirname / 'main_lcp.bdf').exists())
            self.assertTrue(pressure_filename.exists())
            csv_filename = dirname / 'lcp.force_moment.csv'
            self.assertTrue(csv_filename.exists())

            # Node,Fx,Fy,Fz,Mx,My,Mz (trailing comma, no header)
            with open(csv_filename, 'r') as csv_file:
                rows = [line.strip().rstrip(',').split(',')
                        for line in csv_file if line.strip()]
            data = np.array(rows, dtype='float64')

            # the loads deck reads back with the same cards
            loads_model = BDF(log=log)
            loads_model.read_bdf(str(pressure_filename), xref=False, punch=True)
            self.assertEqual(sorted(loads_model.loads), [2, 3])

            # total: LOAD 6 = 1.0*FORCE(2) + 1.0*MOMENT(3), and SUBCASE 4 uses it
            self.assertEqual(sorted(loads_model.load_combinations), [6])
            load = loads_model.load_combinations[6][0]
            self.assertEqual(load.scale, 1.0)
            self.assertEqual(list(load.scale_factors), [1.0, 1.0])
            self.assertEqual(list(load.load_ids), [2, 3])
            main_text = (dirname / 'main_lcp.bdf').read_text()
            self.assertIn('SUBCASE 4\n  SUBTITLE = Total (Force + Moment)\n  LOAD = 6\n', main_text)
            self.assertLess(main_text.index('SUBCASE 3'), main_text.index('SUBCASE 4'))

        np.testing.assert_array_equal(data[:, 0], [2, 3])
        np.testing.assert_allclose(data[0, 1:], [0., 0., 2., -1/3, -1/3, 0.],
                                   rtol=1e-6, atol=1e-9)
        np.testing.assert_allclose(data[1, 1:], [-2., 0., 2., -1/3, 2/3, -1/3],
                                   rtol=1e-6, atol=1e-9)

    # ------------------------------------------------------------------
    # load_control_nids groups (wing/tail/...)
    # ------------------------------------------------------------------
    def test_lcp_groups(self):
        """
        each group maps only its own panels to only its own nodes

        node 1 (3,1,1) is nearest the tilted tri and node 2 (1,1,0) is
        nearest the flat tri; the groups force the opposite assignment
        """
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(3., 1., 1.), (1., 1., 0.)])
        groups = {
            'wing': {'regions': [1], 'nodes': [1]},
            'tail': {'regions': [2], 'nodes': [2]},
        }
        out_model, bals, text = self._map_groups(
            structure_model, _two_region_cart3d(log), groups)
        self.assertNotIn('missing regions', text)
        self.assertNotIn('missing nodes', text)
        self.assertNotIn('LCP aero regions', text)
        forces, moments = self._cards(out_model)
        self.assertEqual(sorted(forces), [1, 2])

        # the first FORCE/MOMENT card of each group is labeled with the group
        with tempfile.TemporaryDirectory() as dirname:
            bdf_filename = Path(dirname) / 'loads.bdf'
            out_model.write_bdf(bdf_filename)
            bdf_text = bdf_filename.read_text()
        for word in ('Force', 'Moment'):
            for group_name in groups:
                self.assertEqual(bdf_text.count(f'${word}; group={group_name}; qinf='), 1,
                                 (word, group_name))

        # flat tri -> node 1; d=(-11/6, -1/6, -1/2); d x (0,0,2)
        np.testing.assert_allclose(forces[1], [0., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(moments[1], [-1/3, 11/3, 0.], atol=1e-12)
        # tilted tri -> node 2; d=(11/6, -1/6, 5/6); d x (-2,0,2)
        np.testing.assert_allclose(forces[2], [-2., 0., 2.], atol=1e-12)
        np.testing.assert_allclose(moments[2], [-1/3, -16/3, -1/3], atol=1e-12)

        # one balance per group, against that group's aero panels
        self.assertEqual(sorted(bals), ['tail', 'wing'])
        np.testing.assert_allclose(bals['wing']['aero_f'], [0., 0., 2.], **LOG_TOL)
        np.testing.assert_allclose(bals['wing']['aero_m'], [5/3, -7/3, 0.], **LOG_TOL)
        np.testing.assert_allclose(bals['tail']['aero_f'], [-2., 0., 2.], **LOG_TOL)
        np.testing.assert_allclose(bals['tail']['aero_m'], [5/3, -22/3, 5/3], **LOG_TOL)
        for bal in bals.values():
            np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
            np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)

    def test_lcp_groups_match_separate_calls(self):
        """grouped == one call per component (aero/structure cut by hand)"""
        log = SimpleLogger(level='warning')
        xyz_wing = (3., 1., 1.)
        xyz_tail = (1., 1., 0.)
        groups = {
            'wing': {'regions': [1], 'nodes': [1]},
            'tail': {'regions': [2], 'nodes': [2]},
        }
        out_model, unused_bals, unused_text = self._map_groups(
            _lcp_structure(log, [xyz_wing, xyz_tail]), _two_region_cart3d(log), groups)
        forces, moments = self._cards(out_model)

        # separate 'total' calls; the single node is nid=1 in each
        out_wing, unused_bal = self._map(
            _lcp_structure(log, [xyz_wing]), _two_region_cart3d(log, keep_region=1))
        out_tail, unused_bal = self._map(
            _lcp_structure(log, [xyz_tail]), _two_region_cart3d(log, keep_region=2))
        forces_wing, moments_wing = self._cards(out_wing)
        forces_tail, moments_tail = self._cards(out_tail)

        np.testing.assert_allclose(forces[1], forces_wing[1], atol=1e-12)
        np.testing.assert_allclose(moments[1], moments_wing[1], atol=1e-12)
        np.testing.assert_allclose(forces[2], forces_tail[1], atol=1e-12)
        np.testing.assert_allclose(moments[2], moments_tail[1], atol=1e-12)

    def test_lcp_groups_shared_node(self):
        """a node in two groups gets one FORCE/MOMENT per group; the sum is the total"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(10., 0., 0.)])
        groups = {
            'wing': {'regions': [1], 'nodes': [1]},
            'tail': {'regions': [2], 'nodes': [1]},
        }
        out_model, bals, unused_text = self._map_groups(
            structure_model, _two_region_cart3d(log), groups)
        force_cards = out_model.loads[2]
        moment_cards = out_model.loads[3]
        self.assertEqual([load.node for load in force_cards], [1, 1])
        self.assertEqual([load.node for load in moment_cards], [1, 1])
        total_f = sum(load.scaled_vector for load in force_cards)
        total_m = sum(load.scaled_vector for load in moment_cards)
        # same as test_lcp_all_panels_to_one_node
        np.testing.assert_allclose(total_f, [-2., 0., 4.], atol=1e-12)
        np.testing.assert_allclose(total_m, [10/3, 91/3, 5/3], atol=1e-12)
        self.assertEqual(sorted(bals), ['tail', 'wing'])

    def test_lcp_groups_missing_ids_dropped(self):
        """unknown regions/nodes are dropped with a warning; the rest still maps"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(3., 1., 1.), (1., 1., 0.)])
        groups = {
            'wing': {'regions': [1, 99], 'nodes': [1, 1000]},
            'tail': {'regions': [2], 'nodes': [2]},
        }
        out_model, bals, text = self._map_groups(
            structure_model, _two_region_cart3d(log), groups)
        self.assertIn("LCP group='wing' has 1 missing regions; [99]", text)
        self.assertIn("LCP group='wing' has 1 missing nodes; [1000]", text)
        forces, unused_moments = self._cards(out_model)
        self.assertEqual(sorted(forces), [1, 2])
        np.testing.assert_allclose(forces[1], [0., 0., 2.], atol=1e-12)
        self.assertEqual(sorted(bals), ['tail', 'wing'])

    def test_lcp_groups_unassigned_and_overlapping_regions_warn(self):
        """a region in no group or in two groups is logged"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(3., 1., 1.), (1., 1., 0.)])

        # region 2 in no group -> its load is dropped
        groups = {'wing': {'regions': [1], 'nodes': [1, 2]}}
        out_model, bals, text = self._map_groups(
            structure_model, _two_region_cart3d(log), groups)
        self.assertIn('LCP aero regions=[2] are not in any group', text)
        forces, unused_moments = self._cards(out_model)
        np.testing.assert_allclose(sum(forces.values()), [0., 0., 2.], atol=1e-12)

        # region 1 in both groups -> counted twice
        groups = {
            'wing': {'regions': [1, 2], 'nodes': [1]},
            'tail': {'regions': [1], 'nodes': [2]},
        }
        unused_out, unused_bals, text = self._map_groups(
            _lcp_structure(log, [(3., 1., 1.), (1., 1., 0.)]),
            _two_region_cart3d(log), groups)
        self.assertIn('LCP aero regions=[1] are in more than one group', text)

    def test_lcp_groups_empty_raises(self):
        """a group with nothing left after the cleanup can't be mapped"""
        log = SimpleLogger(level='error')

        # no regions left -> raises
        bad_regions_groups = [
            {'wing': {'regions': [], 'nodes': [1]}},
            {'wing': {'regions': [99], 'nodes': [1]}},
        ]
        for groups in bad_regions_groups:
            with self.assertRaises(ValueError, msg=str(groups)):
                self._map_groups(_lcp_structure(log, [(1., 1., 0.)]),
                                 _two_region_cart3d(log), groups)

        # no nodes left -> the group is skipped with a message (no raise);
        # nothing is mapped, so there are no cards
        bad_nodes_groups = [
            {'wing': {'regions': [1], 'nodes': []}},
            {'wing': {'regions': [1], 'nodes': [1000]}},
        ]
        for groups in bad_nodes_groups:
            out_model, bals, text = self._map_groups(
                _lcp_structure(log, [(1., 1., 0.)]),
                _two_region_cart3d(log), groups)
            self.assertIn('skipping group wing because nodes are empty', text, msg=str(groups))
            self.assertEqual(out_model.loads, {}, msg=str(groups))
            self.assertEqual(out_model.load_combinations, {}, msg=str(groups))
            self.assertEqual(bals, {}, msg=str(groups))

    def test_lcp_total_load_card(self):
        """the LOAD card combines force_sid and moment_sid; load_sid must be unique"""
        log = SimpleLogger(level='warning')
        structure_model = _lcp_structure(log, [(10., 0., 0.)])
        out_model, unused_bal = self._map(structure_model, _two_tri_cart3d(log))
        self.assertEqual(sorted(out_model.load_combinations), [6])
        load = out_model.load_combinations[6][0]
        self.assertEqual(list(load.load_ids), [2, 3])
        self.assertEqual(list(load.scale_factors), [1.0, 1.0])

        for load_sid in [2, 3]:
            with self.assertRaises(ValueError, msg=str(load_sid)):
                pressure_map_to_lcp_model(
                    _two_tri_cart3d(log), _lcp_structure(log, [(10., 0., 0.)]),
                    np.array([], dtype='int32'), np.zeros(3), {},
                    qinf=PANEL_QINF, force_sid=2, moment_sid=3, load_sid=load_sid)

    # ------------------------------------------------------------------
    # aero centroid vs average node location
    # ------------------------------------------------------------------
    def test_lcp_check_group_location(self):
        """
        aero: centroids (0,0,0) A=1 and (4,0,0) A=3
          -> area-weighted centroid (3,0,0); size = bbox diagonal = 4
        """
        aero_centroid = np.array([[0., 0., 0.], [4., 0., 0.]])
        aero_area = np.array([1., 3.])

        # nodes average to (3,0,0) -> distance 0 -> PASS
        log, msgs = _capture_log()
        distance, size = _check_lcp_group_location(
            log, 'fin', aero_centroid, aero_area,
            np.array([[3., 1., 0.], [3., -1., 0.]]), 0.25)
        self.assertEqual(distance, 0.)
        self.assertEqual(size, 4.)
        text = '\n'.join(msgs)
        self.assertIn('PASS', text)
        self.assertNotIn('FAIL', text)

        # node on the other side of the body: distance 2 -> 0.5 sizes -> FAIL
        log, msgs = _capture_log()
        distance, size = _check_lcp_group_location(
            log, 'fin', aero_centroid, aero_area,
            np.array([[3., -2., 0.]]), 0.25)
        self.assertAlmostEqual(distance, 2.)
        self.assertEqual(size, 4.)
        text = '\n'.join(msgs)
        self.assertIn("LCP group='fin' location", text)
        self.assertIn('FAIL: the nodes are 0.5 group sizes', text)

        # same nodes with a looser tolerance -> PASS
        log, msgs = _capture_log()
        _check_lcp_group_location(
            log, 'fin', aero_centroid, aero_area,
            np.array([[3., -2., 0.]]), 0.6)
        self.assertNotIn('FAIL', '\n'.join(msgs))

    def test_lcp_check_group_location_single_panel(self):
        """one panel has no bbox; the size falls back to sqrt(area)"""
        log, msgs = _capture_log()
        distance, size = _check_lcp_group_location(
            log, 'tab', np.array([[1., 2., 3.]]), np.array([4.]),
            np.array([[1., 2., 4.]]), 0.25)
        self.assertAlmostEqual(distance, 1.)
        self.assertEqual(size, 2.)
        self.assertIn('FAIL', '\n'.join(msgs))

    def test_lcp_groups_location_warns(self):
        """
        pressure_map_to_lcp_model logs the check per group; a group whose
        nodes are on the far side of the model is flagged, the loads are
        still mapped
        """
        log = SimpleLogger(level='warning')
        # wing (region 1, flat tri) node is near; tail (region 2) node is far
        structure_model = _lcp_structure(log, [(7/6, 5/6, 1/2), (17/6, -20., 5/6)])
        groups = {
            'wing': {'regions': [1], 'nodes': [1]},
            'tail': {'regions': [2], 'nodes': [2]},
        }
        out_model, bals, text = self._map_groups(
            structure_model, _two_region_cart3d(log), groups)
        self.assertIn("LCP group='wing' location", text)
        self.assertIn("LCP group='tail' location", text)
        wing_line = text.split("LCP group='wing' location")[1].split('LCP group=')[0]
        tail_line = text.split("LCP group='tail' location")[1].split('LCP group=')[0]
        self.assertIn('PASS', wing_line)
        self.assertIn('FAIL', tail_line)
        forces, unused_moments = self._cards(out_model)
        self.assertEqual(sorted(forces), [1, 2])
        for bal in bals.values():
            np.testing.assert_allclose(bal['struct_f'], bal['aero_f'], **LOG_TOL)
            np.testing.assert_allclose(bal['struct_m'], bal['aero_m'], **LOG_TOL)


class TestLogAeroRegions(unittest.TestCase):
    """
    log_aero_force_moment / get_aero_model with regions_to_include/remove
    (_two_region_cart3d at qinf=2: region 1 F=(0,0,2); region 2 F=(-2,0,2))
    """
    def test_log_aero_force_moment_regions(self):
        log = SimpleLogger(level='error')
        model = _two_region_cart3d(log)
        force, unused_moment = log_aero_force_moment(model, 'cart3d', qinf=PANEL_QINF)
        np.testing.assert_allclose(force, [-2., 0., 4.], atol=1e-12)

        force, unused_moment = log_aero_force_moment(
            model, 'cart3d', qinf=PANEL_QINF, regions_to_remove=[2])
        np.testing.assert_allclose(force, [0., 0., 2.], atol=1e-12)

        force, unused_moment = log_aero_force_moment(
            model, 'cart3d', qinf=PANEL_QINF, regions_to_include=[2])
        np.testing.assert_allclose(force, [-2., 0., 2.], atol=1e-12)

        # include and remove the same region -> nothing left
        log, msgs = _capture_log()
        model = _two_region_cart3d(log)
        out = log_aero_force_moment(
            model, 'cart3d', qinf=PANEL_QINF,
            regions_to_include=[2], regions_to_remove=[2])
        self.assertIsNone(out)
        self.assertIn('no aero elements left', '\n'.join(msgs))

    def test_get_aero_model_logs_filtered_sum(self):
        """get_aero_model's logged sum skips regions_to_remove; the model isn't filtered"""
        log, msgs = _capture_log()
        model = _two_region_cart3d(log)
        model_out, unused_variables = get_aero_model(
            model, 'cart3d', qinf=PANEL_QINF, regions_to_remove=[2], log=log)
        text = '\n'.join(msgs)
        self.assertIn('aero force/moment sum (1 of 2 aero elements; '
                      'regions_to_include=None regions_to_remove=[2])', text)
        self.assertIn('F = [0, 0, 2]', text)
        self.assertEqual(len(model_out.elements), 2)


def _two_region_cart3d(log: SimpleLogger, keep_region: int=0) -> Cart3D:
    """
    _two_tri_cart3d with the flat tri in region 1 and the tilted tri in region 2

    keep_region : int
        0 -> both tris; 1/2 -> only that tri
    """
    model = _two_tri_cart3d(log)
    model.regions = np.array([1, 2], dtype='int32')
    if keep_region:
        ikeep = np.where(model.regions == keep_region)[0]
        model.elements = model.elements[ikeep, :]
        model.regions = model.regions[ikeep]
        model.loads = {'Cp': model.loads['Cp'][ikeep]}
    return model


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
