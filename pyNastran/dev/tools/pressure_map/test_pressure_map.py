from pathlib import Path
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
    pressure_map_to_structure_model)
from pyNastran.dev.tools.pressure_map.setup_aero import (
    get_aero_model, get_aero_pressure_centroid)

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_DIR = PKG_PATH / '..' / 'models'
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
        aero_model.loads['Cp'] = cp

        pressure_map(
            aero_model,
            bdf_model,
            # eids_structure=np.array([]),
            # eid_csv_filename='',
            eid_load_id=-1,
            aero_format=aero_format,
            map_type='pressure',
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
        vrt_filename = MODEL_DIR / 'bwb' / 'bwb-saero.vrt'
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
        log.info('running 1st prssure map')
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

        log.info('running 2nd prssure map')
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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
