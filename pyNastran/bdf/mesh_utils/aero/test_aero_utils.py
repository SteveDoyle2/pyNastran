"""various mesh_utils tests"""
import os
import copy
import unittest
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.aero.deform_aero_spline import (
    deform_aero_spline, deform_aero_spline_from_files)
from pyNastran.bdf.mesh_utils.aero.export_caero_mesh import export_caero_mesh, get_skj, rodriguez_rotate
from pyNastran.bdf.mesh_utils.aero.map_aero_model import map_aero_model
from pyNastran.bdf.mesh_utils.aero.map_pressure_to_caero import map_caero

from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.cmd_line.bdf_cmd_line import cmd_line

from pyNastran.bdf.mesh_utils.normals import get_normals_at_nodes, get_normals_at_elements
import pyNastran.bdf.mesh_utils.add_remove_mesh


TEST_DIR = (Path(__file__) / '..').resolve()
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()
BWB_PATH = MODEL_PATH / 'bwb'

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)
DIRNAME = Path(os.path.dirname(__file__))


class TestMeshUtilsAero(unittest.TestCase):
    def test_bwb_caero_map(self):
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        log = SimpleLogger(level='warning')
        map_caero(bdf_filename, log=log)

    def test_map_aero_model(self):
        """tests ``map_aero_model``"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        bdf_filename_out = BWB_PATH / 'bwb_saero_mapped.bdf'
        log = SimpleLogger(level='warning')
        model_old = read_bdf(bdf_filename, log=log)
        get_normals_at_nodes(model_old)
        get_normals_at_elements(model_old)

        copy.deepcopy(model_old)
        model_new = read_bdf(bdf_filename, xref=False, log=log)

        # make the results garbage
        for set1 in model_new.sets.values():
            set1.ids = [-1]

        map_aero_model(model_old, model_new, bdf_filename_out,
                       remove_new_aero_cards=True)

    def test_deform_aero_spline(self):
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        structure_model = read_bdf(bdf_filename, log=log, debug=False)
        displacement_aero = deform_aero_spline(
            structure_model,
            aero_model=None,
            nids=None,
            xyz_cid0=None,
            displacement0=None)

    def test_deform_aero_spline_from_files(self):
        dirname = MODEL_PATH / 'bwb'
        bdf_filename = dirname / 'bwb_saero.bdf'
        op2_filename = dirname / 'bwb_saero.op2'
        bdf_filename_out = dirname / 'bwb_saero_spline.bdf'
        op2_filename_out = dirname / 'bwb_saero_spline.op2'
        if op2_filename.exists():
            deform_aero_spline_from_files(
                bdf_filename, op2_filename,
                bdf_filename_out=bdf_filename_out,
                op2_filename_out=op2_filename_out,
            )

    def test_export_caero_mesh_caero1_wkk(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=None)
        model.bdf_filename = 'test'
        p1 = [0., 0., 0.]
        p4 = [0., 1., 0.]
        eid = 10
        pid = 11
        igroup = 1
        model.add_caero1(eid, pid, igroup,
                   p1=p1, x12=1.0,
                   p4=p4, x43=1.0,
                   nspan=2,
                   nchord=2)
        model.add_paero1(pid)
        model.add_aero(velocity=0., cref=1.0, rho_ref=1.0)
        name = 'WKK'
        form = 'rectangular'
        tin = tout = 1
        nrows = 8
        ncols = 1
        GCj = [1, 1, 1, 1, 1, 1, 1, 1]
        GCi = [1, 2, 3, 4, 5, 6, 7, 8]
        Real = [1., 2., 3., 4., 5., 6., 7., 8.]
        model.add_dmi(name, form,
                      tin, nrows, ncols,
                      GCj, GCi,
                      Real, Complex=None, tout=tout, comment='wkk')

    def test_export_caero_mesh_caero1_w2gj_dmik(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.bdf_filename = 'caero1_w2gj'
        p1 = [0., 0., 0.]
        p4 = [0., 1., 0.]
        eid = 10
        pid = 11
        igroup = 1
        npanels = 0
        caero = model.add_caero1(
            eid, pid, igroup,
            p1=p1, x12=1.0,
            p4=p4, x43=1.0,
            nspan=2, nchord=2)
        npointsi, npanelsi = caero.get_panel_npoints_nelements()
        assert npanelsi == 4, npanelsi
        npanels += npanelsi

        model.add_paero1(pid)
        model.add_aero(velocity=0., cref=1.0, rho_ref=1.0)
        # name = 'W2GJ'
        name = 'FA2J'
        form = 'column'
        #tin = tout = 1
        # nrows = 8
        # ncols = 1
        # GCj = [1, 1, 1, 1, 1, 1, 1, 1]
        # GCi = [1, 2, 3, 4, 5, 6, 7, 8]
        Real = [1., 2., 3., 4., 5., 6., 7., 8.]
        real_array = np.ones((len(Real), 1))
        dmij = model.add_dense_dmijk(
            'DMIJ', name, real_array, form,
            validate=True)
        model.dmij[name] = dmij
        with self.assertRaises(IndexError):
            export_caero_mesh(model, is_aerobox_model=True)
        model.log.error('IndexError :(')
        # save_load_deck(model, run_remove_unused=False,
        #                run_mirror=False)

    def test_export_caero_mesh_caero1_wkk2(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.bdf_filename = 'test'
        p1 = [0., 0., 0.]
        p4 = [0., 1., 0.]
        eid = 10
        pid = 11
        igroup = 1
        npanels = 0
        caero = model.add_caero1(
            eid, pid, igroup,
            p1=p1, x12=1.0,
            p4=p4, x43=1.0,
            nspan=2, nchord=2)
        npointsi, npanelsi = caero.get_panel_npoints_nelements()
        assert npanelsi == 4, npanelsi
        npanels += npanelsi

        model.add_paero1(pid)
        model.add_aero(velocity=0., cref=1.0, rho_ref=1.0)
        name = 'WKK'
        form = 'diagonal'
        #tin = tout = 1
        # nrows = 8
        # ncols = 1
        # GCj = [1, 1, 1, 1, 1, 1, 1, 1]
        # GCi = [1, 2, 3, 4, 5, 6, 7, 8]
        Real = [1., 2., 3., 4., 5., 6., 7., 8.]
        real_array = np.ones((len(Real), 1))
        model.add_dense_dmi(name, real_array, form, validate=True)
        export_caero_mesh(model, is_aerobox_model=True, )
        save_load_deck(model, run_remove_unused=False,
                       run_mirror=False)

    def test_export_caero_mesh_caero5_wtfact(self):
        """tests multiple ``bdf`` tools"""
        path = MODEL_PATH / 'aero'
        bdf_filename = str(path / 'ha145z.bdf')
        #bdf export_caero_mesh IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--aerobox] [--pid PID]\n'
        #is_aerobox_model : bool; default=True
        #    True : write the aerobox as CQUAD4s
        #    False : write the macro elements as CQUAD4s
        #pid_method : str; default='aesurf'
        #    'aesurf' : write the referenced AESURF as the property ID
        #               main structure will be pid=1
        #    'caero' : write the CAERO1 as the property id
        #    'paero' : write the PAERO1 as the property id
        log = SimpleLogger(level='warning')
        model = read_bdf(bdf_filename, log=log)
        skj = get_skj(model, percent_location=25)
        # print(skj.tolist())
        # assert skj
        tin = tout = 'float32'
        nrows = 1
        GCj = [101]
        reals = [np.radians(5.)]
        model.add_dmi_w2gj(tin, nrows, reals, tout=tout, GCj=GCj)
        #print(model.caeros)
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'ha145z.aesurf_aeroboxes.bdf',
                '--pid', 'aesurf', '--aerobox']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'ha145z.aesurf.bdf',
                '--pid', 'aesurf']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'ha145z.caero.bdf', '--pid', 'caero']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'ha145z.paero.bdf',
                '--pid', 'paero']
        cmd_line(argv=argv, log=log, quiet=True)

    def test_export_caero_mesh_w2gj(self):
        path = MODEL_PATH / 'aero'
        bdf_filename = str(path / 'cpmopt.bdf')
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'cpmopt.paero.bdf',
                '--pid', 'caero', '--aerobox', '--skip_zero_check']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'cpmopt.paero.bdf',
                '--pid', 'caero', '--skip_zero_check']
        cmd_line(argv=argv, quiet=True)

    def test_export_caero_mesh(self):
        """tests multiple ``bdf`` tools"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        argv = ['bdf', 'export_caero_mesh', str(bdf_filename), '-o', 'caero_no_sub.bdf']
        with self.assertRaises(SystemExit):
            cmd_line(argv=argv[:1], quiet=True)
        with self.assertRaises(SystemExit):
            cmd_line(argv=argv[:2], quiet=True)
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_aesurf.bdf',
                '--aerobox', '--pid', 'aesurf']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_caero.bdf',
                '--aerobox', '--pid', 'caero']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_paero.bdf',
                '--aerobox', '--pid', 'paero']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero.bdf',
                '--aerobox']
        cmd_line(argv=argv, quiet=True)

        #bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]
        argv = ['bdf', 'mirror', 'caero.bdf', '-o', 'caero2.bdf',
                '--plane', 'xz', '--tol', '1e-5']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'equivalence', 'caero2.bdf', '0.001', '-o', 'caero3.bdf']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'merge', 'caero2.bdf', 'caero2.bdf', '-o', 'caero3_merged.bdf']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'renumber', 'caero3.bdf', 'caero4.bdf', '--size', '8']
        cmd_line(argv=argv, quiet=True)

        #bdf transform IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--shift XYZ]
        argv = ['bdf', 'transform', 'caero4.bdf', '-o', 'caero5.bdf', '--shift', '0,0,20.']
        cmd_line(argv=argv, quiet=True)

        #'  bdf convert IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
        argv = ['bdf', 'convert', 'caero5.bdf',
                '-o', 'caero6.bdf',
                '--in_units', 'in,lbm', '--out_units', 'ft,lbm']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'scale', 'caero6.bdf',
                #'-o', 'caero6.bdf',
                '--length', '0.5', '--time', '1.', '--mass', str(0.5**3.)]
        cmd_line(argv=argv, quiet=True)

        os.remove('caero.bdf')
        os.remove('caero2.bdf')
        os.remove('caero3.bdf')
        os.remove('caero3_merged.bdf')
        os.remove('caero4.bdf')
        os.remove('caero5.bdf')
        os.remove('caero6.bdf')
        #os.remove('caero5.scaled.bdf')
        os.remove('caero6.scaled.bdf')
        os.remove('caero_aesurf.bdf')
        os.remove('caero_caero.bdf')
        os.remove('caero_paero.bdf')
        os.remove('caero_no_sub.bdf')


class TestRodriguezRotate(unittest.TestCase):
    """Tests for rodriguez_rotate: Rodrigues' rotation about a coord axis through coord.origin."""

    def _make_coord(self, origin, zaxis, xzplane):
        """Create a CORD2R and cross-reference it."""
        model = BDF(debug=False)
        model.add_cord2r(1, origin=origin, zaxis=zaxis, xzplane=xzplane)
        model.cross_reference()
        return model.coords[1]

    def test_rotate_90_about_z_at_origin(self):
        """Rotate (1,0,0) by 90 deg about z-axis at (0,0,0) -> (0,1,0)."""
        coord = self._make_coord([0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        xyz = np.array([[1., 0., 0.]])
        # iaxis=2 -> k-axis = z
        result = rodriguez_rotate(xyz, np.pi / 2, coord, iaxis=2)
        np.testing.assert_allclose(result, [[0., 1., 0.]], atol=1e-14)

    def test_rotate_180_about_y_at_origin(self):
        """Rotate (1,0,0) by 180 deg about y-axis at (0,0,0) -> (-1,0,0)."""
        coord = self._make_coord([0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        xyz = np.array([[1., 0., 0.]])
        # iaxis=1 -> j-axis = y
        result = rodriguez_rotate(xyz, np.pi, coord, iaxis=1)
        np.testing.assert_allclose(result, [[-1., 0., 0.]], atol=1e-14)

    def test_rotate_about_offset_origin(self):
        """Rotate about z-axis at (1,0,0): point (2,0,0) rotated 90 deg -> (1,1,0).

        The point is 1 unit from the axis in the x-direction. After 90 deg
        rotation about z through (1,0,0), it should move to (1,1,0).
        """
        coord = self._make_coord([1., 0., 0.], [1., 0., 1.], [2., 0., 0.])
        xyz = np.array([[2., 0., 0.]])
        result = rodriguez_rotate(xyz, np.pi / 2, coord, iaxis=2)
        np.testing.assert_allclose(result, [[1., 1., 0.]], atol=1e-14)

    def test_rotate_point_on_axis_unchanged(self):
        """Point on the rotation axis should not move."""
        coord = self._make_coord([1., 2., 3.], [1., 2., 4.], [2., 2., 3.])
        xyz = np.array([[1., 2., 5.]])  # on z-axis through origin
        result = rodriguez_rotate(xyz, np.radians(45), coord, iaxis=2)
        np.testing.assert_allclose(result, [[1., 2., 5.]], atol=1e-14)

    def test_rotate_multiple_points(self):
        """Vectorized rotation of multiple points about x-axis at origin."""
        coord = self._make_coord([0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        xyz = np.array([
            [0., 1., 0.],
            [0., 0., 1.],
            [5., 1., 0.],
        ])
        # 90 deg about x: (0,1,0)->(0,0,1), (0,0,1)->(0,-1,0), (5,1,0)->(5,0,1)
        result = rodriguez_rotate(xyz, np.pi / 2, coord, iaxis=0)
        expected = np.array([
            [0., 0., 1.],
            [0., -1., 0.],
            [5., 0., 1.],
        ])
        np.testing.assert_allclose(result, expected, atol=1e-14)

    def test_zero_angle_identity(self):
        """Zero rotation angle should return the original points."""
        coord = self._make_coord([3., 4., 5.], [3., 4., 6.], [4., 4., 5.])
        xyz = np.array([[10., 20., 30.], [1., 2., 3.]])
        result = rodriguez_rotate(xyz, 0.0, coord, iaxis=1)
        np.testing.assert_allclose(result, xyz, atol=1e-14)

    def test_full_rotation_returns_to_start(self):
        """360 deg rotation should return to the starting point."""
        coord = self._make_coord([1., 1., 1.], [1., 1., 2.], [2., 1., 1.])
        xyz = np.array([[3., 4., 5.], [-1., 2., 7.]])
        result = rodriguez_rotate(xyz, 2 * np.pi, coord, iaxis=2)
        np.testing.assert_allclose(result, xyz, atol=1e-13)

    def test_invalid_iaxis_raises(self):
        """Invalid iaxis should raise RuntimeError."""
        coord = self._make_coord([0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        xyz = np.array([[1., 0., 0.]])
        with self.assertRaises(RuntimeError):
            rodriguez_rotate(xyz, 0.5, coord, iaxis=3)


class TestGetSKJAgainstNastran(unittest.TestCase):
    """Validates get_skj() against Nastran-exported SKJ from DMAP ALTER.

    The OP4 file was generated by running bwb_saero_trim_export.bdf with
    a PFAERO ALTER that exports AJJ, SKJ, D1JK after AIC computation.
    See models/bwb/bwb_saero_trim_export.bdf for the ALTER.

    SKJ is (2*N_panels, N_panels) where:
      - Even rows (2j): panel area (force per unit Cp)
      - Odd rows (2j+1): panel area * moment_arm (moment per unit Cp)
      - moment_arm = chord * 0.25 (quarter-chord)

    Tolerances:
      - Force rows (areas): rtol=1e-10 (identical geometry)
      - Moment rows (area*arm): rtol=1e-10 (identical geometry)
    """

    OP4_FILE = BWB_PATH / 'bwb_aestatrs_matrices.op4'

    @unittest.skipIf(
        not (BWB_PATH / 'bwb_aestatrs_matrices.op4').exists(),
        'Nastran OP4 export not available')
    def test_skj_force_rows_match_nastran(self):
        """Force rows (areas) must match Nastran's SKJ to machine precision."""
        """Moment rows (area*arm) must match Nastran's SKJ."""
        from pyNastran.op4.op4 import read_op4

        log = SimpleLogger(level='warning')
        model = read_bdf(BWB_PATH / 'bwb_saero.bdf', xref=True, log=log)
        skj_ours = get_skj(model, percent_location=25)

        matrices = read_op4(str(self.OP4_FILE), debug=False)
        skj_nastran = np.asarray(matrices['SKJ'].data)

        assert skj_ours.shape == skj_nastran.shape, (
            f'Shape mismatch: ours={skj_ours.shape} nastran={skj_nastran.shape}')

        # Force rows (even indices) = panel areas
        np.testing.assert_allclose(
            skj_ours[0::2, :], skj_nastran[0::2, :],
            rtol=1e-10, atol=1e-14,
            err_msg='SKJ force rows (panel areas) differ from Nastran')

        # Moment rows (odd indices) = area * arm
        np.testing.assert_allclose(
            skj_ours[1::2, :], skj_nastran[1::2, :],
            rtol=1e-10, atol=1e-14,
            err_msg='SKJ moment rows (area*arm) differ from Nastran')

    def test_skj_bwb_shape_and_structure(self):
        """SKJ for the BWB model has 225 panels with banded structure."""
        log = SimpleLogger(level='warning')
        model = read_bdf(BWB_PATH / 'bwb_saero.bdf', xref=True, log=log)
        skj = get_skj(model, percent_location=25)

        assert skj.shape == (450, 225)

        # Exactly 2 nonzeros per column
        for j in range(225):
            nz = np.count_nonzero(skj[:, j])
            assert nz == 2, f'Column {j} has {nz} nonzeros, expected 2'

        # All force rows (areas) positive
        areas = np.array([skj[2*j, j] for j in range(225)])
        assert np.all(areas > 0)

        # All moment arms non-negative (chord * 0.25 >= 0)
        arms = np.array([skj[2*j+1, j] / skj[2*j, j] for j in range(225)])
        assert np.all(arms >= 0)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
