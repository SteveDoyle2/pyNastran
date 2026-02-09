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
from pyNastran.bdf.mesh_utils.aero.export_caero_mesh import export_caero_mesh
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
                      tin, tout, nrows, ncols,
                      GCj, GCi,
                      Real, Complex=None, comment='wkk')

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
        tin = tout = 'float32'
        nrows = 1
        GCj = [101]
        reals = [np.radians(5.)]
        model.add_dmi_w2gj(tin, tout, nrows, reals, GCj=GCj)
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
                '--pid', 'caero', '--aerobox']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', path / 'cpmopt.paero.bdf',
                '--pid', 'caero']
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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
