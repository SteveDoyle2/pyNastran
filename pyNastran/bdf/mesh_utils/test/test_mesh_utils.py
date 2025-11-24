"""various mesh_utils tests"""
import os
import copy
import unittest
from pathlib import Path
from io import StringIO

from docopt import DocoptExit
import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.map_pressure_to_caero import map_caero

from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.split_elements import split_line_elements
from pyNastran.bdf.mesh_utils.pierce_shells import (
    pierce_shell_model) #, quad_intersection, triangle_intersection)
from pyNastran.bdf.mesh_utils.mirror_mesh import (
    write_bdf_symmetric, bdf_mirror, bdf_mirror_plane)
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  #mass_properties_breakdown
from pyNastran.bdf.mesh_utils.make_half_model import make_half_model
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.bdf.mesh_utils.utils import cmd_line, CMD_MAPS
from pyNastran.bdf.mesh_utils.find_closest_nodes import find_closest_nodes
from pyNastran.bdf.mesh_utils.find_coplanar_elements import find_coplanar_triangles
from pyNastran.bdf.mesh_utils.force_to_pressure import force_to_pressure
from pyNastran.bdf.mesh_utils.free_edges import free_edges, non_paired_edges
from pyNastran.bdf.mesh_utils.get_oml import get_oml_eids
from pyNastran.bdf.mesh_utils.breakdowns import (
    get_mass_breakdown, get_area_breakdown, get_length_breakdown,
    get_volume_breakdown, get_thickness_breakdown,
    get_material_mass_breakdown_table, get_property_mass_breakdown_table)
from pyNastran.bdf.mesh_utils.map_aero_model import map_aero_model

from pyNastran.bdf.mesh_utils.mesh import create_structured_cquad4s, create_structured_chexas
from pyNastran.bdf.mesh_utils.cmd_line.bdf_merge import cmd_line_merge
from pyNastran.bdf.mesh_utils.dvxrel import get_dvprel_ndarrays
from pyNastran.bdf.mesh_utils.rbe_tools import (
    merge_rbe2, rbe3_to_rbe2, rbe2_to_rbe3)

import pyNastran.bdf.mesh_utils.add_remove_mesh
import pyNastran.bdf.mesh_utils.bdf_remove_comments
import pyNastran.bdf.mesh_utils.cleanup_model
import pyNastran.bdf.mesh_utils.normals
import pyNastran.bdf.mesh_utils.add_remove_mesh


TEST_DIR = (Path(__file__) / '..').resolve()
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()
BWB_PATH = MODEL_PATH / 'bwb'

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)
DIRNAME = Path(os.path.dirname(__file__))


class TestRbeTools(unittest.TestCase):

    def test_rbe2_to_rbe3_to_rbe2(self):
        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.add_grid(5, [0., 0., 0.])

        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_conm2(1, 10, 1.0)
        model.add_conm2(2, 11, 1.0)
        model.add_conm2(3, 12, 1.0)
        model.add_rbe2(1, 10, '123456', [1, 2])
        model.add_rbe2(2, 11, '123456', [2, 3])
        model.add_rbe2(3, 12, '123456', [4, 5])
        bdf_filename1 = DIRNAME / 'rbe3_test_1.bdf'
        bdf_filename2 = DIRNAME / 'rbe3_test_2.bdf'
        bdf_filename3 = DIRNAME / 'rbe3_test_3.bdf'
        model.write_bdf(bdf_filename1)
        rbe_eids_to_fix = list(model.rigid_elements)
        rbe2_to_rbe3(model, rbe_eids_to_fix)
        rbe3_to_rbe2(model, rbe_eids_to_fix)
        args = ['bdf', 'rbe3_to_rbe2', str(bdf_filename1), '-o', str(bdf_filename2)]
        cmd_line(args, quiet=True)
        args = ['bdf', 'rbe3_to_rbe2', str(bdf_filename2), '-o', str(bdf_filename3)]
        cmd_line(args, quiet=True)
        for fname in [bdf_filename1, bdf_filename2, bdf_filename3]:
            os.remove(fname)

    def test_merge_rbe2(self):
        model = BDF()
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.add_grid(5, [0., 0., 0.])

        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_conm2(1, 10, 1.0)
        model.add_conm2(2, 11, 1.0)
        model.add_conm2(3, 12, 1.0)
        model.add_rbe2(1, 10, '123456', [1, 2])
        model.add_rbe2(2, 11, '123456', [2, 3])
        model.add_rbe2(3, 12, '123456', [4, 5])
        rbe_eids_to_fix = list(model.rigid_elements)
        assert len(model.nodes) == 8, len(model.nodes)
        assert len(model.rigid_elements) == 3, len(model.rigid_elements)
        assert len(model.masses) == 3, len(model.masses)

        merge_rbe2(model, rbe_eids_to_fix)
        assert len(model.nodes) == 8, len(model.nodes)
        assert len(model.rigid_elements) == 2, len(model.rigid_elements)
        assert len(model.masses) == 2, len(model.masses)


class TestMeshUtilsAero(unittest.TestCase):
    def test_bwb_caero_map(self):
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        map_caero(bdf_filename)

    def test_map_aero_model(self):
        """tests ``map_aero_model``"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        bdf_filename_out = BWB_PATH / 'bwb_saero_mapped.bdf'
        log = SimpleLogger(level='warning')
        model_old = read_bdf(bdf_filename, log=log)
        copy.deepcopy(model_old)
        model_new = read_bdf(bdf_filename, xref=False, log=log)

        # make the results garbage
        for set1 in model_new.sets.values():
            set1.ids = [-1]

        map_aero_model(model_old, model_new, bdf_filename_out,
                       remove_new_aero_cards=True)

    def test_export_caero_mesh_caero1_wkk(self):
        model = BDF(debug=None)
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
        model = BDF(debug=None)
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
        model = read_bdf(bdf_filename, debug=False)
        tin = tout = 'float32'
        nrows = 1
        GCj = [101]
        Real = [np.radians(5.)]
        model.add_dmi_w2gj(tin, tout, nrows, GCj, Real)
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
        cmd_line(argv=argv, quiet=True)

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
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_no_sub.bdf']
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


class TestMeshUtilsCmdLine(unittest.TestCase):
    def test_solid_dof(self):
        bdf_filename = TEST_DIR / 'solid_dof.bdf'
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])
        model.add_grid(4, [0., 0., 1.])
        model.add_ctetra(10, 1, [1, 2, 3, 4])
        model.add_psolid(1, 1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        model.write_bdf(bdf_filename)

        args = ['bdf', 'solid_dof', str(bdf_filename)]
        model, out_nids = cmd_line(args, quiet=True)
        assert len(out_nids) == 4, out_nids

        model.add_grid(5, [0., 0., 1.])
        rbe2 = model.add_rbe2(10, 4, '123456', 5)
        model.write_bdf(bdf_filename)
        model, out_nids = cmd_line(args, quiet=True)
        assert len(out_nids) == 3, out_nids

        model.add_conm2(12, 2, 1.0)
        model.write_bdf(bdf_filename)
        model, out_nids = cmd_line(args, quiet=True)
        assert len(out_nids) == 2, out_nids

        model.add_conrod(100, 1, [1, 2])
        model.write_bdf(bdf_filename)
        model, out_nids = cmd_line(args, quiet=True)
        assert len(out_nids) == 1, out_nids

        os.remove(bdf_filename)

    def test_inclzip_bwb(self):
        """tests ``inclzip``"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        args = ['bdf', 'inclzip', str(bdf_filename)]
        cmd_line(args, quiet=True)

    def test_stats_bwb(self):
        """tests ``stats``"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        args = ['bdf', 'stats', str(bdf_filename)]
        cmd_line(args, quiet=True)

    # def test_free_edges_bwb(self):
    #     """tests ``free_edges``"""
    #     bdf_filename = BWB_PATH / 'bwb_saero.bdf'
    #     args = ['bdf', 'free_edges', str(bdf_filename)]
    #     cmd_line(args, quiet=True)

    def test_bdf_stats(self):
        """tests ```bdf stats```"""
        bdf_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.bdf'
        args = ['bdf', 'stats', str(bdf_filename)]
        cmd_line(args, quiet=True)

    def test_bdf_diff(self):
        """tests ```bdf diff```"""
        bdf_filename1 = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.bdf'
        bdf_filename2 = MODEL_PATH / 'sol_101_elements' / 'mode_solid_shell_bar.bdf'
        args = ['bdf', 'diff', str(bdf_filename1), str(bdf_filename2)]
        cmd_line(args, quiet=True)

    def test_free_edges(self):
        """Finds the free_edges

        4-----3---5
        |   / |
        |  /  |
        | /   |
        1-----2
        """
        log = SimpleLogger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_grid(5, [1., 1., 1.])  # 1,3,5
        model.add_ctria3(1, 1, [1, 2, 3])

        edges1 = free_edges(model, eids=None, maps=None)
        edges2 = non_paired_edges(model, eids=None, maps=None)
        assert edges1 == [(1, 2), (2, 3), (1, 3)]
        assert edges2 == [(1, 2), (2, 3), (1, 3)]

        edges1 = free_edges(model, eids=[1], maps=None)
        edges2 = non_paired_edges(model, eids=[1], maps=None)
        assert edges1 == [(1, 2), (2, 3), (1, 3)]
        assert edges2 == [(1, 2), (2, 3), (1, 3)]

        model.add_ctria3(2, 1, [1, 3, 4])
        edges1 = free_edges(model, eids=None, maps=None)
        edges2 = non_paired_edges(model, eids=None, maps=None)
        assert edges1 == [(1, 2), (2, 3), (3, 4), (1, 4)], edges1
        assert edges2 == [(1, 2), (2, 3), (3, 4), (1, 4)], edges2

        # edge (1, 3) is associated with 3 elements, so it's not free,
        # but it's not paired (an edge with only 2 elements)
        model.add_ctria3(3, 1, [1, 3, 5])
        edges1 = free_edges(model, eids=None, maps=None)
        edges2 = non_paired_edges(model, eids=None, maps=None)
        assert edges1 == [(1, 2), (2, 3),         (3, 4), (1, 4), (3, 5), (1, 5)], edges1
        assert edges2 == [(1, 2), (2, 3), (1, 3), (3, 4), (1, 4), (3, 5), (1, 5)], edges2

        bdf_filename = TEST_DIR / 'test_free_edges.bdf'
        model.write_bdf(bdf_filename)
        args = ['bdf', 'collapse_quads', str(bdf_filename), '--punch', '--size', '16']
        cmd_line(args, quiet=True)

        bdf_filename = TEST_DIR / 'test_free_edges_quad.bdf'
        bdf_filename_out = TEST_DIR / 'test_free_edges_quad.bdf'
        model.add_grid(100, [0., 0., 0.])
        model.add_cquad4(100, 1, [1, 2, 3, 1])
        model.write_bdf(bdf_filename)
        args2 = ['bdf', 'collapse_quads', str(bdf_filename), '--punch',
                '-o', str(bdf_filename_out)]
        cmd_line(args2, quiet=True)

        #bdf flip_shell_normals IN_BDF_FILENAME[-o OUT_BDF_FILENAME] [--punch][--zero_zoffset]
        args3 = ['bdf', 'flip_shell_normals', str(bdf_filename_out), '--punch']
        cmd_line(args3, quiet=True)

        args4 = ['bdf', 'remove_unused', str(bdf_filename_out), '--punch']
        cmd_line(args4, quiet=True)

    def test_free_faces(self):
        """CTETRA10"""
        #bdf free_faces [-d | -l] [-f] [--encoding ENCODE] BDF_FILENAME SKIN_FILENAME
        #with self.assertRaises(SystemExit):
            #cmd_line(argv=['bdf', 'free_faces'])
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        #log = SimpleLogger(level='info', encoding='utf-8')
        skin_filename = DIRNAME / 'skin.bdf'
        cmd_line(argv=['bdf', 'free_faces', str(bdf_filename), str(skin_filename)],
                 quiet=True)
        os.remove(skin_filename)

    def test_exit(self):
        """tests totally failing to run"""
        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf'], quiet=True)

        for cmd in CMD_MAPS:
            with self.assertRaises(SystemExit):
                cmd_line(argv=['bdf', cmd], quiet=True)


class TestMeshUtils(unittest.TestCase):
    """various mesh_utils tests"""
    def test_dvxrel(self):
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])
        model.add_pshell(1, 1, 0.1)
        model.add_pshell(2, 1, 0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)

        oid = 1
        prop_type = 'PSHELL'
        pid = 1
        pid2 = 2
        pname = 'T'
        desvar_ids = [1]
        coeffs = [1.0]
        model.add_dvprel1(oid, prop_type, pid, pname,
                          desvar_ids, coeffs, p_min=0.4, p_max=2.0)
        model.add_dvprel1(oid+1, prop_type, pid2, pname,
                          desvar_ids, coeffs, p_min=0.4, p_max=2.0)
        model.add_desvar(1, 'T1', 1.0, xlb=0.5, xub=1.5)
        properties = np.array([1])
        nelements = len(model.elements)
        get_dvprel_ndarrays(model, nelements, properties)

    def test_coplanar_triangles(self):
        """tests find_coplanar_triangles

        4-----3
        |   / |
        |  /  |
        | /   |
        1-----2
        """
        log = SimpleLogger(level='warning')
        model = BDF(debug=True, log=log, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])
        _junk, coplanar_elements = find_coplanar_triangles(model, eids=None)
        assert coplanar_elements == set([]), coplanar_elements

        model.add_ctria3(2, 1, [2, 3, 1])
        _junk, coplanar_elements = find_coplanar_triangles(model, eids=None)
        assert coplanar_elements == {2}, coplanar_elements

        model.add_ctria3(3, 1, [3, 2, 1])
        _junk, coplanar_elements = find_coplanar_triangles(model, eids=None)
        assert coplanar_elements == {2, 3}, coplanar_elements

        _junk, coplanar_elements = find_coplanar_triangles(model, eids=[2, 3])
        assert coplanar_elements == {3}, coplanar_elements

        model.add_cquad4(4, 1, [1, 2, 3, 4])
        _junk, coplanar_elements = find_coplanar_triangles(model, eids=None)
        assert coplanar_elements == {2, 3}, coplanar_elements

    def test_force_to_pressures(self):
        log = SimpleLogger(level='warning')

        model = BDF(debug=False, log=log, mode='msc')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        # area = 0.5
        model.add_ctria3(1, 1, [1, 2, 3])
        model.add_pshell(1, 1, t=0.1)
        model.add_mat1(1, 3.0e7, None, 0.3)
        force_to_pressure(model, clear_model=False)

        #force = 3.0
        # pressure = 6.0
        model.add_force(1, 1, 1.0, [0., 0., 1.])
        model.add_force(1, 2, 1.0, [0., 0., 1.])
        model.add_force(1, 3, 1.0, [0., 0., 1.])
        force_to_pressure(model, clear_model=True)

    def test_structured_cquads(self):
        """tests create_structured_cquad4s"""
        pid = 42
        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [1., 1., 0.]
        p4 = [0., 1., 0.]
        model = BDF()
        nx = 10
        ny = 20
        create_structured_cquad4s(model, pid, p1, p2, p3, p4, nx, ny, nid=1, eid=1, theta_mcid=0.)

    def test_structured_chexas(self):
        """tests test_structured_chexas"""
        #1U CubeSat is 10 cm, 10 cm, 11.35 cm.
        #2U CubeSat is 10 cm, 10 cm, 22.70 cm.
        #6U CubeSat is 20 cm, 10 cm, 34.05 cm.
        model = BDF(debug=False)
        pid1 = 1
        pid2 = 2

        xmax = 10.
        ymax = 10.
        zmax = 10.
        unused_p1 = [0., 0., 0.]
        unused_p2 = [xmax, 0., 0.]
        unused_p3 = [xmax, ymax, 0.]
        unused_p4 = [0., ymax, 0.]

        unused_p5 = [0., 0., zmax]
        unused_p6 = [xmax, 0., zmax]
        unused_p7 = [xmax, ymax, zmax]
        unused_p8 = [0., ymax, zmax]

        # nnodes
        nx = 20
        ny = 20
        nz = 2

        x = np.linspace(0., xmax, nx)
        y = np.linspace(0., ymax, ny)
        z = np.linspace(0., zmax, nz)
        mid = 1
        E = 69.e9  # Pa
        G = None
        nu = 0.30
        #rho = 2_710  #kg/m^3
        uallow = 240.e6  # MPa
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0,
                       comment='aluminum')
        model.add_psolid(pid1, mid)
        model.add_psolid(pid2, mid)
        create_structured_chexas(model, pid2,
                                 x, y, z, nx, ny, nz, eid=1)


        nodes = []
        cm = 1.1
        model.cross_reference()
        for nid, node in model.nodes.items():
            x, y, z = node.xyz
            if (np.isclose(z, 0.) and
                # corner on the bottom face
                abs(x - xmax / 2.) >= (xmax / 2. - cm) and
                abs(y - ymax / 2.) >= (ymax / 2. - cm)):
                nodes.append(nid)

        constrained_elements = set(model.elements.keys())
        for eid, elem in model.elements.items():
            x, y, z = elem.Centroid()
            if (# corner elements
                abs(x - xmax / 2.) >= (ymax / 2. - cm) and
                abs(y - ymax / 2.) >= (ymax / 2. - cm)):
                elem.pid = 1
                constrained_elements.remove(eid)
        constrained_elements = list(constrained_elements)
        constrained_elements.sort()

        spc_id = 100
        components = '123'
        model.add_spc1(spc_id, components, nodes, comment='floor')

        g = 9.81  # m/s^2
        Nx = [g, 0., 0.]
        Ny = [0., g, 0.]
        Nz = [0., 0., g]
        model.add_grav(1001, -2., Nx, cid=0, mb=0, comment='-2 G Nx')
        model.add_grav(1002, 2., Nx, cid=0, mb=0, comment='+2 G Nx')
        model.add_grav(1003, -2., Ny, cid=0, mb=0, comment='-2 G Ny')
        model.add_grav(1004, 2., Ny, cid=0, mb=0, comment='+2 G Ny')
        model.add_grav(1005, -9., Nz, cid=0, mb=0, comment='-9 G Nz')
        model.add_grav(1006, 2., Nz, cid=0, mb=0, comment='+2 G Nz')

        subcases = model.create_subcases([1, 2, 3, 4, 5, 6])

        options = []
        param_type = 'STRESS-type'
        dconstr_id = 101
        dresp_id = 1001
        response_type = 'STRESS'
        property_type = 'ELEM'
        region = None
        atta = 13
        attb = None  # mode
        comment = 'VM Stress'
        for eid in constrained_elements:
            model.add_dconstr(dconstr_id, dresp_id, lid=-uallow, uid=uallow, lowfq=0.,
                              highfq=1.e20, comment=comment)
            label = f's{eid}'
            atti = [eid]
            model.add_dresp1(dresp_id, label, response_type, property_type,
                            region, atta, attb,
                            atti, validate=True,
                            comment=comment)
            dresp_id += 1
            comment = ''

        for subcase_id, subcase in subcases.items():
            load_id = 1000 + subcase_id
            subcase.add('LOAD', load_id, options, param_type)
            subcase.add('SPC', spc_id, options, param_type)
            #subcase.add('ANALYSIS', 'STATICS', options, 'KEY-type')
            subcase.add('ANALYSIS', 'STATICS', options, param_type)
            subcase.add('DESSUB', dconstr_id, options, param_type)

        bdf_filename = os.path.join(DIRNAME, 'test_structured_chexas.bdf')
        model.write_bdf(bdf_filename)

    def test_bdf_delete_bad_shells(self):
        """tests ```bdf delete_bad_shells```"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        args = ['bdf', 'delete_bad_shells', str(bdf_filename)]
        cmd_line(args, quiet=True)

    def test_bdf_remove_comments(self):
        """tests ```bdf remove_comments```"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        args = ['bdf', 'remove_comments', str(bdf_filename)]
        cmd_line(args, quiet=True)

    def test_bdf_list_conm2(self):
        """tests ```bdf list_conm2```"""
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        args = ['bdf', 'list_conm2', str(bdf_filename)]
        cmd_line(args, quiet=True)


    def test_breakdown_01(self):
        """run the various breakdowns"""
        log = SimpleLogger(level='error')
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'

        model = read_bdf(bdf_filename, log=log)
        get_mass_breakdown(model)
        get_area_breakdown(model)
        get_length_breakdown(model)
        get_volume_breakdown(model)
        get_thickness_breakdown(model)
        get_material_mass_breakdown_table(model)
        get_property_mass_breakdown_table(model)

    def test_merge_01(self):
        """merges multiple bdfs into a single deck"""
        log = SimpleLogger(level='error')
        bdf_filename1 = BWB_PATH / 'bwb_saero.bdf'
        bdf_filename2 = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        bdf_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        bdf_filename4 = os.path.join(MODEL_PATH, 'iSat', 'ISat_Dploy_Sm.dat')
        bdf_filename_out1 = MODEL_PATH / 'bwb' / 'BWBsaero_staticbar_8.out'
        bdf_filename_out2 = MODEL_PATH / 'bwb' / 'BWBsaero_static_bar_16.out'
        bdf_filename_out3 = MODEL_PATH / 'bwb' / 'BWBsaero_staticbar_isat.out'

        bdf_filenames1 = [bdf_filename1, bdf_filename2]
        bdf_filenames2 = [bdf_filename1, bdf_filename2, bdf_filename3, bdf_filename4]

        args = ['bdf', 'merge', '--debug'] + [str(pathi) for pathi in bdf_filenames1]
        cmd_line_merge(args, quiet=True)
        bdf_merge(bdf_filenames1, bdf_filename_out=bdf_filename_out1,
                  renumber=True, encoding=None, size=8, is_double=False,
                  cards_to_skip=None, log=log)
        bdf_merge(bdf_filenames1, bdf_filename_out=bdf_filename_out2,
                  renumber=False, encoding=None, size=16, is_double=False,
                  cards_to_skip=None, log=log)

        bdf_merge(bdf_filenames2, bdf_filename_out=bdf_filename_out3,
                  renumber=False, encoding=None, size=16, is_double=False,
                  cards_to_skip=None, log=log)
        read_bdf(bdf_filename_out1, log=log)
        read_bdf(bdf_filename_out2, log=log)
        read_bdf(bdf_filename_out3, log=log)

        os.remove(bdf_filename_out1)
        os.remove(bdf_filename_out2)
        os.remove(bdf_filename_out3)

    def test_export_mcids(self):
        """creates material coordinate systems"""
        log = SimpleLogger(level='error')
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        csv_filename = MODEL_PATH / 'bwb' / 'mcids.csv'
        export_mcids(bdf_filename, csv_filename,
                     export_xaxis=True, export_yaxis=True,
                     iply=9, log=log, debug=False)

        model = read_bdf(bdf_filename, xref=False, debug=False)
        model.safe_cross_reference()
        #os.remove('mcids.csv')

        argv = ['bdf', 'export_mcids', str(bdf_filename), '-o', str(csv_filename),
                '--iplies', '0,1,2,3,4,5,6,7,8,9,10', '--no_x', '--no_y']
        with self.assertRaises(DocoptExit):
            # can't define both --no_x and --no_y
            cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_mcids', str(bdf_filename), '-o', str(csv_filename),
                '--iplies', '0,1,2,3,4,5,6,7,8,9', '--no_x']
        cmd_line(argv=argv, quiet=True)

        eids = [1204, 1211]
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=True, export_yaxis=True,
                     iply=9, log=log, debug=False)
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=True, export_yaxis=False,
                     iply=9, log=log, debug=False)
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=False, export_yaxis=True,
                     iply=9, log=log, debug=False)
        with self.assertRaises(AssertionError):
            # export_xaxis and export_yaxis can't both be False
            export_mcids(model, csv_filename=None, eids=eids,
                         export_xaxis=False, export_yaxis=False,
                         iply=9)

        with self.assertRaises(RuntimeError):
            # no iply=10
            export_mcids(model, csv_filename, eids=eids,
                         export_xaxis=True, export_yaxis=True,
                         iply=10)
        #os.remove('mcids.csv')

    def test_split_cbars_by_pin_flag_1(self):
        """null pin flag test"""
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        split_cbars_by_pin_flag(bdf_filename, pin_flags_filename='pin_flags.csv',
                                bdf_filename_out='pin_flags.bdf', debug=False)
        os.remove('pin_flags.csv')
        os.remove('pin_flags.bdf')

        argv = ['bdf', 'split_cbars_by_pin_flags', bdf_filename,
                '-o', 'pin_flags.bdf', '-p', 'pin_flags.csv']
        cmd_line(argv=argv, quiet=True)
        os.remove('pin_flags.csv')
        os.remove('pin_flags.bdf')

    def test_split_cbars_by_pin_flag_2(self):
        """real pin flag test"""
        lines = [
            'SOL 101\n',
            'CEND\n',
            'SUBCASE 10\n',
            '    LOAD = 10\n',
            '    SPC = 123456\n',
            '    DISP(PLOT) = ALL\n',
            '    STRESS(PLOT) = ALL\n',
            'BEGIN BULK\n',
            'ENDDATA',
        ]
        model = BDF(debug=False)
        with self.assertRaises(NotImplementedError):
            model.read_bdf(bdf_filename=lines, validate=True, xref=True,
                           punch=False, read_includes=True,
                           encoding=None)

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [2., 0., 0.])
        model.add_grid(4, [3., 0., 0.])

        pid = 1000
        mid = 1000
        Type = 'BAR'
        dim = [1., 2.]
        model.add_pbarl(pid, mid, Type, dim)
        E = 3.0e7
        G = 3.0e6
        nu = None
        model.add_mat1(mid, E, G, nu)

        x = [0., 1., 0.]
        g0 = None
        model.add_cbar(1, pid, [1, 2], x, g0, offt='GGG', pa=0, pb=0,
                       wa=None, wb=None, comment='reaction')
        model.add_cbar(2, pid, [2, 3], x, g0, offt='GGG', pa=0, pb=456,
                       wa=None, wb=None, comment='End B')
        model.add_cbar(3, pid, [3, 4], x, g0, offt='GGG', pa=456, pb=0,
                       wa=None, wb=None, comment='End A')
        sid = 10
        node = 4
        mag = 1.
        xyz = [1., 1., 0.]
        model.add_force(sid, node, mag, xyz)
        model.add_spc1(123456, '123456', 1)
        model.validate()

        bdf_file = StringIO()
        bdf_file.writelines(lines)
        bdf_file.seek(0)
        model.read_bdf(bdf_filename=bdf_file, validate=True, xref=False,
                       punch=False, read_includes=True,
                       encoding=None)
        #model.write_bdf('spike.bdf')

        split_cbars_by_pin_flag(model, pin_flags_filename='pin_flags.csv',
                                bdf_filename_out='pin_flags.bdf', debug=False)
        os.remove('pin_flags.csv')
        os.remove('pin_flags.bdf')


    def test_split_line_elements(self):
        """tests split_line_elements"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])

        pid = 1000
        mid = 1000
        Type = 'BAR'
        dim = [1., 2.]
        model.add_pbarl(pid, mid, Type, dim)
        E = 3.0e7
        G = 3.0e6
        nu = None
        model.add_mat1(mid, E, G, nu)

        x = [0., 1., 0.]
        g0 = None
        #model.add_cbar(1, pid, 1, 2, x, g0, offt='GGG', pa=0, pb=0,
                       #wa=None, wb=None, comment='reaction')
        #model.add_cbar(2, pid, 2, 3, x, g0, offt='GGG', pa=0, pb=456,
                       #wa=None, wb=None, comment='End B')
        #model.add_cbar(3, pid, 3, 4, x, g0, offt='GGG', pa=456, pb=0,
                       #wa=None, wb=None, comment='End A')
        #eids = [1, 2, 3]

        nids = [1, 2]
        model.add_cbar(1, pid, nids, x, g0, offt='GGG', pa=456, pb=5,
                       wa=None, wb=None, comment='End A')
        model.add_cbeam(2, 2000, nids, x, g0, offt='GGG', bit=None, pa=456,
                        pb=5, wa=None, wb=None, sa=0,
                        sb=0, comment='')
        A = 42.
        model.add_conrod(3, mid, nids, A)
        model.add_prod(4000, mid, A)
        model.add_crod(4, 4000, nids)

        Type = 'ROD'
        xxb = [0.]
        dims = [[1.]]
        model.add_pbeaml(2000, mid, Type, xxb, dims)
        eids = [1, 2, 3, 4]
        split_line_elements(model, eids, neids=10,
                            eid_start=101, nid_start=101)
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        #print(bdf_file.getvalue())

    def test_shells_add(self):
        """
        tests differential mass and material coordinate systems
        on CQUAD4/CTRIA3 elements

        """
        pid = 10
        mid1 = 100
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(10, pid, [1, 2, 3, 4])
        model.add_ctria3(11, pid, [1, 2, 3])

        mids = [100, 100, 100]
        thicknesses = [0.1, 0.1, 0.1]
        model.add_pcomp(pid, mids, thicknesses, thetas=[0., 45., 90.], souts=None,
                        nsm=0., sb=0., ft=None,
                        tref=0., ge=0., lam=None,
                        z0=None, comment='')

        pid = 11
        model.add_ctria3(12, pid, [1, 2, 3], theta_mcid=45., zoffset=0.,
                         tflag=0, T1=0.1, T2=0.1, T3=0.1,  # absolute - mass=0.1*0.5=0.05
                         comment='')
        model.add_ctria3(13, pid, [1, 2, 3], theta_mcid=1, zoffset=0.,
                         tflag=0, T1=0.1, T2=0.1, T3=0.1,  # absolute
                         comment='')

        model.add_cquad4(14, pid, [1, 2, 3, 4], theta_mcid=45., zoffset=0.,
                         tflag=0, T1=0.1, T2=0.1, T3=0.1, T4=0.1,  # absolute
                         comment='')
        model.add_cquad4(15, pid, [1, 2, 3, 4], theta_mcid=1, zoffset=0.,
                         tflag=1, T1=0.1, T2=0.1, T3=0.1, T4=0.1,  # relative
                         comment='')

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        model.add_cord2r(1, origin, zaxis, xzplane, rid=0)
        model.add_pshell(pid, mid1=mid1, t=2.)

        e11 = 1.0
        e22 = 2.0
        nu12 = 0.3
        model.add_mat8(mid1, e11, e22, nu12, rho=1.0)
        model.validate()

        model.cross_reference()
        model.pop_xref_errors()

        mass = mass_properties(model, element_ids=13)[0]
        bdf_file = StringIO()
        model.write_bdf(bdf_file)
        model.uncross_reference()
        model.cross_reference()
        model.pop_xref_errors()

        assert np.allclose(mass, 0.05), mass # t=0.1; A=0.5; nsm=0.; mass=0.05

        mass = mass_properties(model, element_ids=14)[0]
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)
        assert np.allclose(mass, 0.1), mass # t=0.1; A=1.0; nsm=0.; mass=0.1

        csv_filename = DIRNAME / 'mcids.csv'
        export_mcids(model, csv_filename=csv_filename, eids=[12, 13],
                     export_xaxis=True, export_yaxis=True,
                     iply=0)
        #with open(csv_filename, 'r') as csv_file:
            #lines = csv_file.readlines()
            #assert len(lines) > 0, 'lines=%s' % lines
            #for line in lines:
                #print(line.rstrip())
        #print('-------------')
        export_mcids(model, csv_filename=csv_filename, eids=[14, 15],
                     export_xaxis=True, export_yaxis=True,
                     iply=0)
        model.uncross_reference()
        model.safe_cross_reference()
        model.uncross_reference()
        os.remove(csv_filename)
        #bdf_file = model.write_bdf(bdf_file)

        model2 = BDF(debug=False)
        model2.read_bdf(bdf_file, punch=True)
        #with open(csv_filename, 'r') as csv_file:
            #lines = csv_file.readlines()
            #assert len(lines) > 0, 'lines=%s' % lines
            #for line in lines:
                #print(line.rstrip())

    def test_mirror(self):
        """tests bdf mirroring"""
        log = SimpleLogger(level='error')
        pid_pshell = 10
        pid_psolid = 11
        mid1 = 100
        model = BDF(log=log) # (log=log)
        model.add_grid(1, [10., 10., 10.])
        model.add_grid(2, [11., 10., 10.])
        model.add_grid(3, [11., 11., 10.])
        model.add_grid(4, [10., 11., 10.])

        model.add_grid(5, [10., 10., 11.])
        model.add_grid(6, [11., 10., 11.])
        model.add_grid(7, [11., 11., 11.])
        model.add_grid(8, [10., 11., 11.])

        nodes = [1, 4]
        components = ['123', '123']
        coefficients = [1.0, 1.0]
        mpc = model.add_mpc(42, nodes, components, coefficients, comment='')
        mpc.validate()

        model.add_cquad4(1, pid_pshell, [1, 2, 3, 4]) # mass=1
        model.add_ctria3(2, pid_pshell, [1, 2, 3]) # mass=0.5
        model.add_conrod(3, mid1, [1, 3], A=1.0, j=0.0, c=0.0, nsm=0.0)

        #model.add_ctetra(4, pid_psolid, [1, 2, 3, 5])
        # penta
        # pyram
        #model.add_chexa(7, pid_psolid, [1, 2, 3, 4, 5, 6, 7, 8])

        model.add_pshell(pid_pshell, mid1=mid1, t=1.)
        model.add_psolid(pid_psolid, mid1)
        E = 1.0
        G = None
        nu = 0.3
        model.add_mat1(mid1, E, G, nu, rho=1.0)
        model.validate()
        model.cross_reference()
        mass1, unused_cg1, unused_inertia1 = mass_properties(model)

        # mirror_model=None -> new model
        #
        # just a cord2r
        #   y+ right
        plane = np.array([
            [0., 0., 0.],
            [0., 0., 1.],
            [1., 0., 0.],
        ])
        assert len(model.mpcs) == 1, model.mpcs
        assert len(model.mpcs[42]) == 1, model.mpcs[42]
        model, unused_mirror_model, unused_nid_offset, unused_eid_offset = bdf_mirror_plane(
            model, plane, mirror_model=None, log=None, debug=True,
            use_nid_offset=False)
        #for nid, node in sorted(mirror_model.nodes.items()):
            #print(nid, node.xyz)

        assert len(model.mpcs) == 1, model.mpcs
        assert len(model.mpcs[42]) == 2, model.mpcs[42]

        out_filename = DIRNAME / 'sym.bdf'
        write_bdf_symmetric(model, out_filename=out_filename, encoding=None, size=8,
                            is_double=False,
                            enddata=None,
                            close=True, plane='xz') # +y/-y
        # ----------------------------------------------
        # validate
        model2 = read_bdf(out_filename, log=log)
        assert len(model2.nodes) == 16, model2.nodes
        mass2, cg2, unused_inertia2 = mass_properties(model2)
        #print('cg1=%s cg2=%s' % (cg1, cg2))
        assert np.allclose(mass1*2, mass2), 'mass1=%s mass2=%s' % (mass1, mass2)
        assert np.allclose(cg2[1], 0.), 'cg2=%s stats=%s' % (cg2, model2.get_bdf_stats())

        assert len(model2.mpcs) == 1, model2.mpcs
        assert len(model2.mpcs[42]) == 4, model2.mpcs[42]
        os.remove(out_filename)

    def test_mirror_tetra(self):
        """tests mirroring a chexa"""
        model = BDF(debug=False, log=None, mode='msc')
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [0., 1, 0.])
        model.add_grid(14, [0., 0., 6.])

        pid = 20
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3

        eid_tetra = 10
        nids = [11, 12, 13, 14]
        model.add_ctetra(eid_tetra, pid, nids, comment='')
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        model.validate()
        model.cross_reference()

        bdf_mirror(model, plane='xz', log=None, debug=True)
        model.cross_reference()
        #for eid, elem in model.elements.items():
            #print(eid, elem.Volume())
        x = 1

    def test_mirror_hexa(self):
        """tests mirroring a chexa"""
        model = BDF(debug=False, log=None, mode='msc')

        xmax = 5.
        ymax = 10.
        zmax = 20.
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [xmax, 0., 0.])
        model.add_grid(3, [xmax, ymax, 0.])
        model.add_grid(4, [0., ymax, 0.])
        model.add_grid(5, [0., 0., zmax])
        model.add_grid(6, [xmax, 0., zmax])
        model.add_grid(7, [xmax, ymax, zmax])
        model.add_grid(8, [0., ymax, zmax])

        pid = 20
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3

        eid_hexa = 11
        nids = [1, 2, 3, 4, 5, 6, 7, 8]
        model.add_chexa(eid_hexa, pid, nids, comment='')
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        model.validate()
        model.cross_reference()

        bdf_mirror(model, plane='xz', log=None, debug=True)
        model.cross_reference()
        #for eid, elem in model.elements.items():
            #print(eid, elem.Volume())
        x = 1

    def test_mirror_bwb(self):
        """mirrors the BDF (we care about the aero cards)"""
        log = SimpleLogger(level='warning')
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        model = bdf_mirror(bdf_filename, plane='xz', log=log)[0]
        model.uncross_reference()
        model.cross_reference()
        make_half_model(model, plane='xz', zero_tol=1e-12)
        #model.validate()

    def test_pierce_model(self):
        """tests pierce_shell_model"""
        log = SimpleLogger(level='error')
        pid = 10
        mid1 = 100
        model = BDF(log=log)

        # intersects (min)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, pid, [1, 2, 3, 4])

        # intersects (max)
        model.add_grid(5, [0., 0., 1.])
        model.add_grid(6, [1., 0., 1.])
        model.add_grid(7, [1., 1., 1.])
        model.add_grid(8, [0., 1., 1.])
        model.add_cquad4(2, pid, [5, 6, 7, 8])

        # intersects (mid)
        model.add_grid(9, [0., 0., 0.5])
        model.add_grid(10, [1., 0., 0.5])
        model.add_grid(11, [1., 1., 0.5])
        model.add_grid(12, [0., 1., 0.5])
        model.add_cquad4(3, pid, [9, 10, 11, 12])

        # doesn't intersect
        model.add_grid(13, [10., 0., 0.])
        model.add_grid(14, [11., 0., 0.])
        model.add_grid(15, [11., 1., 0.])
        model.add_grid(16, [10., 1., 0.])
        model.add_cquad4(4, pid, [13, 14, 15, 16])

        model.add_pshell(pid, mid1=mid1, t=2.)

        E = 1.0
        G = None
        nu = 0.3
        model.add_mat1(mid1, E, G, nu, rho=1.0)
        model.validate()

        model.cross_reference()

        xyz_points = [
            [0.4, 0.6, 0.],
            [-1., -1, 0.],
        ]
        pierce_shell_model(model, xyz_points)

    #def test_intersect(self):
        #p0 = np.array([0,0,0], 'd')
        #p1 = np.array([1,0,0], 'd')
        #p2 = np.array([0,1,0], 'd')
        #p3 = np.array([1,1,0], 'd')

        #v = np.array([0,0,-1], 'd')
        #for i in range(10):
            #for j in range(10):
                #p = np.array([i*.2-.5, j*.2-.5, 1.], 'd')
                #print(i, j, p,
                      #triangle_intersection(p, v, p0, p1, p2),
                      #quad_intersection(p, v, p0, p1, p3, p2))

    def test_get_oml_eids(self):
        bdf_filename = BWB_PATH / 'bwb_saero.bdf'
        eid_start = 10144
        model, eids_oml = get_oml_eids(
            bdf_filename, eid_start, theta_tol=30.,
            is_symmetric=True, consider_flippped_normals=True)
        assert len(eids_oml) == 169, len(eids_oml)
        eids_expected = {
            10250, 10251, 10252, 10253, 10254, 10255, 10256, 10257, 10258, 10259, 10260,
            10261, 10262, 10263, 10264, 10265, 10266, 10267, 10268, 10269, 10270, 10271,
            10272, 10273, 10274, 10275, 10276, 10277, 10278, 10279, 10280, 10281, 10282,
            10283, 10284, 10285, 10286, 10287, 10288, 10289, 10290, 10291, 10292, 10293,
            10294, 10295, 10296, 10297, 10298, 10299, 10300, 10301, 10302, 10303, 10304,
            10305, 10306, 10307, 10308, 10309, 10310, 10311, 10312, 10313, 10314, 10315,
            10316, 10317, 10318, 10319, 10320, 10321, 10322, 10323, 10324, 10325, 10326,
            10327, 10328, 10329, 10330, 10331, 10332, 10333, 10334, 10335, 10336, 10337,
            10338, 10339, 10340, 10341, 10342, 10343, 10344, 10345, 10346, 10347, 10348,
            10349, 10350, 10351, 10352, 10353, 10354, 10355, 10356, 10357, 10358, 10359,
            10360, 10361, 10362, 10363, 10364, 10365, 10144, 10145, 10146, 10147, 10148,
            10149, 10150, 10151, 10152, 10153, 10154, 10155, 10156, 10157, 10158, 10159,
            10160, 10161, 10162, 10163, 10164, 10165, 10166, 10167, 10168, 10169, 10196,
            10197, 10198, 10199, 10200, 10201, 10202, 10203, 10204, 10205, 10206, 10207,
            10208, 10209, 10210, 10211, 10212, 10213, 10214, 10215, 10216, 10217, 10218,
            10219, 10220, 10221, 10222}
        assert eids_oml == eids_expected


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
