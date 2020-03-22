"""various mesh_utils tests"""
import os
import unittest
from io import StringIO

from docopt import DocoptExit
import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.split_elements import split_line_elements
from pyNastran.bdf.mesh_utils.pierce_shells import (
    pierce_shell_model) #, quad_intersection, triangle_intersection)
from pyNastran.bdf.mesh_utils.mirror_mesh import (
    write_bdf_symmetric, bdf_mirror, bdf_mirror_plane)
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  #mass_properties_breakdown
from pyNastran.bdf.mesh_utils.make_half_model import make_symmetric_model
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.bdf.mesh_utils.utils import cmd_line
from pyNastran.bdf.mesh_utils.find_closest_nodes import find_closest_nodes
from pyNastran.bdf.mesh_utils.find_coplanar_elements import find_coplanar_triangles
from pyNastran.bdf.mesh_utils.force_to_pressure import force_to_pressure
from pyNastran.bdf.mesh_utils.free_edges import free_edges, non_paired_edges
from pyNastran.bdf.mesh_utils.get_oml import get_oml_eids

from pyNastran.bdf.mesh_utils.mesh import create_structured_cquad4s, create_structured_chexas

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)
DIRNAME = os.path.dirname(__file__)

class TestMeshUtils(unittest.TestCase):
    """various mesh_utils tests"""

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

    def test_free_edges(self):
        """Finds the free_edges

        4-----3---5
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
        model.add_grid(5, [1., 1., 1.])  # 1,3,5
        model.add_ctria3(1, 1, [1, 2, 3])
        edges1 = free_edges(model, eids=None, maps=None)
        edges2 = non_paired_edges(model, eids=None, maps=None)
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

    def test_free_faces(self):
        """CTETRA10"""
        #bdf free_faces [-d | -l] [-f] [--encoding ENCODE] BDF_FILENAME SKIN_FILENAME
        #with self.assertRaises(SystemExit):
            #cmd_line(argv=['bdf', 'free_faces'])

        bdf_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        #log = get_logger(log=None, level='info', encoding='utf-8')
        skin_filename = os.path.join(DIRNAME, 'skin.bdf')
        cmd_line(argv=['bdf', 'free_faces', bdf_filename, skin_filename], quiet=True)
        os.remove(skin_filename)

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
        model = BDF()
        pid1 = 1
        pid2 = 2

        xmax = 10.
        ymax = 10.
        zmax = 10.
        p1 = [0., 0., 0.]
        p2 = [xmax, 0., 0.]
        p3 = [xmax, ymax, 0.]
        p4 = [0., ymax, 0.]

        p5 = [0., 0., zmax]
        p6 = [xmax, 0., zmax]
        p7 = [xmax, ymax, zmax]
        p8 = [0., ymax, zmax]

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
        rho = 2_710  #kg/m^3
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

    def test_merge_01(self):
        """merges multiple bdfs into a single deck"""
        log = SimpleLogger(level='error')
        bdf_filename1 = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        bdf_filename2 = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        bdf_filename3 = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.bdf')
        bdf_filename4 = os.path.join(MODEL_PATH, 'iSat', 'ISat_Dploy_Sm.dat')
        bdf_filename_out1 = os.path.join(MODEL_PATH, 'bwb', 'BWBsaero_staticbar_8.out')
        bdf_filename_out2 = os.path.join(MODEL_PATH, 'bwb', 'BWBsaero_static_bar_16.out')
        bdf_filename_out3 = os.path.join(MODEL_PATH, 'bwb', 'BWBsaero_staticbar_isat.out')

        bdf_filenames1 = [bdf_filename1, bdf_filename2]
        bdf_filenames2 = [bdf_filename1, bdf_filename2, bdf_filename3, bdf_filename4]
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

    def test_exit(self):
        """tests totally failing to run"""
        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'export_caero_mesh'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'convert'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'scale'])

        #with self.assertRaises(SystemExit):
            #cmd_line(argv=['bdf', 'bin'])

        #with self.assertRaises(SystemExit):
            #cmd_line(argv=['bdf', 'filter'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'mirror'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'renumber'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'equivalence'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'free_faces'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'merge'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'export_caero_mesh'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'transform'])

        with self.assertRaises(SystemExit):
            cmd_line(argv=['bdf', 'export_mcids'])

    def test_export_caero_mesh(self):
        """tests multiple ``bdf`` tools"""
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_no_sub.bdf']
        with self.assertRaises(SystemExit):
            cmd_line(argv=argv[:1])
        with self.assertRaises(SystemExit):
            cmd_line(argv=argv[:2])

        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_aesurf.bdf', '--subpanels',
                '--pid', 'aesurf']
        cmd_line(argv=argv, quiet=True)
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_caero.bdf', '--subpanels',
                '--pid', 'caero']
        cmd_line(argv=argv, quiet=True)
        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero_paero.bdf', '--subpanels',
                '--pid', 'paero']
        cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_caero_mesh', bdf_filename, '-o', 'caero.bdf', '--subpanels']
        cmd_line(argv=argv, quiet=True)

        #bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]
        argv = ['bdf', 'mirror', 'caero.bdf', '-o', 'caero2.bdf', '--plane', 'xz', '--tol', '1e-5']
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

    def test_export_mcids(self):
        """creates material coordinate systems"""
        log = SimpleLogger(level='error')
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        csv_filename = os.path.join(MODEL_PATH, 'bwb', 'mcids.csv')
        export_mcids(bdf_filename, csv_filename,
                     export_xaxis=True, export_yaxis=True,
                     iply=9, log=log, debug=False)

        model = read_bdf(bdf_filename, xref=False, debug=False)
        model.safe_cross_reference()
        #os.remove('mcids.csv')

        argv = ['bdf', 'export_mcids', bdf_filename, '-o', csv_filename,
                '--iplies', '0,1,2,3,4,5,6,7,8,9,10', '--no_x', '--no_y']
        with self.assertRaises(DocoptExit):
            # can't define both --no_x and --no_y
            cmd_line(argv=argv, quiet=True)

        argv = ['bdf', 'export_mcids', bdf_filename, '-o', csv_filename,
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

        csv_filename = 'mcids.csv'
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
        model, unsed_mirror_model, unused_nid_offset, unused_eid_offset = bdf_mirror_plane(
            model, plane, mirror_model=None, log=None, debug=True,
            use_nid_offset=False)
        #for nid, node in sorted(mirror_model.nodes.items()):
            #print(nid, node.xyz)


        out_filename = os.path.join(DIRNAME, 'sym.bdf')
        write_bdf_symmetric(model, out_filename=out_filename, encoding=None, size=8,
                            is_double=False,
                            enddata=None,
                            close=True, plane='xz') # +y/-y
        model2 = read_bdf(out_filename, log=log)
        assert len(model2.nodes) == 16, model2.nodes
        mass2, cg2, unused_inertia2 = mass_properties(model2)
        #print('cg1=%s cg2=%s' % (cg1, cg2))
        assert np.allclose(mass1*2, mass2), 'mass1=%s mass2=%s' % (mass1, mass2)
        assert np.allclose(cg2[1], 0.), 'cg2=%s stats=%s' % (cg2, model2.get_bdf_stats())
        os.remove(out_filename)

    def test_mirror2(self):
        """mirrors the BDF (we care about the aero cards)"""
        log = SimpleLogger(level='warning')
        bdf_filename = os.path.join(MODEL_PATH, 'bwb', 'bwb_saero.bdf')
        model = bdf_mirror(bdf_filename, plane='xz', log=log)[0]
        model.uncross_reference()
        model.cross_reference()
        make_symmetric_model(model, plane='xz', zero_tol=1e-12)
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


class TestEquiv(unittest.TestCase):

    def test_closest(self):
        """Finds the closest nodes to specified points"""
        log = SimpleLogger(level='error')
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID,2,,0.,0.,0.5\n'
            'GRID,3,,0.,0.,0.51\n'
            'GRID,10,,0.,0.,1.\n'
            'GRID,11,,0.,0.,1.\n'
            'CTRIA3,1,1,1,2,11\n'
            'CTRIA3,3,1,2,3,11\n'
            'CTRIA3,4,1,1,2,10\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = os.path.join(DIRNAME, 'nonunique.bdf')
        #bdf_filename_out = os.path.join(DIRNAME, 'unique.bdf')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename)
        out = model.get_displacement_index_xyz_cp_cd()
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = out

        nids = nid_cp_cd[:, 0]
        xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nids, icp_transform,
                cid=0)
        nodes_xyz = xyz_cid0

        #'GRID,2,,0.,0.,0.5\n'
        #'GRID,3,,0.,0.,0.51\n'
        xyz_compare = np.array([
            [0., 0., 0.52]])
        nids_close = find_closest_nodes(
            nodes_xyz, nids,
            xyz_compare, neq_max=2, tol=None,
            msg='')
        nids_expected = [
            [3, 2],
        ]
        assert np.array_equal(nids_close, nids_expected), f'A: nids_close={nids_close} nids_expected={nids_expected}'
        # --------------------------------------------------------------
        #'GRID,3,,0.,0.,0.51\n'
        nids_close = find_closest_nodes(
            nodes_xyz, nids,
            xyz_compare, neq_max=1, tol=None,
            msg='')
        nids_expected = [3]
        assert np.array_equal(nids_close, nids_expected), f'B: nids_close={nids_close} nids_expected={nids_expected}'

    def test_eq1(self):
        """Collapse nodes 2 and 3; consider 1-3"""
        log = SimpleLogger(level='error')
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID,2,,0.,0.,0.5\n'
            'GRID,3,,0.,0.,0.51\n'
            'GRID,10,,0.,0.,1.\n'
            'GRID,11,,0.,0.,1.\n'
            'CTRIA3,1,1,1,2,11\n'
            'CTRIA3,3,1,2,3,11\n'
            'CTRIA3,4,1,1,2,10\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = os.path.join(DIRNAME, 'nonunique.bdf')
        bdf_filename_out = os.path.join(DIRNAME, 'unique.bdf')

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=3, skip_cards=['CTRIA3'])
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 2, 10], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=3, skip_cards=['CTRIA3'])
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 2, 10], node_ids

        os.remove(bdf_filename)

    def test_eq2(self):
        r"""
          5
        6 *-------* 40
          | \     |
          |   \   |
          |     \ |
          *-------* 3
          1       20

        """
        log = SimpleLogger(level='error')
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1, , 0.,   0.,   0.\n'
            'GRID,20,, 1.,   0.,   0.\n'
            'GRID,3, , 1.01, 0.,   0.\n'
            'GRID,40,, 1.,   1.,   0.\n'
            'GRID,5, , 0.,   1.,   0.\n'
            'GRID,6, , 0.,   1.01, 0.\n'
            'CTRIA3,1, 100,1,20,6\n'
            'CTRIA3,10,100,3,40,5\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = os.path.join(DIRNAME, 'nonunique.bdf')
        bdf_filename_out = os.path.join(DIRNAME, 'unique.bdf')

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=4)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 40], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=4)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 40], node_ids

        tol = 0.009
        # Don't collapse anything because the tolerance is too small
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=6)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 20, 40], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=6)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 20, 40], node_ids

        tol = 0.2
        node_set = [2, 3]
        # Node 2 is not defined, so crash
        with self.assertRaises(RuntimeError):
            # node 2 is not defined because it should be node 20
            bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                                  renumber_nodes=False, neq_max=4, xref=True,
                                  node_set=node_set, crash_on_collapse=False,
                                  log=log, debug=False)

        tol = 0.2
        node_list = [20, 3]
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_list, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids

        tol = 0.2
        node_list = [20, 3]
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_list, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids

        tol = 0.2
        node_set = {20, 3}
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids

        tol = 0.2
        aset = np.array([20, 3, 4], dtype='int32')
        bset = np.array([20, 3], dtype='int32')

        node_set = np.intersect1d(aset, bset)
        assert len(node_set) > 0, node_set
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False, debug=False,
                              method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False, debug=False,
                              method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 5, 6, 40], node_ids


    def test_eq3(self):
        """node_set=None"""
        log = SimpleLogger(level='error')
        lines = [
            '$pyNastran: version=msc',
            '$pyNastran: punch=True',
            '$pyNastran: encoding=ascii',
            '$NODES',
            '$ Nodes to merge:',
            '$ 5987 10478',
            '$   GRID        5987           35.46     -6.      0.',
            '$   GRID       10478           35.46     -6.      0.',
            '$ 5971 10479',
            '$   GRID        5971           34.92     -6.      0.',
            '$   GRID       10479           34.92     -6.      0.',
            '$ 6003 10477',
            '$   GRID        6003             36.     -6.      0.',
            '$   GRID       10477             36.     -6.      0.',
            'GRID        5971           34.92     -6.      0.',
            'GRID        5972           34.92-5.73333      0.',
            'GRID        5973           34.92-5.46667      0.',
            'GRID        5987           35.46     -6.      0.',
            'GRID        5988           35.46-5.73333      0.',
            'GRID        5989           35.46-5.46667      0.',
            'GRID        6003             36.     -6.      0.',
            'GRID        6004             36.-5.73333      0.',
            'GRID        6005             36.-5.46667      0.',
            'GRID       10476             36.     -6.    -1.5',
            'GRID       10477             36.     -6.      0.',
            'GRID       10478           35.46     -6.      0.',
            'GRID       10479           34.92     -6.      0.',
            'GRID       10561           34.92     -6.    -.54',
            '$ELEMENTS_WITH_PROPERTIES',
            'PSHELL         1       1      .1',
            'CQUAD4      5471       1    5971    5987    5988    5972',
            'CQUAD4      5472       1    5972    5988    5989    5973',
            'CQUAD4      5486       1    5987    6003    6004    5988',
            'CQUAD4      5487       1    5988    6004    6005    5989',
            'PSHELL        11       1      .1',
            'CTRIA3      9429      11   10561   10476   10478',
            'CTRIA3      9439      11   10478   10479   10561',
            'CTRIA3      9466      11   10476   10477   10478',
            '$MATERIALS',
            'MAT1           1      3.              .3',
        ]
        bdf_filename = os.path.join(DIRNAME, 'nonunique2.bdf')
        bdf_filename_out = os.path.join(DIRNAME, 'unique2.bdf')

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(lines))

        tol = 0.01
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=11)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [5971, 5972, 5973, 5987, 5988, 5989, 6003, 6004, 6005, 10476, 10561], node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=11)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [5971, 5972, 5973, 5987, 5988, 5989, 6003, 6004, 6005, 10476, 10561], node_ids

        os.remove(bdf_filename)

    def test_eq4(self):
        r"""
          5
        6 *-------* 40
          | \     |
          |   \   |
          |     \ |
          *-------* 3
          1       20

        """
        log = SimpleLogger(level='error')
        msg = 'CEND\n'
        msg += 'BEGIN BULK\n'
        msg += 'GRID,1, , 0.,   0.,   0.\n'
        msg += 'GRID,20,, 1.,   0.,   0.\n'
        msg += 'GRID,3, , 1.01, 0.,   0.\n'

        msg += 'GRID,41,, 1.,   1.,   0.\n'  # eq
        msg += 'GRID,4,, 1.,   1.,   0.\n'  # eq
        msg += 'GRID,40,, 1.,   1.,   0.\n'  # eq
        msg += 'GRID,4,, 1.,   1.,   0.\n'  # eq

        msg += 'GRID,5, , 0.,   1.,   0.\n'
        msg += 'GRID,6, , 0.,   1.01, 0.\n'
        msg += 'CTRIA3,1, 100,1,20,6\n'
        msg += 'CTRIA3,10,100,3,40,5\n'
        msg += 'PSHELL,100,1000,0.1\n'
        msg += 'MAT1,1000,3.0,, 0.3\n'
        msg += 'ENDDATA'
        bdf_filename = os.path.join(DIRNAME, 'nonunique.bdf')
        bdf_filename_out = os.path.join(DIRNAME, 'unique.bdf')

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        node_set = [4, 40, 41]
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='old')

        model = save_check_nodes(bdf_filename_out, log, nnodes=6)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1, 3, 4, 5, 6, 20], node_ids
        os.remove(bdf_filename)

    def test_eq5(self):
        log = SimpleLogger(level='info')
        bdf_filename_out = 'eq5.bdf'

        model = BDF(debug=True, log=log, mode='msc')
        for nid in range(1, 11):
            model.add_grid(nid, [0., 0., 0.])
        eid = 1
        k = 2.0
        nids = [10, None]
        model.add_celas2(eid, k, nids, c1=2, c2=0, ge=0., s=0., comment='')

        tol = 1.0
        node_set = None
        eid = 1
        k = 2.0
        nids = [10, None]
        bdf_equivalence_nodes(model, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=6)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [5, 6, 7, 8, 9, 10], node_ids
        del model

        model2 = BDF(debug=True, log=log, mode='msc')
        for nid in range(1, 11):
            model2.add_grid(nid, [0., 0., 0.])
        model2.add_celas2(eid, k, nids, c1=2, c2=0, ge=0., s=0., comment='')
        bdf_equivalence_nodes(model2, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=True, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=1)
        node_ids = list(sorted(model.nodes))
        assert node_ids == [1], node_ids


def save_check_nodes(bdf_filename, log, nnodes, skip_cards=None):
    model = BDF(log=log, debug=False)
    model.disable_cards(skip_cards)
    model.read_bdf(bdf_filename)
    msg = 'nnodes=%s\n' % len(model.nodes)
    for nid, node in sorted(model.nodes.items()):
        msg += 'nid=%s xyz=%s\n' % (nid, node.xyz)
    assert len(model.nodes) == nnodes, msg
    os.remove(bdf_filename)
    return model


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
