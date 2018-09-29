"""various mesh_utils tests"""
from __future__ import print_function
import os
import unittest
from codecs import open

from six import StringIO
import numpy as np
#import pyNastran
#from pyNastran.bdf.bdf import BDF

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, CORD2R
from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.split_elements import split_line_elements
from pyNastran.bdf.mesh_utils.pierce_shells import (
    pierce_shell_model) #, quad_intersection, triangle_intersection)
from pyNastran.bdf.mesh_utils.mirror_mesh import (
    write_bdf_symmetric, bdf_mirror, make_symmetric_model)
from pyNastran.bdf.mesh_utils.cut_model_by_plane import (
    cut_edge_model_by_coord, cut_face_model_by_coord, connect_face_rows)
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.utils.log import SimpleLogger

# not tested
from pyNastran.bdf.mesh_utils.mesh import create_structured_cquad4s

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)


class TestMeshUtils(unittest.TestCase):

    def test_eq1(self):
        """
        Collapse nodes 2 and 3; consider 1-3
        """
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
        bdf_filename = 'nonunique.bdf'
        bdf_filename_out = 'unique.bdf'

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False)

        # model = BDF(debug=False)
        # model.read_bdf(bdf_filename_out)
        # assert len(model.nodes) == 3, len(model.nodes)

        os.remove(bdf_filename)
        os.remove(bdf_filename_out)

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
        bdf_filename = 'nonunique.bdf'
        bdf_filename_out = 'unique.bdf'

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False)

        model = BDF(log=log, debug=False)
        model.read_bdf(bdf_filename_out)

        msg = 'nnodes=%s\n' % len(model.nodes)
        for nid, node in sorted(model.nodes.items()):
            msg += 'nid=%s xyz=%s\n' % (nid, node.xyz)

        assert len(model.nodes) == 4, msg
        #os.remove(bdf_filename)
        os.remove(bdf_filename_out)

        tol = 0.009
        # Don't collapse anything because the tolerance is too small
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False)
        model = BDF(log=log, debug=False)
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 6, len(model.nodes)
        os.remove(bdf_filename_out)

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
                              log=log, debug=False)
        model = BDF(log=log, debug=False)
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 5, len(model.nodes)
        os.remove(bdf_filename_out)

        tol = 0.2
        node_set = set([20, 3])
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False)
        model = BDF(log=log, debug=False)
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 5, len(model.nodes)
        os.remove(bdf_filename_out)

        tol = 0.2
        aset = np.array([20, 3, 4], dtype='int32')
        bset = np.array([20, 3], dtype='int32')

        node_set = np.intersect1d(aset, bset)
        assert len(node_set) > 0, node_set
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False, debug=False)
        model = BDF(debug=False)
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 5, len(model.nodes)
        os.remove(bdf_filename_out)


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
        bdf_filename = 'nonunique2.bdf'
        bdf_filename_out = 'unique2.bdf'

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(lines))

        tol = 0.01
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False)

        model = BDF(debug=False)
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 11, len(model.nodes)

        os.remove(bdf_filename)
        os.remove(bdf_filename_out)

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
        bdf_filename = 'nonunique.bdf'
        bdf_filename_out = 'unique.bdf'

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        node_set = [4, 40, 41]
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False)

        model = BDF(log=log, debug=False)
        model.read_bdf(bdf_filename_out)
        nids = model.nodes.keys()
        assert len(model.nodes) == 6, 'nnodes=%s nodes=%s' % (len(model.nodes), nids)
        assert 1 in nids, nids
        assert 20 in nids, nids
        assert 3 in nids, nids
        assert 4 in nids, nids
        assert 5 in nids, nids
        assert 6 in nids, nids
        assert 40 not in nids, nids
        assert 41 not in nids, nids
        #print(nids)
        os.remove(bdf_filename)
        os.remove(bdf_filename_out)

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
            export_mcids(model, csv_filename=None, eids=eids,
                         export_xaxis=False, export_yaxis=False,
                         iply=9)

        with self.assertRaises(RuntimeError):
            export_mcids(model, csv_filename, eids=eids,
                         export_xaxis=True, export_yaxis=True,
                         iply=10)

    def test_split_cbars_by_pin_flag_1(self):
        """null pin flag test"""
        bdf_filename = os.path.join(MODEL_PATH, 'sol_101_elements', 'static_solid_shell_bar.bdf')
        split_cbars_by_pin_flag(bdf_filename, pin_flags_filename='pin_flags.csv',
                                bdf_filename_out='pin_flags.bdf', debug=False)
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

        mass = model.mass_properties(element_ids=13)[0]
        bdf_file = StringIO()
        model.write_bdf(bdf_file)
        model.uncross_reference()
        model.cross_reference()
        model.pop_xref_errors()

        assert np.allclose(mass, 0.05), mass # t=0.1; A=0.5; nsm=0.; mass=0.05

        mass = model.mass_properties(element_ids=14)[0]
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
        mass1, unused_cg1, unused_inertia1 = model.mass_properties()

        out_filename = 'sym.bdf'
        write_bdf_symmetric(model, out_filename=out_filename, encoding=None, size=8,
                            is_double=False,
                            enddata=None,
                            close=True, plane='xz') # +y/-y
        model2 = read_bdf(out_filename, log=log)
        assert len(model2.nodes) == 16, model2.nodes
        mass2, cg2, unused_inertia2 = model2.mass_properties()
        #print('cg1=%s cg2=%s' % (cg1, cg2))
        assert np.allclose(mass1*2, mass2), 'mass1=%s mass2=%s' % (mass1, mass2)
        assert np.allclose(cg2[1], 0.), 'cg2=%s stats=%s' % (cg2, model2.get_bdf_stats())
        os.remove('sym.bdf')

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

    def test_cut_shell_model_1(self):
        """tests pierce_shell_model"""
        model, nodal_result = _cut_shell_model_quads()
        coord = CORD2R(1, rid=0, origin=[0.5, 0., 0.], zaxis=[0.5, 0., 1], xzplane=[1.5, 0., 0.],
                       comment='')
        model.coords[1] = coord
        tol = 2.
        #-------------------------------------------------------------------------
        local_points_array, global_points_array, result_array = cut_edge_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5)
        assert len(result_array) == 16, len(result_array)

        geometry_array, result_array = cut_face_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5)
        assert result_array is None, len(result_array) # no quad support

    def test_cut_shell_model_2(self):
        """tests pierce_shell_model"""
        tol = 2.
        coord = CORD2R(1, rid=0, origin=[0.5, 0., 0.], zaxis=[0.5, 0., 1], xzplane=[1.5, 0., 0.],
                       comment='')
        model, nodal_result = _cut_shell_model_quads()
        #-------------------------------------------------------------------------
        # triangles
        elements2 = {}
        neids = len(model.elements)
        for eid, elem in model.elements.items():
            elem_a, elem_b = elem.split_to_ctria3(eid, eid + neids)
            elements2[elem_a.eid] = elem_a
            elements2[elem_b.eid] = elem_b
        model.elements = elements2
        model.coords[1] = coord
        model.write_bdf('tris.bdf')

        #print('----------------------------')
        local_points_array, global_points_array, result_array = cut_edge_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5, csv_filename='cut_edge_2.csv')
        assert len(result_array) == 20, len(result_array)

        geometry_arrays, result_arrays = cut_face_model_by_coord(
            model, coord, tol, nodal_result,
            plane_atol=1e-5, csv_filename='cut_face_2.csv')
        assert len(result_arrays[0]) == 8, len(result_arrays)
        os.remove('tris.bdf')


    def test_connect_face_rows(self):
        # in order
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 4, 5],
            [5, 5, 6],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, geometry_arrays2, results_arrays2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 4]]), 'iedges=%s' % iedges
        #-----------------------------------------------------------------------

        # out of order
        geometry_array = np.array([
            [1, 1, 2], # 0
            [2, 4, 5], # 3
            [3, 5, 6], # 4
            [4, 3, 4], # 2
            [5, 2, 3], # 1
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        iedges, geometry_arrays2, results_arrays2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 4, 3, 1, 2]]), 'iedges=%s' % iedges
        #print(geometry_array2)

        #-----------------------------------------------------------------------

        # in order, two blocks
        #print('*****************')
        geometry_array = np.array([
            # block 1
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],

            # block 2
            [10, 10, 20],
            [20, 20, 30],
            [30, 30, 40],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, geometry_array2, results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2], [3, 4, 5]]), 'iedges=%s' % iedges

    def test_connect_face_rows_ring_1(self):
        # in order, one ring
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 1, 4],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, geometry_array2, results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 0]]), 'iedges=%s' % iedges

    def test_connect_face_rows_ring_2(self):
        # in order, two rings
        geometry_array = np.array([
            [1, 1, 2],
            [2, 2, 3],
            [3, 3, 4],
            [4, 1, 4],

            [10, 10, 20],
            [20, 20, 30],
            [30, 30, 40],
            [40, 10, 40],
        ])
        nedges = geometry_array.shape[0]
        results_array = np.arange(0, nedges)
        #print(results_array)
        iedges, geometry_array2, results_array2 = connect_face_rows(
            geometry_array, results_array, skip_cleanup=False)
        assert np.array_equal(iedges, [[0, 1, 2, 3, 0], [4, 5, 6, 7, 4]]), 'iedges=%s' % iedges

def _cut_shell_model_quads():
    """helper method"""
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

    #xyz_points = [
        #[0.4, 0.6, 0.], [-1., -1, 0.],]

    #tol = 2.
    nodal_result = np.linspace(0., 1., num=16)
    return model, nodal_result

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
