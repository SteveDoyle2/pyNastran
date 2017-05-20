from __future__ import print_function
import os
import unittest
from codecs import open as codec_open

from six import StringIO, iteritems
import numpy as np
#import pyNastran
#from pyNastran.bdf.bdf import BDF

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import get_bad_shells
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.split_elements import split_line_elements
from pyNastran.utils.log import SimpleLogger

# testing these imports are up to date
from pyNastran.bdf.mesh_utils.utils import *

pkg_path = pyNastran.__path__[0]

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)

log = SimpleLogger(level='error')
class TestMeshUtils(unittest.TestCase):

    def test_quad_180_01(self):
        r"""
        Identify a 180+ degree quad

        y
        ^         4
        |       / |
        |     /   |
        |   /     |
        | /       |
        /         |
        1------2  |----> x
                \ |
                 \|
                  3
        """
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID,2,,1.,0.,0.\n'
            'GRID,3,,2.,-1.,0.\n'
            'GRID,4,,2., 1.,0.\n'

            'CQUAD4,100,1, 1,2,3,4\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = 'cquad4.bdf'
        with codec_open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, log=log, xref=True)
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
        nid_map = {}
        for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
            #xyz = node.get_position()
            #xyz_cid0[i, :] = xyz
            nid_map[nid] = i
        eids_to_delete = get_bad_shells(model, xyz_cid0, nid_map, max_theta=180.,
                                        max_skew=1000., max_aspect_ratio=1000.)
        assert eids_to_delete == [100], eids_to_delete
        os.remove(bdf_filename)

    def test_eq1(self):
        """
        Collapse nodes 2 and 3; consider 1-3
        """
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

        with codec_open(bdf_filename, 'w') as bdf_file:
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

        with codec_open(bdf_filename, 'w') as bdf_file:
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
        for nid, node in sorted(iteritems(model.nodes)):
            msg += 'nid=%s xyz=%s\n' % (nid, node.xyz)

        assert len(model.nodes) == 4, msg
        # os.remove(bdf_filename)
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

        with codec_open(bdf_filename, 'w') as bdf_file:
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

        with codec_open(bdf_filename, 'w') as bdf_file:
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

    def test_fix_bad_quads(self):
        """split high interior angle quads"""
        msg = [
            'SOL 101',
            'CEND',
            'BEGIN BULK',
            'GRID,1,,0.,0.,0.',
            'GRID,2,,1.,0.,0.',
            'GRID,3,,1.,1.,0.',
            'GRID,4,,0.,1.,0.',
            'GRID,5,,1.,1.,0.00001',
            'GRID,6,,0.,0.,0.00001',
            'CQUAD4,1, 2, 1,2,3,4',
            'CQUAD4,2, 2, 1,2,3,5',
            'CQUAD4,3, 2, 1,6,3,5',
        ]
        bdf_filename = 'fix_bad_quads.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(msg))
        model = read_bdf(bdf_filename, log=log, xref=False, debug=False)
        model.cross_reference(xref=True, xref_elements=False,
                              xref_nodes_with_elements=False,
                              xref_properties=False,
                              xref_masses=False,
                              xref_materials=False,
                              xref_loads=False,
                              xref_constraints=False,
                              xref_aero=False,
                              xref_sets=False,
                              xref_optimization=False)
        convert_bad_quads_to_tris(model, min_edge_length=0.01)
        #for eid, elem in sorted(iteritems(model.elements)):
            #print(elem)
        assert model.card_count['CQUAD4'] == 2, model.card_count
        assert model.card_count['CTRIA3'] == 1, model.card_count
        os.remove(bdf_filename)

    def test_renumber_01(self):
        """renumbers a deck in a couple ways"""
        bdf_filename = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWB_saero.bdf'))
        bdf_filename_out1 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWB_saero1.out'))
        bdf_filename_out2 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWB_saero2.out'))
        bdf_filename_out3 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWB_saero3.out'))
        model = bdf_renumber(bdf_filename, bdf_filename_out1, size=8,
                             is_double=False, starting_id_dict=None,
                             round_ids=False, cards_to_skip=None)

        model = read_bdf(bdf_filename, log=log)
        bdf_renumber(model, bdf_filename_out2, size=16, is_double=False,
                     starting_id_dict={
                         'eid' : 1000, 'pid':2000, 'mid':3000,
                         'spc_id' : 4000,},
                     round_ids=False, cards_to_skip=None)
        bdf_renumber(bdf_filename, bdf_filename_out3, size=8,
                     is_double=False, starting_id_dict=None,
                     round_ids=True, cards_to_skip=None)
        read_bdf(bdf_filename_out1, log=log)
        read_bdf(bdf_filename_out2, log=log)
        read_bdf(bdf_filename_out3, log=log)

    def test_merge_01(self):
        """merges multiple bdfs into a single deck"""
        #log = SimpleLogger(level='info')
        bdf_filename1 = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'bwb', 'BWB_saero.bdf'))
        bdf_filename2 = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'sol_101_elements', 'static_solid_shell_bar.bdf'))
        bdf_filename3 = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
        bdf_filename4 = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'iSat', 'ISat_Dploy_Sm.dat'))

        bdf_filename_out1 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWBsaero_staticbar_8.out'))
        bdf_filename_out2 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWBsaero_static_bar_16.out'))
        bdf_filename_out3 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'BWBsaero_staticbar_isat.out'))

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
        bdf_filename = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'bwb', 'BWB_saero.bdf'))
        csv_filename = os.path.abspath(os.path.join(
            pkg_path, '..', 'models', 'bwb', 'mcids.csv'))
        export_mcids(bdf_filename, csv_filename,
                     export_xaxis=True, export_yaxis=True,
                     iply=9)

        model = read_bdf(bdf_filename, xref=False)
        model.safe_cross_reference()

        eids = [1204, 1211]
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=True, export_yaxis=True,
                     iply=9)
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=True, export_yaxis=False,
                     iply=9)
        export_mcids(model, csv_filename=None, eids=eids,
                     export_xaxis=False, export_yaxis=True,
                     iply=9)
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
        bdf_filename = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'sol_101_elements', 'static_solid_shell_bar.bdf'))
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

        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])
        model.add_grid(3, xyz=[2., 0., 0.])
        model.add_grid(4, xyz=[3., 0., 0.])

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
        model = BDF()
        model.add_grid(1, xyz=[0., 0., 0.])
        model.add_grid(2, xyz=[1., 0., 0.])

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
        f = StringIO()
        model.write_bdf(f, close=False)
        #print(f.getvalue())

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
