"""various mesh_utils tests"""
import os
import unittest
from pathlib import Path

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.mesh_utils.find_closest_nodes import find_closest_nodes
from pyNastran.dev.bdf_vectorized3.bdf import BDF, read_bdf
from pyNastran.dev.bdf_vectorized3.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = (PKG_PATH / '..' / 'models').resolve()
BWB_PATH = MODEL_PATH / 'bwb'

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)
DIRNAME = PKG_PATH / 'bdf' / 'mesh_utils' / 'test'
assert DIRNAME.exists(), DIRNAME

class TestEquiv(unittest.TestCase):

    def test_closest(self):
        """Finds the closest nodes to specified points"""
        log = SimpleLogger(level='warning')
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
        bdf_filename = DIRNAME / 'nonunique.bdf'
        #bdf_filename_out = os.path.join(DIRNAME, 'unique.bdf')
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, log=log)
        #out = model.get_displacement_index_xyz_cp_cd()
        #icd_transform, icp_transform, xyz_cp, nid_cp_cd = out

        #nids = nid_cp_cd[:, 0]
        nids = model.grid.node_id
        xyz_cid0 = model.grid.xyz_cid0()
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
        bdf_filename = DIRNAME / 'nonunique.bdf'
        bdf_filename_out = DIRNAME / 'unique.bdf'
        if bdf_filename_out.exists():
            os.remove(bdf_filename_out)

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = read_bdf(bdf_filename_out, xref=False)
        model = save_check_nodes(bdf_filename_out, log, nnodes=3, skip_cards=['CTRIA3'])
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 2, 10]), node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=3, skip_cards=['CTRIA3'])
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 2, 10]), node_ids

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
        bdf_filename = DIRNAME / 'nonunique.bdf'
        bdf_filename_out = DIRNAME / 'unique.bdf'
        if bdf_filename_out.exists():
            os.remove(bdf_filename_out)

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        node_list = [20, 3]
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_list, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids

        tol = 0.2
        node_list = [20, 3]
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_list, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids

        tol = 0.2
        node_set = {20, 3}
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids

        tol = 0.2
        aset = np.array([20, 3, 4], dtype='int32')
        bset = np.array([20, 3], dtype='int32')

        node_set = np.intersect1d(aset, bset)
        assert len(node_set) > 0, node_set
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              #log=log,
                              debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False, debug=False,
                              method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=5)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 5, 6, 40]), node_ids


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
        bdf_filename = DIRNAME / 'nonunique2.bdf'
        bdf_filename_out = DIRNAME / 'unique2.bdf'

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('\n'.join(lines))

        tol = 0.01
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = save_check_nodes(bdf_filename_out, log, nnodes=11)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [5971, 5972, 5973, 5987, 5988, 5989, 6003, 6004, 6005, 10476, 10561]), node_ids

        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False,
                              log=log, debug=False, method='new')
        model = save_check_nodes(bdf_filename_out, log, nnodes=11)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [5971, 5972, 5973, 5987, 5988, 5989, 6003, 6004, 6005, 10476, 10561]), node_ids

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
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1, , 0.,   0.,   0.\n'
            'GRID,20,, 1.,   0.,   0.\n'
            'GRID,3, , 1.01, 0.,   0.\n'

            'GRID,41,, 1.,   1.,   0.\n'  # eq
            'GRID,4,, 1.,   1.,   0.\n'  # eq
            'GRID,40,, 1.,   1.,   0.\n'  # eq
            'GRID,4,, 1.,   1.,   0.\n'  # eq

            'GRID,5, , 0.,   1.,   0.\n'
            'GRID,6, , 0.,   1.01, 0.\n'
            'CTRIA3,1, 100,1,20,6\n'
            'CTRIA3,10,100,3,40,5\n'
            'PSHELL,100,1000,0.1\n'
            'MAT1,1000,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = DIRNAME / 'nonunique.bdf'
        bdf_filename_out = DIRNAME / 'unique.bdf'
        if bdf_filename_out.exists():
            os.remove(bdf_filename_out)

        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        tol = 0.2
        node_set = [4, 40, 41]
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=False, method='old')
        model = read_bdf(bdf_filename_out)

        model = save_check_nodes(bdf_filename_out, log, nnodes=6)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1, 3, 4, 5, 6, 20]), node_ids
        os.remove(bdf_filename)

    def test_eq5(self):
        log = SimpleLogger(level='info')
        bdf_filename_out = DIRNAME / 'eq5.bdf'
        if bdf_filename_out.exists():
            os.remove(bdf_filename_out)

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
        model = save_check_nodes(bdf_filename_out, log, nnodes=1)
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1]), node_ids
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
        node_ids = model.grid.node_id
        assert np.array_equal(node_ids, [1]), node_ids

    def test_multi_eq1(self):
        model = BDF(debug=False)
        log = model.log
        # group 1
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])

        # group 2
        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        node_set = [
            [1, 2],
            [3, 4, 5],
        ]
        #nid_pairs = [
            #[1, 2],
            #[2, 3],
            #[1, 3],
            #[4, 5],
        #]
        #[(1, 2), (1, 3), (2, 3), (4, 5)]

        tol = 0.1
        bdf_filename_out = DIRNAME / 'multi_eq1.bdf'
        bdf_equivalence_nodes(model, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False,
                              log=log, debug=True, method='new')
        model2 = read_bdf(bdf_filename_out, debug=None)
        assert len(model2.grid) == 3, model2.grid

def save_check_nodes(bdf_filename, log: SimpleLogger, nnodes: int, skip_cards=None):
    model = BDF(log=log, debug=False)
    model.disable_cards(skip_cards)
    model.read_bdf(bdf_filename)

    grid = model.grid
    msg = 'nnodes=%s\n' % len(grid)
    for nid, xyz in zip(grid.node_id, grid.xyz):
        msg += 'nid=%s xyz=%s\n' % (nid, xyz)
    assert len(grid) == nnodes, msg
    os.remove(bdf_filename)
    return model


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
