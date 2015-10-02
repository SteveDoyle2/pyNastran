from __future__ import unicode_literals, print_function
import unittest
#from six import PY2
#from codecs import open as codec_open

import os
from numpy import array, intersect1d
#import pyNastran
#from pyNastran.bdf.bdf import BDF

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdfInterface.dev_utils import bdf_equivalence_nodes


class DevUtils(unittest.TestCase):
    def test_eq1(self):
        """
        Collapse nodes 2 and 3; consider 1-3
        """
        msg = 'CEND\n'
        msg += 'BEGIN BULK\n'
        msg += 'GRID,1,,0.,0.,0.\n'
        msg += 'GRID,2,,0.,0.,0.5\n'
        msg += 'GRID,3,,0.,0.,0.51\n'
        msg += 'GRID,10,,0.,0.,1.\n'
        msg += 'GRID,11,,0.,0.,1.\n'
        msg += 'CTRIA3,1,1,1,2,11\n'
        msg += 'CTRIA3,3,1,2,3,11\n'
        msg += 'CTRIA3,4,1,1,2,10\n'
        msg += 'PSHELL,1,1,0.1\n'
        msg += 'MAT1,1,3.0,, 0.3\n'
        msg += 'ENDDATA'

        bdf_filename = 'nonunique.bdf'
        bdf_filename_out = 'unique.bdf'

        bdf_file = open(bdf_filename, 'wb')
        bdf_file.write(msg)
        bdf_file.close()

        tol = 0.2
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False)

        # model = BDF()
        # model.read_bdf(bdf_filename_out)
        # assert len(model.nodes) == 3, len(model.nodes)

        os.remove(bdf_filename)
        os.remove(bdf_filename_out)

    def test_eq2(self):
        """
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
        msg += 'GRID,40,, 1.,   1.,   0.\n'
        msg += 'GRID,5, , 0.,   1.,   0.\n'
        msg += 'GRID,6, , 0.,   1.01, 0.\n'
        msg += 'CTRIA3,1, 100,1,20,6\n'
        msg += 'CTRIA3,10,100,3,40,5\n'
        msg += 'PSHELL,100,1000,0.1\n'
        msg += 'MAT1,1000,3.0,, 0.3\n'
        msg += 'ENDDATA'
        bdf_filename = 'nonunique.bdf'
        bdf_filename_out = 'unique.bdf'

        bdf_file = open(bdf_filename, 'wb')
        bdf_file.write(msg)
        bdf_file.close()

        tol = 0.2
        # Collapse 5/6 and 20/3; Put a 40 and 20 to test non-sequential IDs
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False)

        model = BDF()
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 4, len(model.nodes)
        # os.remove(bdf_filename)
        os.remove(bdf_filename_out)

        tol = 0.009
        # Don't collapse anything because the tolerance is too small
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False)
        model = BDF()
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 6, len(model.nodes)
        os.remove(bdf_filename_out)

        tol = 0.2
        node_set = [2, 3]
        # Node 2 is not defined, so crash
        with self.assertRaises(AssertionError):
            # node 2 is not defined because it should be node 20
            bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                                  renumber_nodes=False, neq_max=4, xref=True,
                                  node_set=node_set, crash_on_collapse=False)

        tol = 0.2
        node_set = [20, 3]
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False)
        model = BDF()
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 5, len(model.nodes)
        os.remove(bdf_filename_out)

        tol = 0.2
        node_set = set([20, 3])
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False)
        model = BDF()
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 5, len(model.nodes)
        os.remove(bdf_filename_out)

        tol = 0.2
        aset = array([20, 3, 4], dtype='int32')
        bset = array([20, 3], dtype='int32')

        node_set = intersect1d(aset, bset)
        assert len(node_set) > 0, node_set
        # Only collpase 2 nodes
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=node_set, crash_on_collapse=False)
        model = BDF()
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

        bdf_file = open(bdf_filename, 'wb')
        bdf_file.write('\n'.join(lines))
        bdf_file.close()
        tol = 0.01
        bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                              renumber_nodes=False, neq_max=4, xref=True,
                              node_set=None, crash_on_collapse=False)

        model = BDF()
        model.read_bdf(bdf_filename_out)
        assert len(model.nodes) == 11, len(model.nodes)

        os.remove(bdf_filename)
        os.remove(bdf_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
