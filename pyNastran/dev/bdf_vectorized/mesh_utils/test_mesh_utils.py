"""
tests
 - vectorized renumbering
 - shell mesh quality

"""
import os
import unittest

import numpy as np
from cpylog import SimpleLogger
#import pyNastran

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

import pyNastran
from pyNastran.dev.bdf_vectorized.bdf import read_bdf
from pyNastran.dev.bdf_vectorized.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.dev.bdf_vectorized.mesh_utils.delete_bad_elements import get_bad_shells

pkg_path = pyNastran.__path__[0]

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)


class TestMeshUtilsVectorized(unittest.TestCase):
    """runs the vectorized mesh utils tests"""
    def test_renumber_01(self):
        """renumber a bdf"""
        log = SimpleLogger(level='warning')
        bdf_filename = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero.bdf'))
        bdf_filename_out1 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero1.out'))
        bdf_filename_out2 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero2.out'))
        bdf_filename_out3 = os.path.abspath(
            os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero3.out'))
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
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, xref=True)
        max_theta = 180.
        max_skew = 1000.
        max_aspect_ratio = 1000.
        eids_to_delete = get_bad_shells(
            model,
            max_theta=max_theta,
            max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio,
            max_taper_ratio=4.0)
        assert eids_to_delete == [100], eids_to_delete
        os.remove(bdf_filename)

    def test_tri_180_01(self):
        r"""
        Identify a reasonable tri with super tight tolerances

        y
        ^         4
        |       / /
        |     /   /
        |   /    /
        | /      /
        /       /
        1------2-------> x
        """
        msg = (
            'CEND\n'
            'BEGIN BULK\n'
            'GRID,1,,0.,0.,0.\n'
            'GRID,2,,1.,0.,0.\n'
            'GRID,4,,2., 1.,0.\n'

            'CTRIA3,100,1, 1,2,4\n'
            'PSHELL,1,1,0.1\n'
            'MAT1,1,3.0,, 0.3\n'
            'ENDDATA'
        )
        bdf_filename = 'ctria3.bdf'
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write(msg)

        model = read_bdf(bdf_filename, xref=True)
        eids_to_delete = get_bad_shells(model, max_theta=180.,
                                        max_skew=1000., max_aspect_ratio=1000.)
        assert eids_to_delete == [100], eids_to_delete
        os.remove(bdf_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
