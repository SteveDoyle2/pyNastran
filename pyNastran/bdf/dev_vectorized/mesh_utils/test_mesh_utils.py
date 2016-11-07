"""
tests
 - vectorized renumbering
"""
from __future__ import print_function
import os
import unittest

import numpy as np
#import pyNastran

#root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'test', 'unit')

import pyNastran
from pyNastran.bdf.dev_vectorized.bdf import read_bdf
from pyNastran.bdf.dev_vectorized.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.utils.log import SimpleLogger

pkg_path = pyNastran.__path__[0]

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)


class TestMeshUtilsVectorized(unittest.TestCase):

    def test_renumber_01(self):
        log = SimpleLogger(level='warning')
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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
