"""tests remove_unused"""
import os
import unittest

import numpy as np


from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf, CaseControlDeck, PARAM
from pyNastran.bdf.mesh_utils.convert import convert, get_scale_factors
from cpylog import SimpleLogger

pkg_path = pyNastran.__path__[0]

np.set_printoptions(edgeitems=3, infstr='inf',
                    linewidth=75, nanstr='nan', precision=3,
                    suppress=True, threshold=1000, formatter=None)

log = SimpleLogger(level='error')
class TestRemoveUnused(unittest.TestCase):
    """various BDF cleanup tests"""

    def test_remove_bar(self):
        """removes unused data from the bar model"""
        model_path = os.path.join(pkg_path, '..', 'models', 'beam_modes')
        bdf_filename = os.path.join(model_path, 'beam_modes.dat')
        bdf_filename_out = os.path.join(model_path, 'beam_modes_temp.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)

        remove_unused(model)
        #os.remove(bdf_filename_out)

    def _test_remove_isat(self):
        """removes unused data from the isat model"""
        model_path = os.path.join(pkg_path, '..', 'models', 'isat')
        bdf_filename = os.path.join(model_path, 'ISat_Dploy_Sm.dat')
        bdf_filename_out = os.path.join(model_path, 'isat.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)

        remove_unused(model)
        model.write_bdf(bdf_filename_out)
        os.remove(bdf_filename_out)

    def test_remove_bwb(self):
        """removes unused data from the bwb model"""
        bdf_filename = os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_saero.bdf')
        bdf_filename_out = os.path.join(pkg_path, '..', 'models', 'bwb', 'bwb_modes.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)

        remove_unused(model)
        model.write_bdf(bdf_filename_out)
        os.remove(bdf_filename_out)

    def test_remove_sine(self):
        """removes unused data from the sine model"""
        model_path = os.path.join(pkg_path, '..', 'models', 'freq_sine')
        bdf_filename = os.path.join(model_path, 'good_sine.dat')
        bdf_filename_out = os.path.join(model_path, 'sine_modes.bdf')
        model = read_bdf(bdf_filename, log=log, validate=False)

        remove_unused(model)
        model.write_bdf(bdf_filename_out)
        os.remove(bdf_filename_out)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
