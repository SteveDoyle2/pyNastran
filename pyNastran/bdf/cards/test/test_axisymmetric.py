# coding: utf-8
"""
tests aero cards
"""
from __future__ import print_function
import os
import unittest
from six import StringIO
import numpy as np

import pyNastran
from pyNastran.utils.log import SimpleLogger
from pyNastran.bdf.bdf import BDF, CORD2R, BDFCard, SET1, GRID, read_bdf
from pyNastran.bdf.test.test_bdf import run_bdf
from pyNastran.bdf.cards.aero import (
    FLFACT, AEFACT, AEPARM, AERO, AEROS, AESTAT,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    AELIST, FLUTTER, TRIM, CSSCHD, MKAERO1, MKAERO2, GUST, AESURF, AESURFS,
    AELINK, DIVERG, AECOMP,
    SPLINE1, SPLINE2 #, SPLINE3, SPLINE4, SPLINE5
)
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestAxi(unittest.TestCase):
    """
    The cards are:
     * CCONEAX
     * AXIC
    """
    def test_pconeax(self):
        """PCONEAX"""
        model = BDF(debug=False)
        model.add_grid(1)
        pid = 10
        mid1 = 100
        nsm = 0.0
        z1 = 0.
        z2 = 0.
        pconeax = model.add_pconeax(pid, mid1, t1=None, mid2=0, i=None, mid3=None,
                                    t2=None, nsm=nsm, z1=z1,
                                    z2=z2, phi=None,
                                    comment='pconeax')
        model.add_mat1(mid=mid1, E=3.0e7, G=None, nu=0.3)

        nharmonics = 12
        axic = model.add_axic(nharmonics, comment='axic')
        pconeax.raw_fields()
        axic.raw_fields()
        save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

