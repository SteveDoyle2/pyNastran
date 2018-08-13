# coding: utf-8
"""
tests aero cards
"""
from __future__ import print_function
import unittest

import pyNastran
#from pyNastran.utils.log import SimpleLogger
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.test.test_bdf import run_bdf
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
        model.add_grid(1, [0., 0., 0.])
        pid = 100
        mid1 = 1000
        nsm = 0.0
        z1 = 0.
        z2 = 0.
        pconeax = model.add_pconeax(pid, mid1, t1=None, mid2=0, i=None, mid3=None,
                                    t2=None, nsm=nsm, z1=z1,
                                    z2=z2, phi=None,
                                    comment='pconeax')
        model.add_mat1(mid=mid1, E=3.0e7, G=None, nu=0.3)

        eid = 10
        rings = [2, 3]
        cconeax = model.add_cconeax(eid, pid, rings, comment='ccone')

        sid = 42
        ring = 6
        phi = 37.
        temperature = 420.
        tempax = model.add_tempax(sid, ring, phi, temperature, comment='tempax')

        R = 1.2
        z = 3.2
        ringax = model.add_ringax(ring, R, z, ps=None, comment='ringax')

        nid = 7
        pointax = model.add_pointax(nid, ring, phi, comment='pointax')

        nharmonics = 12
        axic = model.add_axic(nharmonics, comment='axic')
        ringax.raw_fields()
        pointax.raw_fields()
        cconeax.raw_fields()
        pconeax.raw_fields()
        axic.raw_fields()
        save_load_deck(model, run_convert=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
