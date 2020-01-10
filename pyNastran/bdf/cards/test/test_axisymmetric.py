# coding: utf-8
"""tests aero cards"""
import unittest

from cpylog import get_logger
from pyNastran.bdf.bdf import BDF
#from pyNastran.bdf.test.test_bdf import run_bdf
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestAxi(unittest.TestCase):
    """
    The cards are:
     * CCONEAX
     * AXIC
     * CQUADX

    """
    def test_cquadx(self):
        log = get_logger(level='warning')
        model = BDF(debug=False, log=log, mode='msc')
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [0., 0., 0.])
        model.add_grid(13, [0., 0., 0.])
        model.add_grid(14, [0., 0., 0.])
        nids = [11, 12, 13, 14,
                None, None, None, None,
                None]
        eid = 10
        pid = 20
        mid = 30
        model.add_cquadx(eid, pid, nids, theta_mcid=10., comment='cquadx_a')

        #PLPLANE or PAXSYMH
        model.add_cquadx(eid+1, pid, nids, theta_mcid=10, comment='cquadx_b')
        model.add_plplane(pid, mid, cid=0, stress_strain_output_location='GRID', comment='plplane')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_mathp(mid, a10=0., a01=0., d1=None, rho=0.,
                        av=0., tref=0., ge=0., na=1, nd=1,
                        a20=0., a11=0., a02=0., d2=0.,
                        a30=0., a21=0., a12=0., a03=0., d3=0.,
                        a40=0., a31=0., a22=0., a13=0., a04=0., d4=0.,
                        a50=0., a41=0., a32=0., a23=0., a14=0., a05=0., d5=0.,
                        tab1=None, tab2=None, tab3=None, tab4=None, tabd=None,
                        comment='mathp')
        #model.add_mathe(mid, model, bulk, rho, texp, mus, alphas, betas, mooney,
                        #sussbat, aboyce, comment='')
        sid = 1
        ring_id = 2
        hid = 0
        scale = 3.
        f_rtz = [0., 1., 2.]
        forceax = model.add_forceax(sid, ring_id, hid, scale, f_rtz, comment='forceax')

        pressure = 2.0
        rid1 = 2
        rid2 = 2
        presax = model.add_presax(sid, pressure, rid1, rid2, phi1=0., phi2=360., comment='prsax')
        forceax.raw_fields()
        presax.raw_fields()

        model.validate()
        model._verify_bdf()
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(model)

    def test_pconeax(self):
        """PCONEAX"""
        log = get_logger(level='warning')
        model = BDF(log=log)
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
        tempax.raw_fields()
        ringax.raw_fields()
        pointax.raw_fields()
        cconeax.raw_fields()
        pconeax.raw_fields()
        axic.raw_fields()
        save_load_deck(model, run_mass_properties=False, run_test_bdf=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
