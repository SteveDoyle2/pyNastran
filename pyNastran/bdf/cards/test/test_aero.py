from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.utils.log import SimpleLogger
from pyNastran.bdf.bdf import BDF, CORD2R, BDFCard
from pyNastran.bdf.cards.aero import (
    FLFACT, AEFACT, AEPARM, AERO, AEROS, CAERO1, CAERO2, CAERO3, PAERO1, PAERO2, PAERO3,
    AELIST, FLUTTER, TRIM, CSSCHD, MKAERO1, MKAERO2, GUST, AESURF,
)

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

comment_bad = 'this is a bad comment'
comment_good = 'this is a good comment\n'
class TestAero(unittest.TestCase):
    """
    The Aero cards are:
     * AEFACT
     * AELINK
     * AELIST
     * AEPARM
     * AESTAT
     * AESURF / AESURFS
     * AERO / AEROS
     * CSSCHD
     * CAERO1 / CAERO2 / CAERO3 / CAERO4 / CAERO5
     * FLFACT
     * FLUTTER
     * GUST
     * MKAERO1 / MKAERO2
     * PAERO1 / PAERO2 / PAERO3
     * SPLINE1 / SPLINE2 / SPLINE4 / SPLINE5
    """
    def test_aefact_1(self):
        """checks the AEFACT card"""
        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model = BDF(debug=None)
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', '98', '.3', '0.7', '1.0']
        model.add_card(data, data[0], comment_good, is_list=True)

        msg = 'this is a bad commentAEFACT        97      .3      .7      1.\n'
        aefact97 = model.aefacts[97]
        aefact98 = model.aefacts[98]
        self.assertTrue(all(aefact97.Di == [.3, .7, 1.0]))
        self.assertTrue(all(aefact98.Di == [.3, .7, 1.0]))

        out = aefact97.write_card(8, None)
        self.assertEqual(msg, out)

        msg = 'this is a good comment\nAEFACT        98      .3      .7      1.\n'
        out = aefact98.write_card(8, None)
        self.assertEqual(msg, out)

        #data = ['AEFACT', 99, .3, 0.7, 1.0, None, 'cat']
        #with self.assertRaises(SyntaxError):
            #model.add_card(data, data[0], comment_good, is_list=True)

        #data = ['AEFACT', 100, .3, 0.7, 1.0, 'cat']
        #with self.assertRaises(SyntaxError):
            #model.add_card(data, data[0], comment_good, is_list=True)

        #data = ['AEFACT', 101, .3, 0.7, 1.0, 2]
        #with self.assertRaises(SyntaxError):
            #model.add_card(data, data[0], comment_good, is_list=True)

        Di = [1., 2., 3.]
        aefact = AEFACT(200, Di, comment='')
        aefact.validate()
        aefact.write_card()

   # def test_aelink_1(self):
    def test_aelist_1(self):
        """checks the AELIST card"""
        model = BDF(debug=None)
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1201, 1202]
        aelist = AELIST(74, elements)
        aelist.validate()
        aelist.write_card()
        aelist75 = model.aelists[75]
        #print(aelist.elements)
        #print(elements)
        self.assertTrue(elements == aelist75.elements)

        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1108, 1202]
        data = ['AELIST', 76, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1108, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        aelist76 = model.aelists[76]
        #print(aelist76 .elements)
        #print(elements)
        self.assertFalse(elements == aelist76.elements)

        elements = list(set(elements))
        elements.sort()
        self.assertTrue(elements == aelist76.elements)

    def test_aeparm_1(self):
        """checks the AEPARM card"""
        aeparm = AEPARM.add_card(BDFCard(['AEPARM', 100, 'THRUST', 'lb']))
        aeparm = AEPARM(100, 'THRUST', 'lb')
        aeparm.validate()
        aeparm.write_card()

   # def test_aestat_1(self):
   # def test_aesurf_1(self):
   # def test_aesurfs_1(self):

    def test_aero_1(self):
        """checks the AERO card"""
        acsid = 0.
        velocity = None
        cref = 1.0
        rho_ref = 1.0
        aero = AERO(acsid, velocity, cref, rho_ref, sym_xz=0., sym_xy=0,
                    comment='aero card')
        aero.validate()
        aero.write_card()
        aero.raw_fields()

        acsid = None

        aero = AERO(acsid, velocity, cref, rho_ref, sym_xz=0., sym_xy=0,
                    comment='aero card')
        aero.validate()
        aero.write_card()
        aero.raw_fields()

    def test_aeros_1(self):
        """checks the AEROS card"""
        acsid = 0.
        #velocity = None
        cref = 1.0
        bref = 2.0
        sref = 100.
        acsid = 0
        rcsid = 0
        aeros = AEROS.add_card(BDFCard(['AERO', acsid, rcsid, cref, bref, sref]))
        aeros = AEROS(cref, bref, sref, acsid, rcsid, sym_xz=0, sym_xy=0,
                      comment='aeros card')
        aeros.validate()
        aeros.write_card()
        aeros.raw_fields()

        acsid = None
        rcsid = None
        sym_xz = None
        sym_xy = None
        aeros = AEROS(cref, bref, sref, acsid, rcsid, sym_xz=sym_xz, sym_xy=sym_xy,
                      comment='aeros card')
        aeros.validate()
        aeros.write_card()
        aeros.raw_fields()


    def test_caero1_1(self):
        """checks the CAERO1/PAERO1/AEROS/AEFACT card"""
        eid = 1
        pid = 10
        cp = 4
        nspan = None
        lspan = 3
        nchord = None
        lchord = 4
        igid = 0
        p1 = [0., 0., 0.]
        x12 = 5.
        p4 = [2., 3., 4.]
        x43 = 1.

        model = BDF()
        caero = CAERO1.add_card(BDFCard(['CAERO1', eid, pid, cp, nspan, nchord, lspan, lchord,
                                         igid, ] + p1 + [x12] + p4 + [x43]))
        caero.validate()
        caero = CAERO1.add_card(BDFCard(['CAERO1', eid, pid, None, nspan, nchord, lspan, lchord,
                                         igid, ] + p1 + [x12] + p4 + [x43]))
        caero.validate()
        caero = CAERO1(eid, pid, cp, nspan, lspan, nchord, lchord, igid, p1,
                       x12, p4, x43, comment='caero1')
        caero.raw_fields()
        caero.validate()
        caero.write_card()
        model.caeros[eid] = caero

        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 1., 0.]
        p4 = [0.1, 1., 0.]
        span = 5
        chord = 10
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, cp, span, chord, igid, p1, p2, p3, p4,
                                 spanwise='y', comment='')
        caeroq.validate()

        span = 0.1
        chord = 0.05
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, cp, span, chord, igid, p1, p2, p3, p4,
                                 spanwise='y', comment='')
        caeroq.validate()


        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 0., 1.]
        p4 = [0.1, 0., 1.]
        span = 0.1
        chord = 0.05
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, cp, span, chord, igid, p1, p2, p3, p4,
                                 spanwise='z', comment='')
        caeroq.validate()

        paero = PAERO1(pid, Bi=None, comment='')
        paero.validate()
        paero.write_card()
        model.paeros[pid] = paero

        coord = CORD2R(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                       comment='')
        coord.validate()
        model.coords[cp] = coord

        acsid = 0.
        #velocity = None
        cref = 1.0
        bref = 2.0
        sref = 100.
        acsid = 0
        rcsid = 0
        aeros = AEROS(cref, bref, sref, acsid, rcsid, sym_xz=0, sym_xy=0,
                      comment='')
        aeros.validate()
        aeros.write_card()
        model.aeros = aeros

        aefact = AEFACT(lspan, [0., 1., 2., 3., 4., 5.])
        aefact.validate()
        model.aefacts[lspan] = aefact

        aefact = AEFACT(lchord, [2., 3., 4., 5., 6., 7.])
        aefact.validate()
        model.aefacts[lchord] = aefact

        paero.cross_reference(model)
        caero.cross_reference(model)
        caero.get_npanel_points_elements()
        caero.get_points()
        caero.panel_points_elements()

        caero.write_card()
        model.uncross_reference()
        model.cross_reference()
        model.uncross_reference()
        #model.safe_cross_reference()
        caero.safe_cross_reference(model)
        caero.panel_points_elements()
        caero.raw_fields()

    def test_caero2_1(self):
        """checks the CAERO2/PAERO2/AERO/AEFACT card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        eid = 1
        pid = 10
        cp = 4
        nsb = 0
        nint = 0

        lsb = 3
        lint = 6
        igid = 0
        p1 = [0., 1., 2.]
        x12 = 10.
        caero = CAERO2.add_card(BDFCard(['CAERO2', eid, pid, cp, nsb, nint,
                                         lsb, lint, igid, ] + p1 + [x12]))
        caero = CAERO2(eid, pid, cp, nsb, nint, lsb, lint, igid, p1, x12,
                       comment='this is a caero')
        caero.validate()
        caero.write_card()

        aefact = AEFACT.add_card(BDFCard(['AEFACT', lint, 0., 1., 2., 3., 4., 5.]))
        aefact = AEFACT(lint, [0., 1., 2., 3., 4., 5.])
        aefact.validate()
        aefact.write_card()
        model.aefacts[lint] = aefact

        orient = 'Z'
        width = 10.
        AR = 2.
        lrsb = 0
        lrib = 3
        lth1 = 0
        lth2 = 0
        thi = [0]
        thn = [0]
        paero = PAERO2.add_card(BDFCard(['PAERO2', pid, orient, width, AR,
                                         lrsb, lrib, lth1, lth2] + thi + thn),
                                comment='paero')
        paero = PAERO2(pid, orient, width, AR, lrsb, lrib, lth1, lth2, thi, thn)
        paero.validate()
        paero.write_card()
        model.paeros[pid] = paero

        coord = CORD2R.add_card(BDFCard(['CORD2R', cp, 0,
                                         0., 0., 0.,
                                         0., 0., 1.,
                                         1., 0., 0.]))
        coord = CORD2R(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                       comment='')
        coord.validate()
        model.coords[cp] = coord

        aefact = AEFACT(lrib, [0., 1., 2., 3., 4., 5.])
        aefact.validate()
        model.aefacts[lrib] = aefact

        acsid = 0
        velocity = None
        cref = 1.0
        rho_ref = 1.0

        aero = AERO.add_card(BDFCard(['AERO', acsid, velocity, cref, rho_ref]))
        aero = AERO(acsid, velocity, cref, rho_ref,
                    comment='')
        aero.validate()
        aero.write_card()
        model.aero = aero

        paero.cross_reference(model)
        caero.cross_reference(model)
        paero.raw_fields()
        caero.raw_fields()
        caero.uncross_reference()
        caero.raw_fields()
        caero.cross_reference(model)
        caero.get_points_elements_3d()

        caero.get_points()
        #caero.get_points_elements_3d()
        xyz, elems = caero.get_points_elements_3d()
        model.uncross_reference()
        model.safe_cross_reference()
        model.uncross_reference()
        model.write_bdf('aero.temp')
        os.remove('aero.temp')

        model.cross_reference()
        model.write_bdf('aero.temp')
        os.remove('aero.temp')

        nsb = 4
        nint = 2
        lsb = None
        lint = None
        caero2 = CAERO2(eid, pid, cp, nsb, nint, lsb, lint, igid, p1, x12,
                       comment='this is a caero')
        caero2.validate()
        caero2.cross_reference(model)
        caero2.write_card()

    def test_caero3_1(self):
        """checks the CAERO3/PAERO3"""
        eid = 100
        pid = 200
        cp = 4
        list_w = 5
        list_c1 = 6
        list_c2 = 7
        p1 = [0., 0., 0.]
        x12 = 10.
        p4 = [5., 10., 0.]
        x43 = 3.

        nbox = 10
        ncontrol_surfaces = 0
        x = None
        y = None

        model = BDF()
        coord = CORD2R.add_card(BDFCard(['CORD2R', cp, 0,
                                         0., 0., 0.,
                                         0., 0., 1.,
                                         1., 0., 0.]))
        coord = CORD2R(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                       comment='')
        coord.validate()
        model.coords[cp] = coord

        paero = PAERO3(pid, nbox, ncontrol_surfaces, x, y)
        model.paeros[pid] = paero

        caero = CAERO3(eid, pid, cp, list_w, list_c1, list_c2, p1, x12, p4,
                       x43, comment='caero3')
        model.caeros[pid] = caero

        caero.write_card()
        caero.cross_reference(model)
        caero.write_card()
        caero.uncross_reference()
        caero.write_card()
   # def test_caero4_1(self):
   # def test_caero5_1(self):

   # def test_paero1_1(self):
   # def test_paero2_1(self):
   # def test_paero3_1(self):
   # def test_paero4_1(self):
   # def test_paero5_1(self):

   # def test_spline1_1(self):
   # def test_spline2_1(self):
   # def test_spline3_1(self):
   # def test_spline4_1(self):
   # def test_spline5_1(self):

    def test_flutter(self):
        """checks the FLUTTER/FLFACT cards"""
        model = BDF()
        sid = 75
        method = 'PKNL'
        idensity = 76
        imach = 77
        ireduced_freq_velocity = 78

        flutter1 = FLUTTER(sid, method, idensity, imach, ireduced_freq_velocity)
        flutter2 = FLUTTER.add_card(BDFCard(['FLUTTER', sid, method, idensity, imach,
                                             ireduced_freq_velocity]), comment='flutter card')
        flutter1.validate()
        flutter1.write_card()
        flutter2.validate()
        flutter2.write_card()
        model.flutters[75] = flutter1

        densities = np.linspace(0., 1.)
        density = FLFACT(idensity, densities)
        model.flfacts[idensity] = density

        machs = np.linspace(0.7, 0.8)
        mach = FLFACT(imach, machs)
        mach = FLFACT.add_card(BDFCard(['FLFACT', imach] + list(machs)), comment='flfact card')
        mach.write_card(size=16)
        model.flfacts[imach] = mach

        velocities = np.linspace(3., 4.)
        velocity = FLFACT(ireduced_freq_velocity, velocities)
        velocity.validate()
        velocity.write_card()
        assert velocity.min() == 3., velocities
        assert velocity.max() == 4., velocities
        model.flfacts[ireduced_freq_velocity] = velocity

        model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()

    def test_mkaero1(self):
        """checks the MKAERO1 card"""
        machs = [0.5, 0.75]
        reduced_freqs = [0.1, 0.2, 0.3, 0.4]
        mkaero = MKAERO1(machs, reduced_freqs, comment='mkaero')
        mkaero.validate()
        mkaero.write_card()
        mkaero = MKAERO1.add_card(BDFCard(
            ['MKAERO', 0.5, 0.75, None, None, None, None, None, None,
             0.1, 0.2, 0.3, 0.4],
        ))

        machs = [0.5, 0.75]
        reduced_freqs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
        mkaero = MKAERO1(machs, reduced_freqs, comment='mkaero')
        mkaero.validate()
        mkaero.write_card()

    def test_mkaero2(self):
        """checks the MKAERO2 card"""
        machs = [0.5, 0.75, 0.8]
        reduced_freqs = [0.1, 0.2, 0.3]
        mkaero = MKAERO2(machs, reduced_freqs, comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()

        machs = [0.5, 0.75]
        reduced_freqs = [0.1, 0.2]
        mkaero = MKAERO2(machs, reduced_freqs, comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()

        mkaero = MKAERO2.add_card(BDFCard(['MKAERO2'] + machs + reduced_freqs), comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()

    def test_gust(self):
        sid = 100
        dload = 200
        wg = 50.
        x0 = 3.
        V = 42.
        gust = GUST(sid, dload, wg, x0, V, comment='gust load')
        gust.validate()
        gust.write_card()

        gust2 = GUST.add_card(BDFCard(['GUST', sid, dload, wg, x0, V]), comment='gust load')
        gust2.validate()
        gust2.write_card()

    def test_cssch(self):
        sid = 10
        aesid = 0
        lalpha = None
        lmach = None
        lschd = None
        csshcd = CSSCHD(sid, aesid, lalpha, lmach, lschd, comment='cssch card')
        csshcd.validate()
        csshcd.write_card()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
