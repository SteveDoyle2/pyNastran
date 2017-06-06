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
from pyNastran.bdf.cards.aero import (
    FLFACT, AEFACT, AEPARM, AERO, AEROS, AESTAT,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    AELIST, FLUTTER, TRIM, CSSCHD, MKAERO1, MKAERO2, GUST, AESURF, AESURFS,
    AELINK, DIVERG, AECOMP,
    SPLINE1, SPLINE2 #, SPLINE3, SPLINE4, SPLINE5
)

ROOTPATH = pyNastran.__path__[0]
#test_path = os.path.join(ROOTPATH, 'bdf', 'cards', 'test')

COMMENT_BAD = 'this is a bad comment'
COMMENT_GOOD = 'this is a good comment\n'
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

    def test_aestat_1(self):
        bdf = BDF(debug=False)
        lines = ['AESTAT  502     PITCH']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = AESTAT.add_card(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_aecomp_1(self):
        """checks the AECOMP card"""
        #sid = 10
        #aesid = 0
        #lalpha = None
        #lmach = None
        #lschd = None

        #sid = 5
        #aesid = 50
        #lalpha = 12
        #lmach = 15
        name = 'WING'
        list_type = 'AELIST' # or SET1, CAEROx
        aelist_ids = [75, 76]

        card = ['AECOMP', name, list_type] + aelist_ids
        bdf_card = BDFCard(card, has_none=True)
        aecomp1 = AECOMP.add_card(bdf_card, comment='aecomp card')
        aecomp1.validate()
        aecomp1.write_card()

        #label = 'ELEV'
        #cid1 = 0
        #alid1 = 37
        #aesurf = AESURF(aesid, label, cid1, alid1)

        #aefact_sid = alid1
        #Di = [0., 0.5, 1.]
        #aefact_elev = AEFACT(aefact_sid, Di)

        #aefact_sid = lalpha
        #Di = [0., 5., 10.]
        #aefact_alpha = AEFACT(aefact_sid, Di)

        #aefact_sid = lmach
        #Di = [0., 0.7, 0.8]
        #aefact_mach = AEFACT(aefact_sid, Di)

        #aefact_sid = lschd
        #Di = [0., 15., 30., 45.]
        #aefact_delta = AEFACT(aefact_sid, Di)

        model = BDF()
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AELIST', 76, 2000, 'THRU', 2010]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        #model.add_aesurf(aesurf)
        #model.add_aefact(aefact_elev)
        #model.add_aefact(aefact_alpha)
        #model.add_aefact(aefact_mach)
        #model.add_aefact(aefact_delta)

        aecomp1.cross_reference(model)
        aecomp1.write_card()
        aecomp1.uncross_reference()
        aecomp1.write_card()

        #-----------
        aecomp2 = AECOMP(name, list_type, aelist_ids, comment='cssch card')
        aecomp2.validate()
        aecomp2.write_card()

        list_type = 'INVALID'
        aecomp3 = AECOMP(name, list_type, aelist_ids, comment='cssch card')
        with self.assertRaises(RuntimeError):
            aecomp3.validate()

    def test_aefact_1(self):
        """checks the AEFACT card"""
        data = ['AEFACT', 97, .3, 0.7, 1.0]
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AEFACT', '98', '.3', '0.7', '1.0']
        model.add_card(data, data[0], COMMENT_GOOD, is_list=True)

        msg = '$this is a bad comment\nAEFACT        97      .3      .7      1.\n'
        aefact97 = model.aefacts[97]
        aefact98 = model.aefacts[98]
        self.assertTrue(all(aefact97.Di == [.3, .7, 1.0]))
        self.assertTrue(all(aefact98.Di == [.3, .7, 1.0]))

        out = aefact97.write_card(8, None)
        self.assertEqual(msg, out)

        msg = '$this is a good comment\nAEFACT        98      .3      .7      1.\n'
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
        #model = BDF()
        #aefact.cross_reference(model)
        #aefact.write_card()
        #aefact.uncross_reference()
        #aefact.write_card()


    def test_aelink_1(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        id = 10
        label = 'CS'
        independent_labels = ['A', 'B', 'C']
        Cis = [1.0, 2.0]
        aelink = AELINK(id, label, independent_labels, Cis, comment='')
        with self.assertRaises(RuntimeError):
            aelink.validate()
        str(aelink)
        aelink.write_card()

        card = ['AELINK', id, label, independent_labels[0], Cis[0],
                independent_labels[1], Cis[1], independent_labels[2]]
        with self.assertRaises(AssertionError):
            model.add_card(card, 'AELINK')

        card = ['AELINK', id, label, independent_labels[0], Cis[0],
                independent_labels[1], Cis[1]]
        model.add_card(card, 'AELINK', comment='cat')



    def test_aelist_1(self):
        """checks the AELIST card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)
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
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)
        aelist76 = model.aelists[76]
        #print(aelist76 .elements)
        #print(elements)
        self.assertFalse(elements == aelist76.elements)

        elements = list(set(elements))
        elements.sort()
        self.assertTrue(elements == aelist76.elements)

        elements = [1000, 1000, 1000, 2000, 1000, 2000]
        aelist = AELIST(75, elements)
        aelist.clean_ids()
        str(aelist.write_card())

    def test_aeparm_1(self):
        """checks the AEPARM card"""
        aeparm = AEPARM.add_card(BDFCard(['AEPARM', 100, 'THRUST', 'lb']),
                                 comment='aeparm_comment')
        aeparm = AEPARM(100, 'THRUST', 'lb', comment='aeparm_comment')
        aeparm.validate()
        aeparm.write_card()

   # def test_aestat_1(self):
   # def test_aesurf_1(self):
    def test_aesurfs_1(self):
        """checks the AESURFS cards"""
        aesid = 6001
        label = 'ELEV'
        list1 = 6002
        list2 = 6003
        card = ['AESURFS', aesid, label, None, list1, None, list2]
        bdf_card = BDFCard(card, has_none=True)
        model = BDF()
        model.add_card(bdf_card, 'AESURFS', comment='aesurfs',
                       is_list=True, has_none=True)
        aesurfs = AESURFS(aesid, label, list1, list2, comment='aesurfs')
        str(aesurfs)
        aesurfs.write_card()

    def test_aero_1(self):
        """checks the AERO card"""
        acsid = 0.
        velocity = None
        cref = 1.0
        rho_ref = 1.0
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=0, sym_xy=0,
                    comment='aero card')
        with self.assertRaises(TypeError):
            aero.validate()

        acsid = 0
        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=0., sym_xy=0,
                    comment='aero card')
        with self.assertRaises(TypeError):
            aero.validate()

        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=0, sym_xy=0.,
                    comment='aero card')
        with self.assertRaises(TypeError):
            aero.validate()

        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=0, sym_xy=0.,
                    comment='aero card')
        with self.assertRaises(TypeError):
            aero.validate()

        aero = AERO(velocity, cref, rho_ref, acsid=None, sym_xz=0, sym_xy=0,
                    comment='aero card')
        aero.validate()
        aero.write_card()
        aero.raw_fields()

        model = BDF()
        aero.cross_reference(model)
        aero.write_card()
        aero.raw_fields()

        aero.uncross_reference()
        aero.write_card()
        aero.raw_fields()

    def test_aeros_1(self):
        """checks the AEROS card"""
        #acsid = 0.
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

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero = CAERO1.add_card(BDFCard(['CAERO1', eid, pid, cp, nspan, nchord, lspan, lchord,
                                         igid, ] + p1 + [x12] + p4 + [x43]))
        caero.validate()
        caero = CAERO1.add_card(BDFCard(['CAERO1', eid, pid, None, nspan, nchord, lspan, lchord,
                                         igid, ] + p1 + [x12] + p4 + [x43]))
        caero.validate()
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        caero.raw_fields()
        caero.validate()
        caero.write_card()
        model.caeros[eid] = caero

        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 1., 0.]
        p4 = [0.1, 1., 0.]
        nspan = 5
        nchord = 10
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, nspan, nchord, igid, p1, p2, p3, p4,
                                 cp=cp, spanwise='y', comment='')
        caeroq.validate()

        span = 0.1
        chord = 0.05
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, span, chord, igid, p1, p2, p3, p4,
                                 cp=cp, spanwise='y', comment='')
        caeroq.validate()


        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 0., 1.]
        p4 = [0.1, 0., 1.]
        span = 0.1
        chord = 0.05
        igid = -1
        caeroq = CAERO1.add_quad(eid, pid, span, chord, igid, p1, p2, p3, p4,
                                 cp=cp, spanwise='z', comment='')
        caeroq.validate()

        paero = PAERO1(pid, Bi=None, comment='')
        paero.validate()
        paero.write_card()
        model.paeros[pid] = paero

        coord = CORD2R(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                       comment='')
        coord.validate()
        model.coords[cp] = coord

        #acsid = 0.
        #velocity = None
        cref = 1.0
        bref = 2.0
        sref = 100.
        acsid = 0
        rcsid = 0
        aeros = model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0,
                       sym_xy=0, comment='aeros')
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
        min_max_eid = caero.min_max_eid
        self.assertEqual(min_max_eid, [1, 26])
        #print('min_eid, max_eid', min_eid, max_eid)

        points = [
            [0., 0., 0.], # p1
            [10., 0., 0.],
            [10., 20., 0.],
            [5., 20., 0.],
        ]
        caero.set_points(points)
        caero.get_points()
        str(caero.write_card())

        nspan = None
        lspan = None
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=None,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(ValueError):
            caero.validate()

        nspan = 5
        lspan = 5
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(ValueError):
            caero.validate()

        nspan = 5
        nchord = None
        lchord = None
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(ValueError):
            caero.validate()

        nchord = 10
        lchord = 10
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(ValueError):
            caero.validate()

        lspan = None
        lchord = None
        nspan = 10
        nchord = 10
        p1 = [0., 0., 0., 0.]
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(AssertionError):
            caero.validate()

        p1 = [0., 0., 0.]
        p4 = [1., 2., 3., 4.]
        caero = CAERO1(eid, pid, igid, p1, x12, p4, x43, cp=cp,
                       nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                       comment='caero1')
        with self.assertRaises(AssertionError):
            caero.validate()


    def test_spline1(self):
        """checks the SPLINE1 card"""
        eid = 1
        caero_id = 1
        box1 = 1
        box2 = 10
        setg = 1
        spline = SPLINE1(eid, caero_id, box1, box2, setg, dz=0., method='IPS',
                         usage='BOTH', nelements=10,
                         melements=10, comment='$ spline1')
        spline.validate()
        spline.write_card(size=8, is_double=False)

    def test_spline2(self):
        """checks the SPLINE2 card"""
        #| SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
        #|         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
        #+---------+------+-------+-------+-------+------+----+------+-----+
        #| SPLINE2 |   5  |   8   |  12   | 24    | 60   | 0. | 1.0  |  3  |
        #|         |  1.  |       |       |       |      |    |      |     |

        cid = 3
        origin = [0., 0., 0.]
        xaxis = [1., 0., 0.]
        xyplane = [0., 1., 0.]
        coord = CORD2R.add_axes(cid, rid=0, origin=origin,
                                xaxis=xaxis, yaxis=None, zaxis=None,
                                xyplane=xyplane, yzplane=None, xzplane=None,
                                comment='comment')
        eid = 8
        pid = 10
        cp = 0
        nsb = 4
        nint = 2
        lsb = None
        lint = None
        p1 = [0., 0., 0.]
        x12 = 42.
        igid = None
        caero2 = CAERO2(eid, pid, igid, p1, x12,
                        cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                        comment='this is a caero')
        #caero = CAERO2(eid, pid, cp, nsb, nint, lsb, lint, igid, p1, x12)

        sid = 60
        ids = [7, 13]
        set_obj = SET1(sid, ids, is_skin=False, comment='set card')
        grid7 = GRID(nid=7, cp=0, xyz=[7., 0., 0.], cd=0, ps='', seid=0, comment='')
        grid13 = GRID(nid=13, cp=0, xyz=[13., 0., 0.], cd=0, ps='', seid=0, comment='')

        model = BDF(log=None)
        model._add_coord_object(coord)
        model._add_caero_object(caero2)
        model._add_set_object(set_obj)
        model._add_node_object(grid7)
        model._add_node_object(grid13)

        eid = 5
        caero = 8
        id1 = 12
        id2 = 24
        setg = 60
        dz = 0.
        dtor = 1.0
        cid = 3
        dthx = 1.
        dthy = None
        usage = None
        card = ['SPLINE2', eid, caero, id1, id2, setg, dz, dtor, cid,
                dthx, dthy, None, usage]

        bdf_card = BDFCard(card, has_none=True)
        spline_a = SPLINE2.add_card(bdf_card, comment='spline2_a')
        spline_a.write_card()

        spline_b = SPLINE2(eid, caero, id1, id2, setg, dz, dtor, cid, dthx,
                           dthy, usage, comment='spline2_b')
        spline_b.validate()
        spline_b.write_card()
        spline_b.cross_reference(model)
        spline_b.write_card()

        #model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()


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
        caero2a = CAERO2.add_card(BDFCard(['CAERO2', eid, pid, cp, nsb, nint,
                                           lsb, lint, igid, ] + p1 + [x12]))

        #---------------
        # nsb=lsb=None=0
        caero2b = CAERO2(eid, pid, igid, p1, x12,
                         cp=cp, nsb=None, nint=None, lsb=None, lint=6,
                         comment='this is a caero')
        with self.assertRaises(ValueError):
            caero2b.validate()

        # nint=lint=None=0
        caero2c = CAERO2(eid, pid, igid, p1, x12,
                         cp=cp, nsb=3, nint=None, lsb=3, lint=None,
                         comment='this is a caero')
        with self.assertRaises(ValueError):
            caero2c.validate()

        # they're all bad?
        caero2e = CAERO2(eid, pid, igid, p1, x12,
                         cp=cp, nsb=0, nint=0, lsb=0, lint=0,
                         comment='this is a caero')
        with self.assertRaises(ValueError):
            caero2e.validate()

        #---------------
        caero2f = model.add_caero2(eid, pid, igid, p1, x12, cp=4, nsb=0, nint=0,
                                   lsb=3, lint=6, comment='this is a caero')

        eid = 200
        caero2g = model.add_caero2(eid, pid, igid, p1, x12, cp=4, nsb=10, nint=7,
                                   lsb=0, lint=0, comment='this is a caero')
        caero2f.validate()
        caero2g.validate()
        #str(caero2f.write_card())

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
        paero2a = PAERO2.add_card(BDFCard(['PAERO2', pid, orient, width, AR,
                                           lrsb, lrib, lth1, lth2] + thi + thn),
                                  comment='paero2')
        paero2a.validate()
        paero2b = model.add_paero2(pid, orient, width, AR, thi, thn,
                                   lrsb=0, lrib=3, lth1=0, lth2=0, comment='paero2')

        pid = 42
        paero2c = model.add_paero2(pid, orient, width, AR, thi, thn,
                                   lrsb=None, lrib=None, lth1=None, lth2=None, comment='')
        paero2b.validate()
        paero2c.validate()
        paero2b.write_card()
        #model.paeros[pid] = paero

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
        aero = AERO(velocity, cref, rho_ref, acsid=acsid,
                    comment='aero')
        aero.validate()
        aero.write_card()
        model.aero = aero

        model.cross_reference()
        model.write_bdf('aero.temp1')

        paero2b.raw_fields()
        caero2f.raw_fields()
        model.uncross_reference()
        model.write_bdf('aero.temp2')

        model.cross_reference()
        model.write_bdf('aero.temp3')

        caero2f.raw_fields()
        caero2f.get_points_elements_3d()
        caero2f.get_points()
        xyz, elems = caero2f.get_points_elements_3d()


        caero2g.get_points()
        caero2g.get_points_elements_3d()
        xyz, elems = caero2g.get_points_elements_3d()

        model.uncross_reference()
        model.safe_cross_reference()
        model.uncross_reference()
        model.write_bdf('aero.temp4')

        model.cross_reference()
        model.write_bdf('aero.temp5')
        os.remove('aero.temp1')
        os.remove('aero.temp2')
        os.remove('aero.temp3')
        os.remove('aero.temp4')
        os.remove('aero.temp5')

        nsb = 4
        nint = 2
        lsb = None
        lint = None
        caero2 = CAERO2(eid, pid, igid, p1, x12,
                        cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                        comment='this is a caero')
        caero2.validate()
        caero2.cross_reference(model)
        caero2.write_card()

        #model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()

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
        x = []
        y = []

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        coord = CORD2R.add_card(BDFCard(['CORD2R', cp, 0,
                                         0., 0., 0.,
                                         0., 0., 1.,
                                         1., 0., 0.]))
        model.add_cord2r(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                         comment='cord2r')
        coord.validate()
        model.coords[cp] = coord

        paero3 = model.add_paero3(pid, nbox, ncontrol_surfaces, x, y,
                                  comment='paero3')
        paero3.validate()
        paero3.raw_fields()

        card = ['CAERO3', 2000, 20001, 0, 22, 33, None, None, None,
                1.0, 0.0, 0., 100., 17., 130., 0., 100.]
        bdf_card = BDFCard(card, has_none=True)
        caero3a = CAERO3.add_card(bdf_card, comment='msg')
        caero3a.validate()
        caero3a.write_card()
        caero3a.raw_fields()

        caero3b = model.add_caero3(eid, pid, list_w,
                                   p1, x12, p4, x43,
                                   cp=cp, list_c1=list_c1, list_c2=list_c2,
                                   comment='caero3')
        caero3b.validate()

        aefact_sid = list_w
        Di = [0., 0.5, 1.]
        model.add_aefact(aefact_sid, Di, comment='aefact')

        aefact_sid = list_c1
        model.add_aefact(aefact_sid, Di, comment='aefact2')

        aefact_sid = list_c2
        model.add_aefact(aefact_sid, Di, comment='aefact2')

        velocity = 100.
        cref = 1.0
        rho_ref = 1.0
        model.add_aero(velocity, cref, rho_ref)
        model.validate()

        caero3b.write_card()
        caero3b.cross_reference(model)
        caero3b.write_card()
        caero3a.raw_fields()
        caero3b.uncross_reference()
        caero3b.write_card()
        caero3a.raw_fields()
        caero3b.safe_cross_reference(model)

        caero3b.get_npanel_points_elements()
        caero3b.get_points()
        caero3b.panel_points_elements()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)
        model2 = read_bdf(bdf_filename, punch=True, debug=False)


    def test_caero4_1(self):
        """checks the CAERO4/PAERO4"""
        model = BDF(debug=False)
        pid = 1001
        docs = []
        caocs = []
        gapocs = []
        paero4 = model.add_paero4(pid, docs, caocs, gapocs,
                                  cla=0, lcla=0, circ=0, lcirc=0,
                                  comment='paero4')
        paero4.validate()
        paero4.raw_fields()

        x1 = 0.
        y1 = 0.
        z1 = 0.
        x12 = 100.
        x4 = 50.
        y4 = 0.
        z4 = 0.
        x43 = 10.

        eid = 1000
        nspan = 4  # number of stations
        lspan = 0  # AEFACT
        cp = 0
        card = ['CAERO4', eid, pid, cp, nspan, lspan, None, None, None,
                x1, y1, z1, x12, x4, y4, z4, x43]

        bdf_card = BDFCard(card, has_none=True)
        caero4a = CAERO4.add_card(bdf_card, comment='msg')
        caero4a.validate()
        caero4a.write_card()
        caero4a.raw_fields()

        p1 = [x1, y1, z1]
        p4 = [x4, y4, z4]
        caero4b = model.add_caero4(eid, pid, p1, x12, p4, x43,
                                   cp=cp, nspan=nspan, lspan=lspan,
                                   comment='caero4b')
        caero4b.validate()
        caero4b.write_card()
        caero4b.raw_fields()

        caero4b.cross_reference(model)
        caero4b.write_card()
        caero4b.raw_fields()
        p1, p2, p3, p4 = caero4b.get_points()

        caero4c = CAERO4(eid, pid, p1, x12, p4, x43,
                         cp=0,nspan=0, lspan=0,
                         comment='caero4c')
        with self.assertRaises(RuntimeError):
            caero4c.validate()


        #model.cross_reference()
        model.uncross_reference()
        #model.safe_cross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)
        model2 = read_bdf(bdf_filename, xref=False, punch=True, debug=False)
        model.safe_cross_reference()

    def test_caero5_1(self):
        """checks the CAERO5/PAERO5"""
        model = BDF(debug=False)
        pid = 6001
        caoci = [0., 0.5, 1.0]
        paero5 = model.add_paero5(pid, caoci,
                                  nalpha=0, lalpha=0, nxis=0, lxis=0,
                                 ntaus=0, ltaus=0, comment='paero5')
        paero5.validate()

        #| PAERO5 | PID   | NALPHA | LALPHA | NXIS    | LXIS  | NTAUS | LTAUS |
        #+--------+-------+--------+--------+---------+-------+-------+-------+
        #|        | CAOC1 | CAOC2  | CAOC3  | CAOC4   | CAOC5 |       |       |
        nalpha = 0
        lalpha = 0
        nxis = 0
        lxis = 0
        ntaus = 0
        ltaus = 0
        card = ['PAERO5', pid, nalpha, lalpha, nxis, lxis, ntaus, ltaus, ] + caoci

        model = BDF(debug=False)
        model.add_card(card, card[0], comment='', is_list=True,
                       has_none=True)
        paero5 = model.paeros[pid]
        paero5.raw_fields()

        model = BDF(debug=None)
        paero5 = model.add_paero5(pid, caoci, nalpha=0, lalpha=0, nxis=0, lxis=0,
                                  ntaus=0, ltaus=0, comment='paero5')
        paero5.validate()
        eid = 6000
        x1 = 0.
        y1 = 0.
        z1 = 0.
        x12 = 1.
        x4 = 0.2
        y4 = 1.
        z4 = 0.
        x43 = 0.8
        p1 = [x1, y1, z1]
        p4 = [x4, y4, z4]
        caero5 = model.add_caero5(eid, pid, p1, x12, p4, x43,
                                  cp=0, nspan=5, lspan=0,
                                  ntheory=0, nthick=0,
                                  comment='msg')
        model.validate()

        lxis = 43
        paero5.lxis = lxis
        aefact_sid = lxis
        Di = [0., 0.5, 1.]
        aefact = model.add_aefact(aefact_sid, Di, comment='aefact')

        ltaus = 44
        paero5.ltaus = ltaus
        aefact_sid = ltaus
        Di = [0., 0.5, 1.]
        aefact = model.add_aefact(aefact_sid, Di, comment='aefact2')

        #caero5.cross_reference(model)
        model.cross_reference()
        npoints, nelements = caero5.get_npanel_points_elements()
        points, elements = caero5.panel_points_elements()
        caero5.write_card()
        #caero5.raw_fields()

        model.uncross_reference()
        caero5.write_card()
        model.cross_reference()
        model.uncross_reference()
        #model.safe_cross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)

        read_bdf(bdf_filename, xref=False, punch=True, debug=False)
        model.safe_cross_reference()


        #caero5.raw_fields()


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

    def test_aesurf_1(self):
        """checks the AESURF/AELIST cards"""
        aesid = 10
        label = 'FLAP'
        cid1 = 0
        aelist_id1 = 10
        cid2 = None
        alid2 = None
        aesurf1 = AESURF(aesid, label, cid1, aelist_id1, cid2, alid2,
                         #eff, ldw,
                         #crefc, crefs, pllim, pulim,
                         #hmllim, hmulim, tqllim, tqulim,
                         comment='aesurf comment')
        aesurf2 = AESURF.add_card(BDFCard(
            [
                'AESURF', aesid, label, cid1, aelist_id1, cid2, alid2,
                #eff, ldw,
                #crefc, crefs, pllim, pulim,
                #hmllim, hmulim, tqllim, tqulim,
            ]), comment='aesurf comment')
        #assert aesurf1 == aesurf2

        cid2 = 1
        coord = CORD2R(cid2, rid=0, origin=[0., 0., 0.],
                       zaxis=[1., 0., 0.], xzplane=[0., 1., 1.], comment='')

        aelist_id1 = 10
        aelist_id2 = 20
        aesurf2 = AESURF.add_card(BDFCard(
            [
                'AESURF', aesid, label, cid1, aelist_id1, cid2, aelist_id2,
                #eff, ldw,
                #crefc, crefs, pllim, pulim,
                #hmllim, hmulim, tqllim, tqulim,
            ]), comment='aesurf comment')

        aesurf1.validate()
        aesurf2.validate()
        log = SimpleLogger(level='warning')

        model = BDF(debug=False)
        model._add_coord_object(coord)
        model._add_aesurf_object(aesurf1)

        elements = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        aelist = AELIST(aelist_id1, elements)
        model._add_aelist_object(aelist)

        elements = [11, 22, 33, 44, 55, 66, 77, 88, 99]
        aelist = AELIST(aelist_id2, elements)
        model._add_aelist_object(aelist)

        aesurf1.cross_reference(model)
        aesurf1.write_card()
        aesurf1.raw_fields()
        aesurf1.uncross_reference()
        aesurf1.write_card()
        aesurf1.cross_reference(model)
        aesurf1.raw_fields()

        aesurf2.cross_reference(model)
        aesurf2.write_card()
        aesurf2.raw_fields()
        aesurf2.uncross_reference()
        aesurf2.write_card()
        aesurf2.cross_reference(model)
        aesurf2.raw_fields()
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()

    def test_flutter(self):
        """checks the FLUTTER/FLFACT cards"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        sid = 75
        method = 'PKNL'
        idensity = 76
        imach = 77
        ivelocity = 78

        # density, mach, velocity
        flutter1 = FLUTTER(sid, method, idensity, imach, ivelocity)
        flutter2 = FLUTTER.add_card(BDFCard(['FLUTTER', sid, method, idensity, imach,
                                             ivelocity]), comment='flutter card')
        assert flutter2.headers == ['density', 'mach', 'velocity'], flutter2.headers

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
        velocity = FLFACT(ivelocity, velocities)
        velocity.validate()
        velocity.write_card()
        assert velocity.min() == 3., velocities
        assert velocity.max() == 4., velocities
        model.flfacts[ivelocity] = velocity

        ikfreq = 79
        kfreqs = np.linspace(0.1, 0.2)
        card = ['FLFACT', ikfreq] + list(kfreqs)
        model.add_card(card, card[0])
        kfreq = model.FLFACT(ikfreq)
        kfreq.validate()
        kfreq.write_card()
        assert kfreq.min() == 0.1, kfreqs
        assert kfreq.max() == 0.2, kfreqs
        model.flfacts[ikfreq] = kfreq

        ikfreq2 = 80
        card = ['FLFACT', ikfreq2, 10., 'THRU', 20., 11]
        model.add_card(card, card[0])
        kfreq = model.FLFACT(ikfreq2)
        kfreq.validate()
        kfreq.write_card()
        assert kfreq.min() == 10., 'min=%s; card=%s factors=%s' % (kfreq.min(), card, kfreq.factors)
        assert kfreq.max() == 20., 'max=%s; card=%s factors=%s' % (kfreq.max(), card, kfreq.factors)
        model.flfacts[ikfreq] = kfreq

        ikfreq3 = 81
        factors = [10., 'THRU', 20., 10]
        kfreq = FLFACT(ikfreq3, factors)
        kfreq.validate()
        kfreq.write_card()
        assert kfreq.min() == 10., 'min=%s; factors=%s' % (kfreq.min(), factors)
        assert kfreq.max() == 20., 'max=%s; factors=%s' % (kfreq.max(), factors)
        model.flfacts[ikfreq] = kfreq

        # density, mach, rfreq
        card = ['FLUTTER', 85, 'KE', idensity, imach, ikfreq]
        model.add_card(card, card[0])
        model.cross_reference()
        flutter = model.Flutter(85)
        assert flutter.headers == ['density', 'mach', 'reduced_frequency'], flutter.headers
        flutter.write_card()

        model.uncross_reference()
        model.safe_cross_reference()

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
        msg = mkaero.write_card()
        lines = msg.strip().split('\n')
        expected = [
            '$mkaero',
            'MKAERO1       .5     .75',
            '              .1      .2      .3      .4      .5      .6      .7      .8',
            'MKAERO1       .5     .75',
            '              .9      1.     1.1',
        ]
        for line1, line2 in zip(lines, expected):
            assert line1 == line2, '\nline=%r\nexpected=%r'%  (line1, line2)

        expected2 = [
            '$mkaero1',
            'MKAERO1       .1      .2      .3      .4      .5      .6      .7      .8',
            '             .01     .02     .03',
            'MKAERO1       .9',
            '             .01     .02     .03',
        ]
        machs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        reduced_freqs = [0.01, 0.02, 0.03]
        model = BDF(debug=False)
        mkaero = model.add_mkaero1(machs, reduced_freqs, comment='mkaero1')
        mkaero.raw_fields()
        msg = mkaero.write_card()
        lines = msg.strip().split('\n')
        for line1, line2 in zip(lines, expected2):
            msg = '\nline    =%r\n' %  str(line1)
            msg += 'expected=%r\n%s' % (str(line2), msg)
            assert line1 == line2, msg

    def test_mkaero2(self):
        """checks the MKAERO2 card"""
        model = BDF()
        machs = [0.5, 0.75, 0.8]
        reduced_freqs = [0.1, 0.2, 0.3]
        mkaero = model.add_mkaero2(machs, reduced_freqs, comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()

        machs = [0.5, 0.75]
        reduced_freqs = [0.1, 0.2]
        mkaero = model.add_mkaero2(machs, reduced_freqs, comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()
        mkaero.raw_fields()

        mkaero = MKAERO2.add_card(BDFCard(['MKAERO2'] + machs + reduced_freqs), comment='mkaero2')
        mkaero.validate()
        mkaero.write_card()

    def test_diverg(self):
        """checks the DIVERG card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        sid = 100
        nroots = 21
        machs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        x0 = 3.
        V = 42.
        diverg = DIVERG(sid, nroots, machs, comment='divergence')
        diverg.validate()
        diverg.write_card()

        diverg = model.add_card(['DIVERG', sid, nroots] + machs, 'DIVERG', comment='divergence')
        #diverg.validate()
        #diverg.write_card()

    def test_trim_01(self):
        """checks the TRIM card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        sid = 100
        mach = 0.75
        q = 100.
        labels = ['ALPHA', 'ALPHA']
        uxs = [10., 20.]
        trim = TRIM(sid, mach, q, labels, uxs)
        with self.assertRaises(RuntimeError):
            trim.validate()

        labels = ['ALPHA']
        uxs = [10., 20.]
        trim = TRIM(sid, mach, q, labels, uxs)
        with self.assertRaises(RuntimeError):
            trim.validate()

        labels = ['ALPHA', 'BETA']
        uxs = [10., 20.]
        trim = TRIM(sid, mach, q, labels, uxs)
        trim.validate()
        trim.write_card()

        labels = ['ALPHA']
        uxs = [10.]
        trim = TRIM(sid, mach, q, labels, uxs, aeqr=3.0, comment='')
        trim.validate()
        trim.write_card()

        labels = ['ALPHA', 'BETA']
        uxs = [10., 20.]
        trim = TRIM(sid, mach, q, labels, uxs, aeqr=3.0, comment='')
        trim.validate()
        trim.write_card()

        model.add_card(['TRIM', sid, mach, q, labels[0], uxs[0]], 'TRIM', comment='$ trim')
        model.validate()
        model._verify_bdf(xref=False)

    def test_trim_02(self):
        """checks the TRIM card with a 2.5g pullup"""
        model = BDF()
        sid = 75
        mach = 0.75
        q = 100.
        labels = ['NZ']
        uxs = [2.5]
        trim = model.add_trim(sid, mach, q, labels, uxs, aeqr=0.0, comment='')
        trim.validate()

    def test_trim_03(self):
        """checks the TRIM card with a 2.5g pullup"""
        model = BDF(debug=None)
        sid = 75
        mach = 0.75
        q = 100.
        labels = ['URDD3', 'PITCH']
        uxs = [2.5, 0.0]
        trim1 = model.add_trim(sid, mach, q, labels, uxs, aeqr=0.0, comment='')

        sid = 76
        labels = ['URDD3', 'URDD5', 'PITCH']
        uxs = [2.5, 0.0, 0.0]
        # good
        trim2 = model.add_trim(sid, mach, q, labels, uxs, aeqr=0.0, comment='trim')

        model.add_aestat(1, 'URDD3', comment='aestat')
        model.add_aestat(2, 'URDD5', comment='aestat')
        model.add_aestat(3, 'PITCH', comment='aestat')
        model.add_aestat(4, 'ANGLEA', comment='aestat')

        #+--------+---------+-----------------------------+
        #| ANGLEA | ur (R2) | Angle of Attack             |
        #| YAW    | ur (R3) | Yaw Rate                    |
        #| SIDES  | ur (R3) | Angle of Sideslip           |
        #+--------+---------+-----------------------------+
        #| ROLL   | ůr (R1) | Roll Rate                   |
        #| PITCH  | ůr (R2) | Pitch Rate                  |
        #+--------+---------+-----------------------------+
        #| URDD1  | ür (T1) | Longitudinal (See Remark 3) |
        #| URDD2  | ür (T2) | Lateral                     |
        #| URDD3  | ür (T3) | Vertical                    |
        #| URDD4  | ür (R1) | Roll                        |
        #| URDD5  | ür (R2) | Pitch                       |
        #| URDD6  | ür (R3) | Yaw                         |
        #+--------+---------+-----------------------------+

        cid1 = 0
        label = 'DELTA'
        aesid = 5
        alid1 = 6
        model.add_aesurf(aesid, label, cid1, alid1)
        suport = model.add_suport([55, 66], ['3', '3'])
        str(suport)
        model.add_aelist(alid1, [100, 101, 102], comment='')
        model.add_grid(55, xyz=[0., 0., 0.])
        model.add_grid(66, xyz=[0., 0., 0.])
        model.validate()

        # why doesn't this work?
        with self.assertRaises(RuntimeError):
            trim1._verify(model.suport, model.suport1, model.aestats, model.aeparams,
                          model.aelinks, model.aesurf, xref=True)
        trim2._verify(model.suport, model.suport1, model.aestats, model.aeparams,
                      model.aelinks, model.aesurf, xref=True)
        model.write_bdf('trim.bdf')
        model2 = read_bdf('trim.bdf', debug=None)
        model2._verify_bdf(xref=True)
        model2.uncross_reference()
        model2._verify_bdf(xref=False)
        model2.cross_reference()
        model2._verify_bdf(xref=True)
        os.remove('trim.bdf')

        model2.uncross_reference()
        model2.safe_cross_reference()

    def test_gust(self):
        """checks the GUST card"""
        sid = 100
        dload = 200
        wg = 50.
        x0 = 3.
        V = 42.
        model = BDF()
        gust = model.add_gust(sid, dload, wg, x0, V=V, comment='gust load')
        gust.validate()
        gust.write_card()

        gust2 = GUST.add_card(BDFCard(['GUST', sid, dload, wg, x0, V]), comment='gust load')
        gust2.validate()
        gust2.write_card()


    def test_csschd(self):
        """checks the CSSCHD card"""
        model = BDF(debug=None)
        #sid = 10
        #aesid = 0
        #lalpha = None
        #lmach = None
        #lschd = None

        sid = 5
        aesid = 50
        lalpha = 12
        lmach = 15
        lschd = 25

        card = ['CSSCHD', sid, aesid, lalpha, lmach, lschd]
        bdf_card = BDFCard(card, has_none=True)
        csshcd_bad = CSSCHD(sid, aesid, lschd, lalpha='lalpha', lmach=4,
                            comment='')
        with self.assertRaises(TypeError):
            csshcd_bad.validate()
        csshcd_bad.lalpha = 4
        csshcd_bad.lmach = 5.0
        with self.assertRaises(TypeError):
            csshcd_bad.validate()
        csshcd_bad.lmach = 5
        csshcd_bad.validate()


        card = ['CSSCHD', sid, aesid, lalpha, lmach, lschd]
        bdf_card = BDFCard(card, has_none=True)
        csshcd1 = CSSCHD.add_card(bdf_card, comment='csschd card')
        csshcd1.validate()
        csshcd1.write_card()

        sid = 6
        csshcd2 = model.add_csschd(sid, aesid, lschd, lalpha=lalpha, lmach=lmach,
                                   comment='csschd card')

        label = 'ELEV'
        cid1 = 0
        alid1 = 37
        aesurf = model.add_aesurf(aesid, label, cid1, alid1, cid2=None, alid2=None,
                                  eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0,
                                  pllim=-np.pi/2., pulim=np.pi/2.,
                                  hmllim=None, hmulim=None,
                                  tqllim=None, tqulim=None, comment='aesurf')

        aelist = model.add_aelist(alid1, [1, 2, 3], comment='')

        aefact_sid = alid1
        Di = [0., 0.5, 1.]
        aefact_elev = model.add_aefact(aefact_sid, Di, comment='aefact')

        aefact_sid = lalpha
        Di = [0., 5., 10.]
        aefact_alpha = model.add_aefact(aefact_sid, Di, comment='aefact')

        aefact_sid = lmach
        Di = [0., 0.7, 0.8]
        aefact_mach = model.add_aefact(aefact_sid, Di, comment='aefact')

        aefact_sid = lschd
        Di = [0., 15., 30., 45.]
        aefact_delta = model.add_aefact(aefact_sid, Di, comment='aefact')

        model.cross_reference()
        csshcd2.write_card()
        #csshcd1.write_card()
        model.uncross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        model.safe_cross_reference()
        bdf_filename.seek(0)

        model2 = read_bdf(bdf_filename, punch=True, debug=False)

        bdf_filename2 = StringIO()
        model.write_bdf(bdf_filename2, size=16, close=False)

        #-----------
        csshcd3 = CSSCHD(sid, aesid, lschd, lalpha=None, lmach=None, comment='cssch card')
        csshcd3.write_card()
        with self.assertRaises(RuntimeError):
            csshcd3.validate()

    def test_monpnt(self):
        model = BDF(debug=False)
        name = 'test'
        label = 'test2'
        axes = '123'
        comp = 'WING'
        xyz = [0., 0., 0.]
        monpnt1 = model.add_monpnt1(name, label, axes, comp, xyz, cp=0,
                                    cd=None, comment='monpnt1')
        monpnt1.raw_fields()
        monpnt1.validate()

        Type = 'CQUAD4'
        table = 'STRESS'
        nddl_item = 42
        eid = 17
        monpnt2 = model.add_monpnt2(name, label, table, Type, nddl_item, eid,
                                    comment='monpnt2')
        monpnt2.raw_fields()
        monpnt2.validate()

        grid_set = 43
        elem_set = 44
        monpnt3 = model.add_monpnt3(name, label, axes, grid_set, elem_set,
                                   xyz, cp=0, cd=None,
                                   xflag=None, comment='monpnt3')
        monpnt3.raw_fields()
        monpnt3.validate()

        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)
        model.uncross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)
        model2 = read_bdf(bdf_filename, punch=True, debug=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
