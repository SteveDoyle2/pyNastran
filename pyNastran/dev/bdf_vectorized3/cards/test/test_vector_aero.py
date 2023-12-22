# coding: utf-8
# pylint: disable=R0914
"""tests aero cards"""
import os
from collections import defaultdict
import unittest
from io import StringIO
from typing import Optional, Any

import numpy as np
from cpylog import SimpleLogger

import pyNastran
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard, read_bdf # CORD2R, SET1,
from pyNastran.dev.bdf_vectorized3.test.test_bdf import run_bdf

from pyNastran.dev.bdf_vectorized3.cards.aero import (
    #AESURF,
    AEFACT, # AELIST, AEPARM,
    #CAERO1, CAERO2, CAERO3, CAERO4, #CAERO5,
    PAERO2, # PAERO1, PAERO3, PAERO4, PAERO5,
    ##AESURF, AESURFS,
    #AELINK,
    #AECOMP,
    ##SPLINE1, SPLINE2, #, SPLINE3, SPLINE4, SPLINE5,
    TRIM, # CSSCHD,
    #build_caero_paneling
)
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck


from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, MKAERO1, MKAERO2
from pyNastran.bdf.cards.aero.static_loads import AESTAT, AEROS # , CSSCHD, TRIM, TRIM2, DIVERG
#from pyNastran.bdf.cards.aero.utils import build_trim_load_cases
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh # build_structure_from_caero

IS_MATPLOTLIB = False
if IS_MATPLOTLIB:
    import matplotlib.pyplot as plt


ROOTPATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(ROOTPATH, '..', 'models')
TEST_PATH = os.path.join(ROOTPATH, 'bdf', 'cards', 'test')

COMMENT_BAD = 'this is a bad comment'
COMMENT_GOOD = 'this is a good comment\n'
RUN_PANELING = False

def build_caero_paneling(model: BDF) -> None:
    all_control_surface_name = None
    caero_control_surfaces = None
    out = None
    return all_control_surface_name, caero_control_surfaces, out

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
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        lines = ['AESTAT  502     PITCH']
        card = model._process_card(lines)
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
        model = BDF()
        model.aecomp.add_card(bdf_card, comment='aecomp card')
        aecomp1 = model.aecompl
        aecomp1.validate()
        aecomp1.write()

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

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AELIST', 76, 2000, 'THRU', 2010]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        #model.add_aesurf(aesurf)
        #model.add_aefact(aefact_elev)
        #model.add_aefact(aefact_alpha)
        #model.add_aefact(aefact_mach)
        #model.add_aefact(aefact_delta)

        #aecomp1.safe_cross_reference(model)
        #aecomp1.uncross_reference()

        #aecomp1.cross_reference(model)
        aecomp1.write()
        #aecomp1.uncross_reference()
        #aecomp1.write_card()

        model.validate()
        save_load_deck(model)

        #-----------
        if RUN_PANELING:
            aecomp2 = AECOMP(model)
            aecomp2.add(name, list_type, aelist_ids, comment='cssch card')
            aecomp2.validate()
            aecomp2.write()

            list_type = 'INVALID'
            aecomp3 = AECOMP(model)
            aecomp3.add(name, list_type, aelist_ids, comment='cssch card')
            #with self.assertRaises(RuntimeError):
                #aecomp3.validate()

            name = 'MYCOMP'
            list_type = 'AELIST'
            lists = 10
            model.add_aecomp(name, list_type, lists)

            lists = 42.0
            aecomp = AECOMP(model)
            aecomp.add(name, list_type, lists)
            with self.assertRaises(TypeError):
                aecomp.validate()


    def test_aefact_1(self):
        """checks the AEFACT card"""
        data = ['AEFACT', 97, .3, 0.7, 1.0]
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        aefact = model.aefact
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)

        data = ['AEFACT', '98', '.3', '0.7', '1.0']
        model.add_card(data, data[0], COMMENT_GOOD, is_list=True)
        model.setup()

        #msg = '$this is a bad comment\n'
        msg = 'AEFACT        97      .3      .7      1.\n'
        aefact97 = aefact.slice_card_by_aefact_id(97)
        aefact98 = aefact.slice_card_by_aefact_id(98)
        self.assertTrue(all(aefact97.fractions == [.3, .7, 1.0]))
        self.assertTrue(all(aefact98.fractions == [.3, .7, 1.0]))

        out = aefact97.write(8, is_double=False)
        self.assertEqual(msg, out)

        #msg = '$this is a good comment\n'
        msg = 'AEFACT        98      .3      .7      1.\n'
        out = aefact98.write(8, None)
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

        fractions = [1., 2., 3.]
        aefact.clear()
        aefact.add(200, fractions, comment='')
        aefact.parse_cards()
        aefact.validate()
        aefact.write()
        #model = BDF()
        #aefact.cross_reference(model)
        #aefact.write_card()
        #aefact.uncross_reference()
        #aefact.write_card()


    def test_aelink_1(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        aelink = model.aelink
        idi = 10
        label = 'CS'
        independent_labels = ['A', 'B', 'C']
        linking_coefficients = [1.0, 2.0]

        aelink = model.aelink
        with self.assertRaises(AssertionError):
            aelink = aelink.add(idi, label, independent_labels,
                                linking_coefficients, comment='')
        #with self.assertRaises(AssertionError):
            #aelink.parse_cards()

        #assert aelink.aelink_id == idi
        #str(aelink)
        #aelink.write()

        aelink.clear()
        card = ['AELINK', idi, label, independent_labels[0], linking_coefficients[0],
                independent_labels[1], linking_coefficients[1], independent_labels[2]]
        with self.assertRaises(AssertionError):
            model.add_card(card, 'AELINK')

        aelink.clear()
        card = ['AELINK', idi, label, independent_labels[0], linking_coefficients[0],
                independent_labels[1], linking_coefficients[1]]
        model.add_card(card, 'AELINK', comment='cat')
        #print(model.aelinks[idi])
        #assert model.aelinks[idi][0].comment == '$cat\n', 'comment=%r' % str(model.aelinks[idi][0].comment)

        #-------------------------------
        aelink.clear()
        idi = 11
        label = 'LABEL'
        independent_labels = ['pig', 'frog', 'dog']
        linking_coefficients = []
        with self.assertRaises(AssertionError):
            aelink.add(idi, label, independent_labels, linking_coefficients)
        #with self.assertRaises(RuntimeError):
            #model.validate()
        aelink.clear()

        linking_coefficients = [1.0, 2.0, 3.0]
        aelink2i = model.add_aelink(idi, label, independent_labels, linking_coefficients)
        model.setup()
        aelink2 = aelink.slice_card_by_id(idi)
        assert np.allclose(aelink2.linking_coefficients, [1., 2., 3.])

        #-------------------------------
        idi = 'ALWAYS'
        label = 'LABEL'
        independent_labels = ['pig', 'frog', 'dog']
        linking_coefficients = [1.0, 2.0, 3.0]
        model.add_aelink(idi, label, independent_labels, linking_coefficients)

        sid = 10
        mach = 0.5
        q = 10.
        labels = ['URDD3']
        uxs = [1.0]
        model.add_trim(sid, mach, q, labels, uxs, aeqr=1.0, trim_type=1, comment='')

        sid = 11
        model.add_trim(sid, mach, q, labels, uxs, aeqr=1.0, trim_type=1, comment='')

        cid1 = 0
        aesid_label_alids = [
            (1, 'CS', 101),
            (2, 'LABEL', 101),
            (3, 'A', 101),
            (4, 'B', 101),
            #(5, 'C', 101),
            (6, 'PIG', 101),
            (7, 'FROG', 101),
            (8, 'DOG', 101),
        ]
        for (aesid, label, aelist_id1) in aesid_label_alids:
            model.add_aesurf(aesid, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
                             eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0,
                             pllim=-np.pi/2., pulim=np.pi/2.,
                             hmllim=None, hmulim=None,
                             tqllim=None, tqulim=None, comment='')
        elements = [100]
        model.add_aelist(101, elements)
        model.validate()
        model.cross_reference()

    def test_aelink_2(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        idi = 31
        label = 'LABEL'
        independent_labels = ['pig', 'frog', 'dog']
        linking_coefficients = [1.0, 2.0, 3.0]
        model.add_aelink(idi, label, independent_labels, linking_coefficients)

        sid = 31
        mach = 0.5
        q = 10.
        labels = ['URDD3']
        uxs = [1.0]
        model.add_trim(sid, mach, q, labels, uxs, aeqr=1.0, trim_type=1, comment='')

        cid1 = 0
        aesid_label_alids = [
            (1, 'LABEL', 101),
            (2, 'PIG', 101),
            (3, 'FROG', 101),
            (4, 'DOG', 101),
        ]
        for (aesid, label, aelist_id1) in aesid_label_alids:
            model.add_aesurf(aesid, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
                             eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0,
                             pllim=-np.pi/2., pulim=np.pi/2.,
                             hmllim=None, hmulim=None,
                             tqllim=None, tqulim=None, comment='')
        elements = [100]
        model.add_aelist(101, elements)
        model.setup()

        save_load_deck(model, run_renumber=False)

    def test_aelist_1(self):
        """checks the AELIST card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        aelist = model.aelist
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)
        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1201, 1202]
        aelisti = model.add_aelist(74, elements)
        model.setup(run_geom_check=False)
        #model.aelist.validate()
        model.aelist.write()
        aelist75 = model.aelist.slice_card_by_aelist_id(75)
        #print(aelist.elements)
        #print(elements)
        assert np.array_equal(elements, aelist75.elements)

        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1108, 1202]
        data = ['AELIST', 76, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1108, 1202]
        model.add_card(data, data[0], COMMENT_BAD, is_list=True)
        model.setup()
        aelist76 = model.aelist.slice_card_by_aelist_id(76)
        #print(aelist76 .elements)
        #print(elements)
        assert not np.array_equal(elements, aelist76.elements)

        elements = list(set(elements))
        elements.sort()
        assert np.array_equal(elements, aelist76.elements)

        elements = [1000, 1000, 1000, 2000, 1000, 2000]
        aelist.clear()
        #aelist = AELIST(model)
        aelist.add(75, elements)
        aelist.parse_cards()
        aelist.clean_ids()
        str(aelist.write())

        elements = 42
        aelist.clear()
        aelist.add(76, elements)

        elements = 42.0
        #with self.assertRaises(TypeError):
            #AELIST(77, elements)

    def test_aeparm_1(self):
        """checks the AEPARM card"""
        aeparm_id = 100
        model = BDF(debug=False)
        aeparm = model.aeparm
        aeparm.add_card(BDFCard(['AEPARM', aeparm_id, 'THRUST', 'lb']),
                        comment='aeparm_comment')

        aeparmi = model.add_aeparm(aeparm_id, 'THRUST', 'lb', comment='aeparm_comment')
        model.setup()
        assert np.allclose(aeparm.aeparm_id, aeparm_id)
        aeparm.validate()
        #aeparm.cross_reference(None)
        #aeparm.uncross_reference()
        #aeparm.safe_cross_reference(None)
        aeparm.write()
        save_load_deck(model)

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

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.add_card(bdf_card, 'AESURFS', comment='aesurfs',
                       is_list=True, has_none=True)
        aesurfs = model.add_aesurfs(aesid, label, list1, list2, comment='aesurfs')
        str(aesurfs)
        aesurfs.write()

        model.add_set1(6002, [1, 2, 3])
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.setup()

        model.validate()
        save_load_deck(model)

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

        assert aero.is_symmetric_xy is False
        assert aero.is_symmetric_xz is False
        assert aero.is_anti_symmetric_xy is False
        assert aero.is_anti_symmetric_xz is False

        #aero.set_ground_effect(True)
        #assert aero.is_symmetric_xy is False
        #assert aero.is_symmetric_xz is False
        #assert aero.is_anti_symmetric_xy is True
        #assert aero.is_anti_symmetric_xz is False

        #aero.set_ground_effect(False)
        #assert aero.is_symmetric_xy is False
        #assert aero.is_symmetric_xz is False
        #assert aero.is_anti_symmetric_xy is False
        #assert aero.is_anti_symmetric_xz is False

        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=1, sym_xy=1,
                    comment='aero card')
        assert aero.is_symmetric_xy is True
        assert aero.is_symmetric_xz is True
        assert aero.is_anti_symmetric_xy is False
        assert aero.is_anti_symmetric_xz is False

        aero = AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=-1, sym_xy=-1,
                    comment='aero card')
        assert aero.is_symmetric_xy is False
        assert aero.is_symmetric_xz is False
        assert aero.is_anti_symmetric_xy is True
        assert aero.is_anti_symmetric_xz is True

        aero.set_ground_effect(True)

    def test_aero_2(self):
        """checks the AERO card"""
        acsid = 0
        velocity = None
        cref = 1.0
        rho_ref = 1.0
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

        #model = BDF()
        #aero.cross_reference(model)
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

        cref = 1
        bref = 2
        sref = 3
        acsid = 42.
        rcsid = 43.
        sym_xz = 44.
        sym_xy = 45.
        aeros = AEROS(cref, bref, sref, acsid, rcsid, sym_xz=sym_xz, sym_xy=sym_xy)
        with self.assertRaises(TypeError):
            aeros.validate()

    def test_caero1_paneling_nspan_nchord_1(self):
        """checks the CAERO1/PAERO1/AEFACT card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero1 = model.caero1
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 10000000
        caeroi = model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                                  cp=0, nspan=3, lspan=0, nchord=2, lchord=0, comment='')
        model.setup()
        npoints, nelements = caero1.get_npanel_points_elements()
        npoints_expected = 12 # 4*3
        nelements_expected = 6 # 2*3

        x, y = caero1.xy
        chord_expected = np.array([0., 0.5, 1.])
        span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        assert np.allclose(x, chord_expected)
        assert np.allclose(y, span_expected)
        assert npoints_expected == npoints
        assert nelements_expected == nelements

    def test_caero1_paneling_nspan_lchord(self):
        """checks the CAERO1/PAERO1/AEFACT card"""
        fig, ax = _setup_aero_plot()

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero1 = model.caero1
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 10000000
        chord_aefact_id = 10000
        model.add_aefact(chord_aefact_id, [0., 0.5, 1.0])
        caeroi = model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                                  cp=0,
                                  nspan=3, lspan=0,
                                  nchord=0, lchord=chord_aefact_id, comment='')
        model.cross_reference()

        npoints, nelements = caero1.get_npanel_points_elements()
        npoints_expected = 12 # 4*3
        nelements_expected = 6 # 2*3
        assert npoints_expected == npoints
        assert nelements_expected == nelements
        #del model.caeros[eid]
        #del model.aefacts[chord_aefact_id]
        #points, elements = caero1.panel_points_elements()
        x, y = caero1.xy
        chord_expected = np.array([0., 0.5, 1.])
        span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        assert np.allclose(x, chord_expected)
        assert np.allclose(y, span_expected)
        if IS_MATPLOTLIB:
            caero1.plot(ax)
            fig.show()

    def test_caero1_paneling_transpose(self):
        fig, ax = _setup_aero_plot()

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero1 = model.caero1
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        #['CAERO1', '2000', '2000', '0', '15', '10', '1', '0', None, '7.314386', '0.', '-0.18288', '1.463854', '8.222755', '1.573341', '-0.18288', '0.365963']
        #card_lines = [
            #'CAERO1,2000,2000,0,15,10,1,0,1',
            #'+,7.314386,0.,-0.18288,1.463854,8.222755,1.573341,-0.18288,0.365963',
        #]
        #model.add_card(card_lines, 'CAERO1', comment='', ifile=None, is_list=False, has_none=True)

        eid = 2000
        #caero = model.caeros[eid]
        #print(caero.get_stats())

        pid = 1
        igroup = 1
        p1 = [7.3, 0., 0.]
        p4 = [8.2, 1.6, 0.]
        x12 = 1.4
        x43 = 0.3
        model.add_paero1(pid, caero_body_ids=None, comment='')

        caeroi = model.add_caero1(
            eid, pid, igroup, p1, x12, p4, x43,
            cp=0, nspan=5, lspan=0, nchord=2, lchord=0, comment='')

        caeroi = model.add_caero1(
            eid+100, pid, igroup, p1, x12, p4, x43,
            cp=0, nspan=5, lspan=0, nchord=2, lchord=0, comment='')
        caero1.parse_cards()
        caero1.validate()
        model.setup()
        x, y = caero1.xy
        x_expected = np.array([0., 0.5, 1.])
        y_expected = np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])
        assert np.allclose(x, x_expected)
        assert np.allclose(y, y_expected)

        #print(caero.get_stats())
        #caero.cross_reference(model)
        if RUN_PANELING:
            all_control_surface_name, caero_control_surfaces, out = build_caero_paneling(model)
            box_id_to_caero_element_map_expected = {
                2000: np.array([0, 3, 4, 1]),
                2001: np.array([1, 4, 5, 2]),
                2002: np.array([3, 6, 7, 4]),
                2003: np.array([4, 7, 8, 5]),
                2004: np.array([ 6,  9, 10,  7]),
                2005: np.array([ 7, 10, 11,  8]),
                2006: np.array([ 9, 12, 13, 10]),
                2007: np.array([10, 13, 14, 11]),
                2008: np.array([12, 15, 16, 13]),
                2009: np.array([13, 16, 17, 14]),
            }
            for key, data in out.box_id_to_caero_element_map.items():
                assert np.array_equal(data, box_id_to_caero_element_map_expected[key])

        #all_control_surface_name, caero_control_surfaces, out = build_caero_paneling(model)
        if IS_MATPLOTLIB:
            caero1.plot(ax)
            fig.show()
        #x = 1

    def test_caero1_paneling_multi(self):
        """checks the CAERO1/PAERO1/AEFACT card"""
        fig, ax = _setup_aero_plot()

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 1000
        chord_aefact_id = 10000
        model.add_aefact(chord_aefact_id, [0., 0.5, 1.0])
        model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                         cp=0,
                         nspan=3, lspan=0,
                         nchord=0, lchord=chord_aefact_id, comment='')

        eid = 2000
        p1 = [1., 16., 0.]
        p4 = [1., 30., 0.]
        x12 = 1.
        x43 = 1.
        caero1b = model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                                   cp=0,
                                   nspan=3, lspan=0,
                                   nchord=2, lchord=0, comment='')
        model.cross_reference()

        caero1a = model.caero1.slice_card_by_index(0)
        npoints, nelements = caero1a.get_npanel_points_elements()
        npoints_expected = 12 # 4*3
        nelements_expected = 6 # 2*3
        assert npoints_expected == npoints
        assert nelements_expected == nelements
        #del model.caeros[eid]
        #del model.aefacts[chord_aefact_id]
        points, elements = caero1a.panel_points_elements()
        x, y = caero1a.xy
        chord_expected = np.array([0., 0.5, 1.])
        span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        assert np.allclose(x, chord_expected)
        assert np.allclose(y, span_expected)
        if IS_MATPLOTLIB:
            caero1a.plot(ax)
            caero1b.plot(ax)
            fig.show()
            x = 1

    def test_caero1_paneling_nspan_nchord_2(self):
        """checks the CAERO1/PAERO1/AEFACT card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero1 = model.caero1
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        # basic
        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        fig, ax = _setup_aero_plot(fig_id=None)
        unused_paero = model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 1000
        aelist_id = 10
        aesurf_id = 10
        caeroi = model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                                  cp=0, nspan=3, lspan=0, nchord=1, lchord=0, comment='')
        model.setup()
        x, y = caero1.xy
        chord_expected = np.array([0., 1.])
        span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        assert np.allclose(x, chord_expected)
        assert np.allclose(y, span_expected)

        elements = [1001, 1003, 1005]
        unused_aelist = model.add_aelist(aelist_id, elements)

        label = 'FLAP'
        cid1 = 0
        aelist_id1 = aelist_id
        unused_aesurf = model.add_aesurf(
            aesurf_id, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
            eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0, pllim=-np.pi/2., pulim=np.pi/2.,
            hmllim=None, hmulim=None, tqllim=None, tqulim=None, comment='')
        model.setup()
        model.cross_reference()
        points, elements = caero1.panel_points_elements()
        if IS_MATPLOTLIB:
            caero1.plot(ax)
            fig.show()

    def test_caero3_paneling(self):
        """checks the CAERO3/PAERO1/AEFACT card"""
        fig, ax = _setup_aero_plot()
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero3 = model.caero3
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        #model.add_paero1(pid, caero_body_ids=None, comment='')
        ncontrol_surfaces = 0
        nbox = 7
        x = []
        y = []
        model.add_paero3(pid, nbox, ncontrol_surfaces, x, y, comment='')
        eid = 1000
        list_w = None
        caeroi = model.add_caero3(eid, pid, p1, x12, p4, x43,
                                 cp=0, list_w=list_w, list_c1=None, list_c2=None, comment='')
        caero3.validate()
        model.setup()
        #caero.cross_reference(model)
        npoints, nelements = caero3.get_npanel_points_elements()
        npoints_expected = 24 # (2+1)*(7+1)
        nelements_expected = 14 # 2*7

         # hardcoded
        nchord_elements = 2
        nchord_points = nchord_elements + 1

        x2, y2 = caero3.xy
        #chord_expected = np.array([0., 0.5, 1.])
        #span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        #assert np.allclose(x, chord_expected)
        #assert np.allclose(y, span_expected)
        assert npoints_expected == npoints
        assert nelements_expected == nelements
        points, elements = caero3.panel_points_elements()
        if IS_MATPLOTLIB:
            caero.plot(ax)
            fig.show()
        x = 1

    def test_caero4_paneling(self):
        """checks the CAERO4/PAERO4 card"""
        fig, ax = _setup_aero_plot()
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero4 = model.caero4
        paero4 = model.paero4
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        #model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 100
        caero = model.add_caero4(eid, pid, p1, x12, p4, x43,
                                 cp=0, nspan=2, lspan=0, comment='')
        docs = []
        caocs = []
        gapocs = []
        model.add_paero4(pid, docs, caocs, gapocs,
                         cla=0, lcla=0, circ=0, lcirc=0, comment='')
        model.setup()
        paero4.validate()
        caero4.validate()
        #caero.cross_reference(model)
        npoints, nelements = caero4.get_npanel_points_elements()
        npoints_expected = 6 # 4*3
        nelements_expected = 2 # 2*3

        x, y = caero4.xy
        #chord_expected = np.array([0., 0.5, 1.])
        #span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        #assert np.allclose(x, chord_expected)
        #assert np.allclose(y, span_expected)
        assert npoints_expected == npoints
        assert nelements_expected == nelements
        if IS_MATPLOTLIB:
            caero4.plot(ax)
            fig.show()
        x = 1

    def test_caero5_paneling(self):
        """checks the CAERO4/PAERO4 card"""
        fig, ax = _setup_aero_plot()
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero5 = model.caero5
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        #model.add_paero1(pid, caero_body_ids=None, comment='')

        eid = 100
        nspan = 10
        model.add_caero5(eid, pid, p1, x12, p4, x43, cp=0,
                         nspan=nspan, lspan=0, ntheory=0, nthick=0, comment='')
        caoci = []
        model.add_paero5(pid, caoci, nalpha=0, lalpha=0,
                         nxis=0, lxis=0, ntaus=0, ltaus=0, comment='')
        model.paero5.validate()
        model.caero5.validate()
        #model.caero5.cross_reference(model)

        model.setup()
        npoints, nelements = model.caero5.get_npanel_points_elements()
        npoints_expected = (nspan + 1) * 2
        nelements_expected = nspan # 2*1
        assert npoints_expected == npoints
        assert nelements_expected == nelements

        #x, y = caero.xy
        #chord_expected = np.array([0., 0.5, 1.])
        #span_expected = np.array([0., 1 / 3, 2 / 3, 1.])
        #assert np.allclose(x, chord_expected)
        #assert np.allclose(y, span_expected)

        if RUN_PANELING:
            all_control_surface_name, caero_control_surfaces, out = build_caero_paneling(model)
            box_id_to_caero_element_map_expected = {
                100: np.array([0, 2, 3, 1]),
                101: np.array([2, 4, 5, 3]),
                102: np.array([4, 6, 7, 5]),
                103: np.array([6, 8, 9, 7]),
                104: np.array([ 8, 10, 11,  9]),
                105: np.array([10, 12, 13, 11]),
                106: np.array([12, 14, 15, 13]),
                107: np.array([14, 16, 17, 15]),
                108: np.array([16, 18, 19, 17]),
                109: np.array([18, 20, 21, 19]),
            }
            assert len(box_id_to_caero_element_map_expected) == len(out.box_id_to_caero_element_map)
            for key, data in out.box_id_to_caero_element_map.items():
                expected_data = box_id_to_caero_element_map_expected[key]
                assert np.array_equal(data, expected_data)

        points, elements = caero5.panel_points_elements()
        if IS_MATPLOTLIB:
            caero5.plot(ax)
            fig.show()
        x = 1

    def _test_caero7_paneling(self):
        """checks the CAERO7/PAERO7T card"""
        fig, ax = _setup_aero_plot()

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        cref = 1.0
        bref = 1.0
        sref = 1.0
        model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [1., 15., 0.]
        x12 = 1.
        x43 = 1.
        #model.add_paero7

        eid = 100
        chord_aefact_id = 10000
        model.add_aefact(chord_aefact_id, [0., 0.5, 1.0])
        label = 'panel'
        nspan = 2
        nchord = 4
        caero = model.add_caero7(eid, label, p1, x12, p4, x43, cp=0,
                                 nspan=nspan, nchord=nchord, lspan=0,
                                 p_airfoil=None, ztaic=None, comment='')
        model.cross_reference()
        model.validate()

        xy = caero.xy
        npoints, nelements = caero.get_npanel_points_elements()
        npoints_expected = (nspan + 1) * (nchord + 1)
        nelements_expected = nspan * nchord
        #npoints_expected = 15 # 4*3
        #nelements_expected = 8 # 2*3
        assert npoints_expected == npoints
        assert nelements_expected == nelements


        points, elements = caero.panel_points_elements()
        x, y = caero.xy
        chord_expected = np.array([0., 0.25, 0.5, 0.75, 1.])
        span_expected = np.array([0., 0.5, 1.])
        assert np.allclose(x, chord_expected)
        assert np.allclose(y, span_expected)
        if IS_MATPLOTLIB:
            caero.plot(ax)
            fig.show()

        all_control_surface_name, caero_control_surfaces, out = build_caero_paneling(model)
        box_id_to_caero_element_map_expected = {
            100: np.array([0, 5, 6, 1]),
            101: np.array([1, 6, 7, 2]),
            102: np.array([2, 7, 8, 3]),
            103: np.array([3, 8, 9, 4]),

            104: np.array([5, 10, 11, 6]),
            105: np.array([6, 11, 12, 7]),
            106: np.array([7, 12, 13, 8]),
            107: np.array([8, 13, 14, 9]),
        }
        assert len(box_id_to_caero_element_map_expected) == len(out.box_id_to_caero_element_map)
        for box_id, data in out.box_id_to_caero_element_map.items():
            expected_data = box_id_to_caero_element_map_expected[box_id]
            assert np.array_equal(data, expected_data)
        x = 1

    def test_caero1_1(self):
        """checks the CAERO1/PAERO1/AEROS/AEFACT card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero1 = model.caero1
        paero1 = model.paero1
        aefact = model.aefact
        model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                nxref_errors=0, stop_on_xref_error=True)

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

        caero1.add_card(BDFCard(['CAERO1', eid, pid, cp, nspan, nchord, lspan, lchord,
                                 igid, ] + p1 + [x12] + p4 + [x43]))
        caero1.parse_cards()
        caero1.validate()

        eid = 2
        caero1.add_card(BDFCard(['CAERO1', eid, pid, None, nspan, nchord, lspan, lchord,
                                 igid, ] + p1 + [x12] + p4 + [x43]))
        caero1.parse_cards()
        caero1.validate()

        eid = 1
        model.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        caero1.parse_cards()
        #caero1c.raw_fields()
        caero1.validate()
        caero1.write()
        #model.caeros[eid] = caero1c

        eid = 4
        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 1., 0.]
        p4 = [0.1, 1., 0.]
        nspan = 5
        nchord = 10
        igid = -1
        caero1 = model.caero1
        caero1.add_quad(eid, pid, nspan, nchord, igid, p1, p2, p3, p4,
                         cp=cp, spanwise='y', comment='')
        caero1.parse_cards()
        caero1.validate()

        eid = 5
        span = 0.1
        chord = 0.05
        igid = -1
        caero1.add_quad(eid, pid, span, chord, igid, p1, p2, p3, p4,
                         cp=cp, spanwise='y', comment='')
        caero1.parse_cards()
        caero1.validate()

        eid = 6
        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 0., 1.]
        p4 = [0.1, 0., 1.]
        span = 0.1
        chord = 0.05
        igid = -1
        caero1.add_quad(eid, pid, span, chord, igid, p1, p2, p3, p4,
                                  cp=cp, spanwise='z', comment='')
        caero1.parse_cards()
        caero1.validate()
        caero1.flip_normal_by_element_id()
        model.setup()

        model.add_cord2r(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                         comment='')
        model.setup()
        model.coord.validate()

        eid = 7
        p1 = [0., 0., 0.]
        p2 = [1., 0., 0.]
        p3 = [0.2, 0., 1.]
        p4 = [0.1, 0., 1.]
        span = 0.1
        chord = 0.05
        igid = -1
        cp = None
        caero1.add_quad(eid, pid, span, chord, igid,
                        p1, p2, p3, p4,
                        cp=cp, spanwise='z', comment='')
        caero1.get_points()

        # caero1c is set as eid=1
        model.validate()
        # ------------------------------------------------
        eid =  1000
        igroup = 1
        lspan_lchord = 1
        fractions = np.linspace(0., 1., num=11)
        model.add_aefact(lspan_lchord, fractions, comment='')
        caero1i = model.add_caero1(eid, pid, igroup, p1, x12, p4, x43, cp=0,
                                   nspan=0, lspan=lspan_lchord,
                                   nchord=0, lchord=lspan_lchord, comment='')

        paero1i = model.add_paero1(pid, caero_body_ids=None, comment='')
        model.setup()
        paero1.validate()
        paero1.write()
        #model.paeros[pid] = paero

        #acsid = 0.
        #velocity = None
        cref = 1.0
        bref = 2.0
        sref = 100.
        aeros = model.add_aeros(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0,
                                sym_xy=0, comment='aeros')
        aeros.validate()
        aeros.write_card()
        model.aeros = aeros

        aefacti = model.add_aefact(lspan, [0., 1., 2., 3., 4., 5.])
        #aefact.validate()
        #model.aefacts[lspan] = aefact
        aefacti = model.add_aefact(lchord, [2., 3., 4., 5., 6., 7.])
        aefact.validate()
        model.setup()
        #model.aefacts[lchord] = aefact

        #paero.cross_reference(model)
        #caero1c.cross_reference(model)
        caero1.get_npanel_points_elements()
        caero1.get_points()
        caero1.panel_points_elements()

        caero1.write()
        #model.uncross_reference()
        model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()
        xref_errors = defaultdict(list)
        #caero1c.safe_cross_reference(model, xref_errors)
        caero1.panel_points_elements()
        #caero1c.raw_fields()
        min_max_eid = caero1.min_max_eid
        caero1.remove_duplicates(inplace=True)

        caero1a = caero1.slice_card_by_index(0)
        min_max_eid = caero1a.min_max_eid
        self.assertEqual(min_max_eid, (1, 26))
        #print('min_eid, max_eid', min_eid, max_eid)

        RUN_CAERO1 = False
        if RUN_CAERO1:
            points = [
                [0., 0., 0.], # p1
                [10., 0., 0.],
                [10., 20., 0.],
                [5., 20., 0.],
            ]
            caero1.set_points(points)
            caero1.get_points()
        str(caero1.write())

        nspan = None
        lspan = None
        modelh = BDF(log=log)
        caero1h = modelh.caero1
        caero1h_id = modelh.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=None,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modelh.setup()

        with self.assertRaises(ValueError):
            caero1h.validate()

        nspan = 5
        lspan = 5
        modeli = BDF(log=log)
        caero1i = modeli.caero1
        caero1i_id = modeli.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modeli.setup()
        with self.assertRaises(ValueError):
            caero1i.validate()

        nspan = 5
        nchord = None
        lchord = None
        modelj = BDF(log=log)
        caero1j = modelj.caero1
        caero1j_id = modelj.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modelj.setup()
        with self.assertRaises(ValueError):
            caero1j.validate()

        nchord = 10
        lchord = 10
        modelk = BDF(log=log)
        caero1k = modelk.caero1
        caero1k_id = modelk.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modelk.setup()
        with self.assertRaises(ValueError):
            caero1k.validate()

        lspan = None
        lchord = None
        nspan = 10
        nchord = 10
        p1 = [0., 0., 0., 0.]
        modell = BDF(log=log)
        caero1l = modell.caero1
        caero1l_id = modell.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')

        with self.assertRaises(ValueError):
            modell.setup()
        #with self.assertRaises(ValueError):
            #caero1l.validate()

        p1 = [0., 0., 0.]
        p4 = [1., 2., 3., 4.]
        modelm = BDF(log=log)
        caero1m = modelm.caero1
        caero1m_id = modelm.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        with self.assertRaises(ValueError):
            modelm.setup()
        #with self.assertRaises(AssertionError):
            #caero1m.validate()

        p4 = [1., 2., 3.]
        eid = 8
        nspan = 1
        nchord = 1
        lchord = None
        lspan = None
        modeln = BDF(log=log)
        caero1_1by1 = modeln.caero1
        caero1_1by1_id = modeln.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modeln.setup()
        caero1_1by1.validate()
        assert caero1_1by1.shape == (1, 1)
        caero1_1by1.get_points()

        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 1.
        x43 = 1.
        eid = 1
        nspan = 3
        nchord = 2
        lchord = None
        lspan = None
        modelo = BDF(log=log)
        caero1_2x3 = modelo.caero1
        caero1_2x3_id = modelo.add_caero1(
            eid, pid, igid, p1, x12, p4, x43, cp=cp,
            nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
            comment='caero1')
        modelo.setup()
        caero1_2x3.validate()
        assert caero1_2x3.shape == (2, 3), caero1_2x3.shape
        if RUN_CAERO1:
            caero1_2x3._init_ids()
        points = caero1_2x3.get_points()
        assert len(points) == 4
        save_load_deck(model)

    def test_spline1(self):
        """checks the SPLINE1 card"""
        eid = 1
        caero_id = 100
        box1 = 1
        box2 = 10
        setg = 42
        #model = BDF(debug=False)
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        spline1 = model.spline1
        spline2 = model.spline2
        caero1 = model.caero1
        caero2 = model.caero2
        paero1 = model.paero1
        paero2 = model.paero2
        aefact = model.aefact


        #pid = 1
        #igroup = 1
        #p1 = [0., 0., 0.]
        #p4 = [0., 10., 0.]
        #x12 = 1.
        #x43 = 0.5
        #model.add_caero1(eid, pid, igroup, p1, x12, p4, x43,
                         #cp=0, nspan=1, lspan=0, nchord=1, lchord=0, comment='')
        #model.add_paero1(pid, caero_body_ids=None, comment='')

        model.add_spline1(eid, caero_id, box1, box2, setg, dz=0., method='IPS',
                          usage='BOTH', nelements=10,
                          melements=10, comment='$ spline1')
        model.setup(run_geom_check=False)
        spline1.validate()
        spline1.write(size=8, is_double=False)
        #spline.raw_fields()

        #model.splines[eid] = spline

        pid = 10
        igid = 1
        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 4.
        x43 = 3.
        cid = 1
        caero1i = model.add_caero1(caero_id, pid, igid, p1, x12, p4, x43,
                                   cp=cid, nspan=5,
                                   lspan=0, nchord=6, lchord=0,
                                   comment='')
        caero_body_ids = [3]
        unused_paero = model.add_paero1(pid, caero_body_ids=caero_body_ids, comment='')
        origin = None
        zaxis = None
        xzplane = None
        model.add_cord2r(cid, origin, zaxis, xzplane, rid=0, comment='')
        velocity = 0.0
        cref = 1.0
        rho_ref = 1.225
        model.add_aero(velocity, cref, rho_ref, acsid=0, sym_xz=0, sym_xy=0,
                       comment='')

        setg = 42
        ids = [100, 101, 102]
        model.add_set1(setg, ids, is_skin=False, comment='')
        model.add_grid(100, [0., 0., 0.])
        model.add_grid(101, [0., 0., 0.])
        model.add_grid(102, [0., 0., 0.])

        #------------------
        # CAERO2
        eid = 3
        caero = 3
        id1 = 21
        id2 = 35
        setg = 43
        spline2i = model.add_spline2(eid, caero, id1, id2, setg, dz=0.0, dtor=1.0, cid=1,
                                    #dthx=None, dthy=None,
                                    usage='BOTH', comment='')
        spline2.validate()

        pid = 3
        caero2i = model.add_caero2(caero, pid, igid, p1, x12, cp=1, nsb=4,
                                  nint=4, lsb=0, lint=0, comment='')
        caero2.parse_cards()
        caero2.validate()

        orient = 'ZY'
        width = 1.0
        AR = 2.0
        thi = []
        thn = []
        paero2i = model.add_paero2(pid, orient, width, AR, thi, thn,
                                   lrsb=10, lrib=None, lth=None, comment='')
        paero2.parse_cards()
        paero2.validate()

        sid = 10
        Di = [0., 1.0, 2.0, 3.0, 0.]
        aefacti = model.add_aefact(sid, Di, comment='')
        aefact.parse_cards()
        aefact.validate()

        model.add_set1(setg, ids, is_skin=False, comment='')

        model.cross_reference(model)
        caero1.panel_points_elements()
        caero2.get_points_elements_3d()
        save_load_deck(model)


    def test_spline2(self):
        """checks the SPLINE2 card"""
        #+---------+------+-------+-------+-------+------+----+------+-----+
        #| SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
        #|         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
        #+---------+------+-------+-------+-------+------+----+------+-----+
        #| SPLINE2 |   5  |   8   |  12   | 24    | 60   | 0. | 1.0  |  3  |
        #|         |  1.  |       |       |       |      |    |      |     |
        #+---------+------+-------+-------+-------+------+----+------+-----+

        cid = 3
        origin = [0., 0., 0.]
        xaxis = [1., 0., 0.]
        xyplane = [0., 1., 0.]
        model = BDF(debug=True, log=None, mode='msc')
        spline2 = model.spline2
        caero2 = model.caero2
        coord = model.coord.add_axes(
            cid, 'R', rid=0, origin=origin,  # cord2r
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
        igid = 1 # None
        caero2_id = model.add_caero2(eid, pid, igid, p1, x12,
                                     cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                                     comment='this is a caero')
        model.setup()
        caero2.validate()
        #caero = CAERO2(eid, pid, cp, nsb, nint, lsb, lint, igid, p1, x12)

        sid = 60
        ids = [7, 13]
        set_obj = model.add_set1(sid, ids, is_skin=False, comment='set card')
        model.setup()

        log = SimpleLogger(level='warning')
        #model = BDF(log=log)
        #add_methods = model._add_methods
        #add_methods._add_coord_object(coord)
        #add_methods._add_caero_object(caero2)
        #add_methods._add_set_object(set_obj)
        model.add_grid(7, [7., 0., 0.], cp=0, cd=0, ps='', seid=0, comment='')
        model.add_grid(13, [13., 0., 0.], cp=0, cd=0, ps='', seid=0, comment='')
        #add_methods._add_node_object(grid7)
        #add_methods._add_node_object(grid13)

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
        model.spline2.add_card(bdf_card, comment='spline2_a')
        model.setup()
        #spline_a.write()
        #spline_a.raw_fields()

        spline_b = model.add_spline2(eid, caero, id1, id2, setg,
                                     dz=dz,
                                     dtor=dtor, cid=cid,
                                     #dthx=dthx,
                                     #dthy=dthy,
                                     usage=usage, comment='spline2_b')
        model.setup()
        #spline_b.validate()
        spline2.write()

        #model.add_set1(sid+1, ids, is_skin=True)
        save_load_deck(model, run_test_bdf=False, run_renumber=False)

    def test_caero2_1(self):
        """checks the CAERO2/PAERO2/AERO/AEFACT card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        caero2 = model.caero2
        eid = 1
        pid = 10
        cp = 4
        nsb = 0
        nint = 0

        lsb = 3
        lint = 6
        igroup = 1
        p1 = [0., 1., 2.]
        x12 = 10.
        caero2.add_card(BDFCard(['CAERO2', eid, pid, cp, nsb, nint,
                                 lsb, lint, igroup, ] + p1 + [x12]))

        #---------------
        # nsb=lsb=None=0
        caero2.add(eid, pid, igroup, p1, x12,
                    cp=cp, nsb=None, nint=None, lsb=None, lint=6,
                    comment='this is a caero')
        caero2.parse_cards()
        #model.setup()
        with self.assertRaises(ValueError):
            caero2.validate()

        # nint=lint=None=0
        caero2.clear()
        caero2.add(eid, pid, igroup, p1, x12,
                   cp=cp, nsb=3, nint=None, lsb=3, lint=None,
                   comment='this is a caero')
        caero2.parse_cards()
        with self.assertRaises(ValueError):
            caero2.validate()

        # they're all bad?
        caero2.clear()
        caero2.add(eid, pid, igroup, p1, x12,
                    cp=cp, nsb=0, nint=0, lsb=0, lint=0,
                    comment='this is a caero')
        caero2.parse_cards()
        with self.assertRaises(ValueError):
            caero2.validate()

        #---------------
        caero2a = model.add_caero2(eid, pid, igroup, p1, x12, cp=4, nsb=0, nint=0,
                                   lsb=3, lint=6, comment='this is a caero')
        caero2.clear()

        eid = 200
        caero2b = model.add_caero2(eid, pid, igroup, p1, x12, cp=4, nsb=10, nint=7,
                                   lsb=0, lint=0, comment='this is a caero')
        caero2.parse_cards()
        caero2.validate()
        str(caero2)
        #str(caero2f.write_card())

        aefact = AEFACT(model)
        #aefact = model.aefact
        aefact.add_card(BDFCard(['AEFACT', lint, 0., 1., 2., 3., 4., 5.]))
        lrib = 3
        aefact.add_card(BDFCard(['AEFACT', lrib, 0., 1., 2.]))
        aefact.validate()
        aefact.write()
        model.aefact = aefact
        model.setup()

        orient = 'Z'
        width = 10.
        AR = 2.
        lrsb = 0 # aefact
        lrib = 3 # aefact
        lth1 = 0
        lth2 = 0
        lth = [lth1, lth2]
        thi = [0, 0, 0]
        thn = [0, 0, 0]
        paero2a = PAERO2(model)
        paero2a.add_card(BDFCard(['PAERO2', pid, orient, width, AR,
                                  lrsb, lrib] + lth + thi + thn),
                         comment='paero2')
        paero2a.parse_cards()
        paero2a.validate()
        model.add_paero2(pid, orient, width, AR, thi, thn,
                         lrsb=0, lrib=3, lth=lth, comment='paero2')
        paero2b = model.paero2
        paero2b.parse_cards()
        paero2b.validate()

        pid = 42
        model.add_paero2(pid, orient, width, AR, thi, thn,
                         lrsb=None, lrib=None, lth=None, comment='')
        paero2c = model.paero2
        paero2c.parse_cards()
        paero2c.validate()
        paero2b.write()
        #model.paeros[pid] = paero

        model.add_cord2r(
            cp, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        model.coord.validate()
        #model.coords[cp] = coord

        #model.add_aefact(lrib, [0., 1., 2., 3., 4., 5.])
        aefact = model.aefact
        aefact.validate()

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

        #paero2b.raw_fields()
        #caero2f.raw_fields()
        #model.uncross_reference()
        model.write_bdf('aero.temp2')

        model.cross_reference()
        model.write_bdf('aero.temp3')

        #caero2f.raw_fields
        paero2 = model.paero2
        caero2.get_points_elements_3d()
        caero2.get_points()
        unused_xyz, unused_elems = caero2.get_points_elements_3d()


        #caero2g.get_points()
        #caero2g.get_points_elements_3d()
        #unused_xyz, unused_elems = caero2g.get_points_elements_3d()

        #model.uncross_reference()
        #model.safe_cross_reference()
        #model.uncross_reference()
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
        caero2i = model.add_caero2(
            eid, pid, igroup, p1, x12,
            cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
            comment='this is a caero')
        caero2.validate()
        model.setup()
        #caero2.cross_reference(model)
        caero2.write()

        #model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()

        caero2_set_points = model.add_caero2(
            eid, pid, igroup, p1, x12,
            cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint)
        p1 = [0., 0., 0.]
        p2 = [1., 2., 3.]
        #caero2_set_points.set_points([p1, p2])
        #assert np.allclose(caero2_set_points.x12, 1.), caero2_set_points.x12
        save_load_deck(model)

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
        caero3 = model.caero3

        model.add_cord2r(
            cp, origin=[0., 0., 0.],
            zaxis=[0., 0., 1.],
            xzplane=[1., 0., 0.], rid=0)
        origin = None
        zaxis = None
        xzplane = None
        model.add_cord2r(cp, origin, zaxis, xzplane, rid=0, comment='cord2r')
        model.setup()
        model.coord.validate()
        #model.coords[cp] = coord

        model.add_paero3(pid, nbox, ncontrol_surfaces, x, y,
                         comment='paero3')
        model.setup()
        model.paero3.validate()
        #paero3.raw_fields()

        card = ['CAERO3', 2000, 20001, 0, 22, 33, None, None, None,
                1.0, 0.0, 0., 100., 17., 130., 0., 100.]
        bdf_card = BDFCard(card, has_none=True)
        caero3.add_card(bdf_card, comment='msg')
        model.setup()
        caero3.validate()
        caero3.write()
        #caero3a.raw_fields()
        caero3.clear()

        caero3i = model.add_caero3(eid, pid,
                                   p1, x12, p4, x43,
                                   cp=cp, list_w=list_w, list_c1=list_c1, list_c2=list_c2,
                                   comment='caero3')
        caero3.validate()
        model.setup()

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
        model.setup()
        model.validate()

        caero3.write()
        #caero3b.cross_reference(model)
        #caero3b.write()
        #caero3a.raw_fields()
        #caero3b.uncross_reference()
        #caero3b.write()
        #caero3a.raw_fields()

        #xref_errors = defaultdict(list)
        #caero3.safe_cross_reference(model, xref_errors)

        npoints_expected = 33
        nelements_expected = 20
        npoints, nelements = caero3.get_npanel_points_elements()
        assert npoints == npoints_expected
        assert nelements == nelements_expected

        caero3.get_points()
        caero3.panel_points_elements()
        if IS_MATPLOTLIB:
            fig, ax = _setup_aero_plot()
            caero3.plot(ax)
            fig.show()

        model.get_bdf_stats()
        save_load_deck(model, run_mirror=False)

    def test_paero3(self):
        """checks the PAERO3"""
        # +--------+------+------+-------+------+-----+------+------+------+
        # |    1   |   2  |   3  |   4   |   5  |  6  |   7  |   8  |  9   |
        # +========+======+======+=======+======+=====+======+======+======+
        # | PAERO3 |  PID | NBOX | NCTRL |      |  X5 |  Y5  |  X6  |  Y6  |
        # +--------+------+------+-------+------+-----+------+------+------+
        # |        |  X7  |  Y7  |   X8  |  Y8  |  X9 |  Y9  |  X10 |  Y10 |
        # +--------+------+------+-------+------+-----+------+------+------+
        # |        |  X11 |  Y11 |  X12  |  Y12 |     |      |      |      |
        # +--------+------+------+-------+------+-----+------+------+------+
        # | PAERO3 | 2001 |  15  |   1   |      | 0.  |  65. |      |      |
        # +--------+------+------+-------+------+-----+------+------+------+
        # |        |  78. |  65. |  108. |  65. | 82. | 97.5 | 112. | 97.5 |
        # +--------+------+------+-------+------+-----+------+------+------+
        # |        |  86. | 130. |  116. | 130. |     |      |      |      |
        # +--------+------+------+-------+------+-----+------+------+------+
        fields = ['PAERO3', 2001, 15, 1, None, 0., 65., None, None,
                  78., 65., 108., 65., 82., 97.5, 112., 97.5,
                  86., 130., 116., 130.]
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        model.add_card(fields, fields[0])
        model.setup()
        paero3 = model.paero3.slice_card_by_property_id(2001)
        #assert paero.npoints == 8, paero.npoints
        #paero.raw_fields()
        paero3.write()

    def test_paero4(self):
        """checks the PAERO4"""
        # +--------+------+-------+--------+-------+-------+--------+--------+--------+
        # |    1   |   2  |   3   |   4    |   5   |   6   |    7   |   8    |    9   |
        # +========+======+=======+========+=======+=======+========+========+========+
        # | PAERO4 | PID  | CLA   |  LCLA  |  CIRC | LCIRC |  DOC1  |  CAOC1 | GAPOC1 |
        # +--------+------+-------+--------+-------+-------+--------+--------+--------+
        # |        | DOC2 | CAOC2 | GAPOC2 |  DOC3 | CAOC3 | GAPOC3 |  etc.  |        |
        # +--------+------+-------+--------+-------+-------+--------+--------+--------+
        # | PAERO4 | 6001 |   1   |   501  |   0   |   0   |   0.0  |   0.0  |   0.0  |
        # +--------+------+-------+--------+-------+-------+--------+--------+--------+
        # |        | 0.50 |  0.25 |  0.02  |  0.53 |  0.24 |   0.0  |        |        |
        # +--------+------+-------+--------+-------+-------+--------+--------+--------+
        pid = 6001
        cla = 1
        lcla = 501
        circ = 0
        lcirc = 0
        dcg1 = [0., 0., 0.]
        dcg2 = [0.5, 0.25, 0.02]
        dcg3 = [0.53, 0.24, 0.]
        card = ['PAERO4', pid, cla, lcla, circ, lcirc] + dcg1 + dcg2 + dcg3

        model = BDF(debug=False)
        bdf_card = BDFCard(card, has_none=True)
        paero4 = model.paero4.add_card(bdf_card, comment='msg')
        str(paero4)
        #paero4.cross_reference(None)
        save_load_deck(model)

    def test_caero4_1(self):
        """checks the CAERO4/PAERO4"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        pid = 1001
        docs = []
        caocs = []
        gapocs = []
        model.add_paero4(pid, docs, caocs, gapocs,
                                  cla=0, lcla=0, circ=0, lcirc=0,
                                  comment='paero4')
        paero4 = model.paero4
        paero4.validate()
        #paero4.raw_fields()

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
        model.add_card(bdf_card, 'CAERO4', comment='msg')
        caero4a = model.caero4

        model.setup()
        caero4a.validate()
        npoints, nelements = caero4a.get_npanel_points_elements()
        assert npoints == 10, npoints
        assert nelements == 4, nelements
        caero4a.write()
        #caero4a.raw_fields()

        #caero4a.cross_reference(model)
        #points, elements = caero4a.panel_points_elements()
        #del points, elements

        p1 = [x1, y1, z1]
        p4 = [x4, y4, z4]
        model.add_caero4(eid, pid, p1, x12, p4, x43,
                                   cp=cp, nspan=nspan, lspan=lspan,
                                   comment='caero4b')
        caero4b = model.caero4
        caero4b.validate()
        caero4b.write()
        #caero4b.raw_fields()

        #caero4b.cross_reference(model)
        #caero4b.write()
        #caero4b.raw_fields()
        points, elements = caero4b.panel_points_elements()
        del points, elements

        p1, unused_p2, unused_p3, p4 = caero4b.get_points()

        caero4c = model.caero4
        caero4c.add(eid, pid, p1, x12, p4, x43,
                    cp=0, nspan=0, lspan=0,
                    comment='caero4c')
        caero4c.parse_cards()
        with self.assertRaises(RuntimeError):
            # nspan=lspan=0
            caero4c.validate()


        model.setup()
        #model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)
        model2 = read_bdf(bdf_filename, xref=False, punch=True, debug=False)
        #model.safe_cross_reference()
        #model2.safe_cross_reference()
        save_load_deck(model)

    def test_caero5_1(self):
        """checks the CAERO5/PAERO5"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        pid = 6001
        caoci = [0., 0.5, 1.0]
        caero5 = model.caero5
        paero5 = model.paero5
        model.add_paero5(pid, caoci,
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
        #paero5 = model.paeros[pid]
        #paero5.raw_fields()

        # ---------------------------------------------
        model = BDF(debug=None)
        paero5 = model.paero5
        model.add_paero5(pid, caoci, nalpha=0, lalpha=0, nxis=0, lxis=0,
                         ntaus=0, ltaus=0, comment='paero5')
        model.setup()
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
        caero5i = model.add_caero5(eid, pid, p1, x12, p4, x43,
                                   cp=0, nspan=5, lspan=0,
                                   ntheory=0, nthick=0,
                                   comment='msg')
        model.setup()
        model.validate()

        lxis = 43
        paero5.lxis[0] = lxis
        aefact_sid = lxis
        Di = [0., 0.5, 1.]
        model.add_aefact(aefact_sid, Di, comment='aefact')

        ltaus = 44
        paero5.ltaus[0] = ltaus
        aefact_sid = ltaus
        Di = [0., 0.5, 1.]
        model.add_aefact(aefact_sid, Di, comment='aefact2')

        #caero5.cross_reference(model)
        model.cross_reference()
        unused_npoints, unused_nelements = model.caero5.get_npanel_points_elements()
        unused_points, unused_elements = model.caero5.panel_points_elements()
        caero5.write()
        #caero5.raw_fields()

        #model.uncross_reference()
        caero5.write()
        model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        bdf_filename.seek(0)

        read_bdf(bdf_filename, xref=False, punch=True, debug=False)
        #model.safe_cross_reference()
        save_load_deck(model, run_renumber=False, run_test_bdf=False)


        #caero5.raw_fields()


   # def test_paero1_1(self):
   # def test_paero2_1(self):
   # def test_paero3_1(self):
   # def test_paero4_1(self):
   # def test_paero5_1(self):

   # def test_spline1_1(self):
   # def test_spline2_1(self):
    def test_spline3(self):
        """checks the SPLINE3 card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        spline3 = model.spline3
        eid = 100
        pid = 10
        igid = 1
        p1 = [0., 0., 0.]
        x12 = x43 = 3.
        p4 = [1., 11., 1.]
        caero = eid
        model.add_caero1(eid, pid, igid, p1, x12, p4, x43,
                         cp=0,
                         nspan=5, lspan=0,
                         nchord=5, lchord=0, comment='')
        model.add_paero1(pid, caero_body_ids=None, comment='')
        model.add_grid(5, [0., 0., 0.])

        spline_id = 101
        box_id = 42
        components = 3
        nids = [5]
        displacement_components = [3]
        coeffs = [1.0]
        model.add_spline3(
            spline_id, caero, box_id, components, nids,
            displacement_components, coeffs, usage='BOTH', comment='spline3')
        spline3 = model.spline3
        model.setup()
        spline3.validate()
        spline3.write()
        #spline3.raw_fields()

        spline_id = 102
        nids = [5, 6]
        displacement_components = [3, 6]
        coeffs = [1.0, 2.0]
        model.add_spline3(
            spline_id, caero, box_id, components, nids,
            displacement_components, coeffs, usage='failed', comment='spline3')
        cref = bref = sref = 1.0
        model.add_aeros(cref, bref, sref)
        model.setup()

        spline3b = spline3.slice_card_by_id(spline_id)
        spline3b.write()
        with self.assertRaises(RuntimeError):
            spline3b.validate()
        spline3b.usage = 'BOTH'
        np.append(spline3b.displacement_components, 1)
        np.append(spline3b.coeffs, 0.1)
        spline3b.validate()
        print(model.spline3.write())

        #del model.splines[spline_id]
        #model._type_to_id_map['SPLINE3'].remove(spline_id)
        model.validate()


        #model.add_grid(5, [0., 0., 0.])
        model.add_grid(6, [1., 0., 0.])
        model.add_grid(7, [1., 1., 0.])
        model.add_grid(42, [0., 1., 0.])

        #spline3.cross_reference(model)
        model.cross_reference()
        #model.pop_xref_errors()
        spline3.write()
        #spline3.raw_fields()
        save_load_deck(model, run_renumber=False)
        #spline3b.eid = 1000

        #spline3b = model.splines[102]
        #del model.splines[102]
        #model._type_to_id_map['SPLINE3'].remove(102)
        #model.splines[1000] = spline3b
        print(spline3b.spline_id)
        spline3b.spline_id[0] = 1000
        spline3.spline_id[1] = 1000
        print(spline3b.spline_id)
        #spline3b.sort()
        #print(model.splines)

        #model.uncross_reference()
        #model.add_grid(42, [0., 1., 0.])
        #spline3b.nodes.append(42)
        #spline3b.displacement_components.append(4)
        #spline3b.coeffs.append(0.5)
        #spline3b.validate()
        model.cross_reference()

        eid = 100
        pid = 100
        mid = 100
        t = 0.1
        G11 = G12 = G13 = G22 = G23 = G33 = 1.0
        model.add_cquad4(eid, pid, [5, 6, 7, 42])
        model.add_pshell(pid, mid1=mid, t=t)
        model.add_mat2(mid, G11, G12, G13, G22, G23, G33,
                       rho=0., a1=None, a2=None, a3=None, tref=0., ge=0.,
                       St=None, Sc=None, Ss=None, mcsid=None, comment='')

        #spline3b.comment = ''
        model.setup()

        print(spline3.write())
        spline3b = spline3.slice_card_by_id(1000, assume_sorted=False)
        print(spline3.write())

        RUN_SPLINE3 = False
        if RUN_SPLINE3:
            spline3.clear()
            lines = spline3b.write().rstrip().split('\n')
            model.add_card(lines, 'SPLINE3', is_list=False)
            model.setup()
            print(spline3.write())
            #spline = model.splines[1000]
            #spline = model.spline3
            assert np.array_equal(spline3b.nodes[0], [5, 6, 7, 42]), spline3b.nodes

        model.cross_reference()
        #model.pop_xref_errors()
        model.validate()
        #spline3.raw_fields()

        # have one earlier...
        #save_load_deck(model, run_renumber=False)

    def test_spline4(self):
        """checks the SPLINE4 card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        eid = 1
        caero = 10
        aelist = 11
        setg = 12
        dz = 0.
        method = 'TPS'
        usage = 'FORCE'
        nelements = 4
        melements = 5
        spline = model.add_spline4(eid, caero, aelist, setg, dz, method, usage,
                                   nelements, melements, comment='spline4')
        spline = model.spline4
        #spline.raw_fields()

        elements = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        model.add_aelist(aelist, elements)

        paero = 20
        igid = 42
        p1 = [0., 0., 0.]
        x12 = 10.
        p4 = [0., 10., 0.]
        x43 = 3.
        model.add_caero1(caero, paero, igid, p1, x12, p4, x43, cp=0, nspan=5,
                         lspan=0, nchord=10, lchord=0,
                         comment='')
        model.add_paero1(paero)

        velocity = None
        cref = 1.0
        rho_ref = 1.0
        model.add_aero(velocity, cref, rho_ref,
                       comment='')

        model.add_set1(setg, [1, 2, 3])
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])

        eid = 2
        setg = 13
        ids = [1, 2, 3]
        model.add_set1(setg, ids)
        model.add_spline4(eid, caero, aelist, setg, dz, method, usage,
                          nelements, melements, comment='spline4')
        #spline = model.splines[eid]
        #del model.splines[eid]
        #spline.cross_reference(model)
        model.setup()

        model.pop_parse_errors()
        #model.pop_xref_errors()
        model.validate()
        save_load_deck(model)

    def test_spline5(self):
        """checks the SPLINE5 card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        eid = 1
        caero = 10
        aelist = 11
        setg = 12
        thx = 7.
        thy = 8.
        #dz = 0.
        #method = 'cat'
        #usage = 'dog'
        #nelements = 4
        #melements = 5
        #dtor = 47
        spline = model.add_spline5(eid, caero, aelist, setg, thx, thy, dz=0., dtor=1.0,
                                   cid=0, usage='BOTH', method='BEAM', ftype='WF2',
                                   rcore=None, comment='spline5')
        #spline.raw_fields()

        elements = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        model.add_aelist(aelist, elements)

        paero = 20
        igid = 42
        p1 = [0., 0., 0.]
        x12 = 10.
        p4 = [0., 10., 0.]
        x43 = 3.
        model.add_caero1(caero, paero, igid, p1, x12, p4, x43, cp=0, nspan=5,
                         lspan=0, nchord=10, lchord=0,
                         comment='')
        model.add_paero1(paero)

        velocity = None
        cref = 1.0
        rho_ref = 1.0
        model.add_aero(velocity, cref, rho_ref,
                       comment='')

        model.add_set1(setg, [1, 2, 3])
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        model.setup()
        model.pop_parse_errors()
        #model.pop_xref_errors()
        model.validate()
        #model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_aesurf_1(self):
        """checks the AESURF/AELIST cards"""
        aesid = 10
        label = 'FLAP'
        cid1 = 0
        aelist_id1 = 10
        cid2 = None
        aelist_id2 = None
        model = BDF(debug=True, log=None, mode='msc')
        aesurf = model.aesurf
        aesurf1 = aesurf.add(aesid, label, cid1, aelist_id1, cid2, aelist_id2,
                             #eff, ldw,
                             #crefc, crefs, pllim, pulim,
                             #hmllim, hmulim, tqllim, tqulim,
                             comment='aesurf comment')
        aesurf2 = aesurf.add_card(BDFCard(
            [
                'AESURF', aesid, label, cid1, aelist_id1, cid2, aelist_id2,
                #eff, ldw,
                #crefc, crefs, pllim, pulim,
                #hmllim, hmulim, tqllim, tqulim,
            ]), comment='aesurf comment')
        #assert aesurf1 == aesurf2

        cid2 = 1
        coord = model.add_cord2r(cid2, rid=0, origin=[0., 0., 0.],
                                 zaxis=[1., 0., 0.], xzplane=[0., 0., 1.], comment='')

        aelist_id1 = 10
        aelist_id2 = 20
        aesurf2 = model.add_aesurf(
            aesid, label, cid1, aelist_id1, cid2, aelist_id2,
            #eff, ldw,
            #crefc, crefs, pllim, pulim,
            #hmllim, hmulim, tqllim, tqulim,
            comment='aesurf comment')
        aesurf1 = aesurf2 = model.aesurf

        aesurf1.validate()
        aesurf2.validate()
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        add_methods = model._add_methods

        elements = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        unused_aelist = model.add_aelist(aelist_id1, elements, comment='')

        elements = [11, 22, 33, 44, 55, 66, 77, 88, 99]
        unused_aelist = model.add_aelist(aelist_id2, elements, comment='')

        aesid += 1
        model.add_aesurf(
            aesid, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
            eff=1.0, ldw='LDW',
            crefc=1.0, crefs=1.2,
            pllim=-np.pi/2, pulim=np.pi/2.,
            hmllim=-42., hmulim=42.,  # hinge moment limits in force/disp
            tqllim=10, tqulim=11,  # TABLEDi deflection limits vs. dynamic pressure
        )
        # lower
        table_id = 10
        x = np.linspace(0.1, 1.)
        y = np.log(np.linspace(1.1, 2.))[::-1]
        model.add_tabled1(table_id, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')

        # upper
        table_id = 11
        y2 = -y
        model.add_tabled1(table_id, x, y2, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')

        model.setup(run_geom_check=True)

        #aesurf1.cross_reference(model)
        aesurf1.write()
        #aesurf1.raw_fields()
        #aesurf1.uncross_reference()
        aesurf1.write()
        #aesurf1.cross_reference(model)
        #aesurf1.raw_fields()

        #aesurf2.cross_reference(model)
        aesurf2.write()
        #aesurf2.raw_fields()
        #aesurf2.uncross_reference()
        aesurf2.write()
        #aesurf2.cross_reference(model)
        #aesurf2.raw_fields()
        #model.cross_reference()
        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_flfact(self):
        """checks the FLFACT card"""
        #list[f1, THRU, fnf, nf, fmid]
            #f1 : float
                #first value
            #THRU : str
                #the word THRU
            #nf : float
                #second value
            #fmid : float; default=(f1 + fnf) / 2.
                #the mid point to bias the array

        sid = 42
        factors = [0.200, 'THRU', 0.100, 11, 0.1333]
        flfact = FLFACT(sid, factors)
        assert len(flfact.factors) == 11
        #print(flfact)

        factors = [0.200, 'THRU', 0.100, 11]
        flfact = FLFACT(sid, factors)
        assert len(flfact.factors) == 11

    def test_flutter(self):
        """checks the FLUTTER/FLFACT cards"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        flutter = model.flutter
        flfact = model.flfact

        sid = 75
        method = 'PKNL'
        idensity = 76
        imach = 77
        ivelocity = 78

        # density, mach, velocity
        flutter1 = model.add_flutter(sid, method, idensity, imach, ivelocity,
                                     imethod='L', nvalue=None,
                                     omax=None, epsilon=1.0e-3)
        card = BDFCard(['FLUTTER', sid, method, idensity, imach, ivelocity])
        flutter2 = FLUTTER.add_card(card, comment='flutter card')
        assert flutter2.headers == ['density', 'mach', 'velocity'], flutter2.headers

        #assert flutter1.get_field(1) == sid, flutter1.get_field(1)
        #assert flutter1.get_field(2) == 'PKNL', flutter1.get_field(2)
        #assert flutter1.get_field(3) == idensity, flutter1.get_field(3)
        #assert flutter1.get_field(4) == imach, flutter1.get_field(4)
        #assert flutter1.get_field(5) == ivelocity, flutter1.get_field(5)
        #assert flutter1.get_field(6) == 'L', flutter1.get_field(6)
        #assert flutter1.get_field(7) is None, flutter1.get_field(7)
        #assert flutter1.get_field(8) == 1.0e-3, flutter1.get_field(8)
        #with self.assertRaises(KeyError):
            #assert flutter1.get_field(9) == 1.0e-3, flutter1.get_field(9)
        #flutter1.validate()
        #flutter1.write()
        #flutter2.validate()
        #flutter2.write()
        flutter.write()

        densities = np.linspace(0., 1.)
        unused_density = model.add_flfact(idensity, densities)

        machs = np.linspace(0.7, 0.8)
        mach = FLFACT(imach, machs)
        card = BDFCard(['FLFACT', imach] + list(machs))
        FLFACT.add_card(card, comment='flfact card')
        model.add_flfact(imach, machs, comment='flfact')
        model.setup(run_geom_check=True)
        flfact.write(size=16)
        flfact.write(size=8)

        velocities = np.linspace(3., 4.)
        model.add_flfact(ivelocity, velocities)
        model.setup(run_geom_check=True)
        velocity = flfact.slice_card_by_flfact_id(ivelocity)
        #velocity.validate()
        velocity.write()
        #assert velocity.min() == 3., velocities
        #assert velocity.max() == 4., velocities
        #model.flfacts[ivelocity] = velocity

        ikfreq = 79
        kfreqs = np.linspace(0.1, 0.2)
        card = ['FLFACT', ikfreq] + list(kfreqs)
        model.add_card(card, card[0])
        model.setup(run_geom_check=True)
        #kfreq = model.FLFACT(ikfreq)
        kfreq = flfact.slice_card_by_flfact_id(ikfreq)
        model.setup(run_geom_check=True)
        #kfreq.validate()
        kfreq.write()
        #assert kfreq.min() == 0.1, kfreqs
        #assert kfreq.max() == 0.2, kfreqs
        #model.flfacts[ikfreq] = kfreq

        ikfreq2 = 80
        card = ['FLFACT', ikfreq2, 10., 'THRU', 20., 11]
        model.add_card(card, card[0])
        model.setup(run_geom_check=True)
        #kfreq = model.FLFACT(ikfreq2)
        kfreq = flfact.slice_card_by_flfact_id(ikfreq2)
        #kfreq.validate()
        kfreq.write()
        #assert kfreq.min() == 10., 'min=%s; card=%s factors=%s' % (kfreq.min(), card, kfreq.factors)
        #assert kfreq.max() == 20., 'max=%s; card=%s factors=%s' % (kfreq.max(), card, kfreq.factors)
        #model.flfacts[ikfreq] = kfreq

        ikfreq3 = 81
        factors = [10., 'THRU', 20., 10]
        kfreq = model.add_flfact(ikfreq3, factors)
        model.setup(run_geom_check=True)
        kfreq = flfact.slice_card_by_flfact_id(ikfreq3)
        #kfreq.validate()
        kfreq.write()
        #assert kfreq.min() == 10., 'min=%s; factors=%s' % (kfreq.min(), factors)
        #assert kfreq.max() == 20., 'max=%s; factors=%s' % (kfreq.max(), factors)
        #model.flfacts[ikfreq] = kfreq
        #kfreq.validate()

        #ikfreq4 = 82
        #kfreq2 = model.add_flfact(ikfreq4, [])
        #with self.assertRaises(ValueError):
            #kfreq2.validate()
        #kfreq2.factors = [1.]
        #kfreq2.validate()
        #kfreq2.write_card()

        # density, mach, rfreq
        card = ['FLUTTER', 85, 'KE', idensity, imach, ikfreq]
        model.add_card(card, card[0])
        model.setup(run_geom_check=True)

        #model.pop_parse_errors()
        model.cross_reference()
        #model.pop_xref_errors()

        #flutter = model.Flutter(85)
        #assert flutter.headers == ['density', 'mach', 'reduced_frequency'], flutter.headers
        flutter.write()
        #flutter.raw_fields()

        #model.uncross_reference()
        #model.safe_cross_reference()
        save_load_deck(model)

    def test_flutter_2(self):
        """validates the FLUTTER card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        flutter = model.flutter
        flfact = model.flfact

        method = 'TEST'
        imethod = 'TEST2'
        sid = 1
        idensity = 10
        imach = 20
        ivelocity = 30
        with self.assertRaises(AssertionError):
            unused_flutter_id = model.add_flutter(
                sid, method, idensity, imach, ivelocity,
                imethod=imethod, nvalue=None,
                omax=None, epsilon=1.0e-3)
        #flutter.parse_cards()
        #with self.assertRaises(ValueError):
            #flutter.validate()

    def test_flutter_3(self):
        """tests the flutter sweeps"""
        alts = np.linspace(-10000., 50000.)[::-1]

        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        sid = 70
        method = 'PKNL'
        density = 71
        mach = 72
        reduced_freq_velocity = 73
        flutter = model.add_flutter(sid, method, density, mach, reduced_freq_velocity)
        RUN_FLUTTER = False
        if RUN_FLUTTER:
            flutter.make_flfacts_alt_sweep_constant_mach(
                model, 0.7, alts, eas_limit=1000.0, alt_units=u'ft',
                velocity_units=u'in/s', density_units=u'slinch/in^3',
                eas_units=u'ft/s')

        sid = 80
        density = 81
        mach = 82
        reduced_freq_velocity = 83
        flutter = model.add_flutter(sid, method, density, mach, reduced_freq_velocity)

        if RUN_FLUTTER:
            alt = 10000.
            machs = np.arange(0.1, 0.8)
            flutter.make_flfacts_mach_sweep(
                model, alt, machs, eas_limit=1000., alt_units='m',
                velocity_units='m/s',
                density_units='kg/m^3',
                eas_units='m/s')

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
        # ----------------------------------------------------------------------
        model = BDF(debug=False)

        machs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        reduced_freqs = [0.01, 0.02, 0.03]
        mkaero = model.add_mkaero1(machs, reduced_freqs, comment='mkaero1')
        mkaero.raw_fields()
        msg = mkaero.write_card()
        lines = msg.strip().split('\n')
        for line1, line2 in zip(lines, expected2):
            msg = '\nline    =%r\n' %  str(line1)
            msg += 'expected=%r\n%s' % (str(line2), msg)
            assert line1 == line2, msg

        mkaerob = MKAERO1([], reduced_freqs)
        with self.assertRaises(ValueError):
            mkaerob.validate()
        with self.assertRaises(ValueError):
            mkaerob.write_card()

        mkaeroc = MKAERO1([0.1, 0.2], [])
        with self.assertRaises(ValueError):
            mkaeroc.validate()
        with self.assertRaises(ValueError):
            mkaeroc.write_card()

        machs = [0.01, 0.02, 0.03]
        reduced_freqs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        mkaero = model.add_mkaero1(machs, reduced_freqs, comment='mkaero1')

        machs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
        reduced_freqs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        mkaero = model.add_mkaero1(machs, reduced_freqs, comment='mkaero1')

        machs_spike = [0.01]
        reduced_freqs_spike = [0.1, 0.2]
        mkaero = model.add_mkaero1(machs_spike, reduced_freqs_spike, comment='mkaero1')

        # TODO: this fails...because it's linked to the first card somehow
        #mkaerod = model.add_mkaero1(machs, [])
        #with self.assertRaises(ValueError):
            #mkaerod.write_card()
        save_load_deck(model)

    def test_mkaero2(self):
        """checks the MKAERO2 card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
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

        # at least one mach
        machs = []
        reduced_freqs = [42.]
        mkaero = MKAERO2(machs, reduced_freqs)
        with self.assertRaises(ValueError):
            mkaero.validate()

        # at least one rfreq
        machs = [0.8]
        reduced_freqs = []
        mkaero = MKAERO2(machs, reduced_freqs)
        with self.assertRaises(ValueError):
            mkaero.validate()

        # should be the same length
        machs = [0.8]
        reduced_freqs = [42., 43.]
        mkaero = MKAERO2(machs, reduced_freqs)
        with self.assertRaises(ValueError):
            mkaero.validate()

        # split the write card method
        machs = [0.1, 0.2, 0.3, 0.4, 0.5]
        reduced_freqs = [1., 2., 3., 4., 5.]
        mkaero = MKAERO2(machs, reduced_freqs)
        mkaero.validate()
        mkaero.write_card()

        mkaerob = model.add_mkaero2([], reduced_freqs)
        with self.assertRaises(ValueError):
            mkaerob.validate()
        with self.assertRaises(ValueError):
            mkaerob.write_card()

        mkaeroc = model.add_mkaero2([0.1, 0.2], [])
        with self.assertRaises(ValueError):
            mkaeroc.validate()
        with self.assertRaises(ValueError):
            mkaeroc.write_card()

        mkaeroc = model.add_mkaero2([], [])
        with self.assertRaises(ValueError):
            mkaeroc.validate()
        with self.assertRaises(ValueError):
            mkaeroc.write_card()

    def _test_diverg(self):
        """checks the DIVERG card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        sid = 100
        nroots = 21
        machs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        diverg = DIVERG(sid, nroots, machs, comment='divergence')
        diverg.validate()
        diverg.write_card()

        diverg = model.add_card(['DIVERG', sid, nroots] + machs, 'DIVERG', comment='divergence')
        model.validate()
        save_load_deck(model)
        #diverg.validate()
        #diverg.write_card()

    def test_trim_01(self):
        """checks the TRIM card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        #model.add_aecompl

        sid = 100
        mach = 0.75
        q = 100.
        labels = ['ALPHA', 'ALPHA']
        uxs = [10., 20.]
        trim1 = TRIM(model)
        trim2 = TRIM(model)
        trim1.add(sid, mach, q, labels, uxs, trim_type=1)
        trim2.add(sid+1, mach, q, labels, uxs, trim_type=2)
        trim1.parse_cards()
        trim2.parse_cards()
        with self.assertRaises(RuntimeError):
            trim1.validate()
        with self.assertRaises(RuntimeError):
            trim2.validate()

        labels = ['ALPHA']
        uxs = [10., 20.]
        trim1 = TRIM(model)
        trim2 = TRIM(model)
        with self.assertRaises(AssertionError):  # len(labels) != len(uxs)
            trim1.add(sid, mach, q, labels, uxs, trim_type=1)
        with self.assertRaises(AssertionError):  # len(labels) != len(uxs)
            trim2.add(sid, mach, q, labels, uxs, trim_type=2)

        labels = ['ALPHA', 'BETA']
        uxs = [10., 20.]
        trim1 = TRIM(model)
        trim2 = TRIM(model)
        trim1.add(sid, mach, q, labels, uxs, trim_type=1)
        trim2.add(sid, mach, q, labels, uxs, trim_type=2)
        trim1.parse_cards()
        trim2.parse_cards()
        trim1.validate()
        trim2.validate()
        trim1.write()
        trim2.write()

        labels = ['ALPHA']
        uxs = [10.]
        trim1 = TRIM(model)
        trim2 = TRIM(model)
        trim1.add(sid, mach, q, labels, uxs, aeqr=3.0, comment='', trim_type=1)
        trim2.add(sid, mach, q, labels, uxs, aeqr=3.0, comment='', trim_type=2)
        trim1.parse_cards()
        trim2.parse_cards()
        trim1.validate()
        trim2.validate()
        trim1.write()
        trim2.write()

        labels = ['ALPHA', 'BETA']
        uxs = [10., 20.]
        trim1 = TRIM(model)
        trim2 = TRIM(model)
        trim1.add(sid, mach, q, labels, uxs, aeqr=3.0, comment='', trim_type=1)
        trim2.add(sid, mach, q, labels, uxs, aeqr=3.0, comment='', trim_type=1)
        trim1.parse_cards()
        trim2.parse_cards()
        trim1.validate()
        trim2.validate()
        trim1.write()
        trim2.write()

        model.add_card(['TRIM', sid, mach, q, labels[0], uxs[0]], 'TRIM', comment='$ trim')
        model.setup()
        model.validate()
        model._verify_bdf(xref=False)
        save_load_deck(model)

    def _test_trim_02(self):
        """checks the TRIM card with a 2.5g pullup"""
        model = BDF()
        sid = 75
        mach = 0.75
        q = 100.
        labels = ['NZ']
        uxs = [2.5]
        trim1 = model.add_trim(sid, mach, q, labels, uxs, aeqr=0.0, comment='')
        trim1.parse_cards()
        trim1.validate()

        trim2 = model.add_trim(sid+1, mach, q, labels, uxs, aeqr=0.0, trim_type=2, comment='')
        trim2.parse_cards()
        trim2.validate()

    def test_trim_03(self):
        """checks the TRIM card with a 2.5g pullup"""
        model = BDF(debug=None)
        sid = 75
        mach = 0.75
        q = 100.
        labels = ['URDD3', 'PITCH']
        uxs = [2.5, 0.0]
        trim1a = TRIM(model)
        trim2a = TRIM(model)
        trim1a.add(sid, mach, q, labels, uxs, aeqr=0.0,
                   trim_type=1, comment='') # 75
        trim2a.add(sid+1, mach, q, labels, uxs, aeqr=0.0,
                   trim_type=2, comment='') # 76
        trim1a.parse_cards()
        trim2a.parse_cards()

        labels = ['URDD3', 'URDD5', 'PITCH']
        uxs = [2.5, 0.0, 0.0]
        # good
        trim1b = TRIM(model)
        trim2b = TRIM(model)
        trim1b.add(sid+2, mach, q, labels, uxs, aeqr=0.0,
                   trim_type=1, comment='trim') # 77
        trim2b.add(sid+3, mach, q, labels, uxs, aeqr=0.0,
                   trim_type=2, comment='trim') # 78
        trim1b.parse_cards()
        trim2b.parse_cards()

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
        model.add_grid(55, [0., 0., 0.])
        model.add_grid(66, [0., 0., 0.])
        #model.add_cord2r(cid, origin, zaxis, xzplane, rid=0, comment='')
        model.setup()
        model.validate()

        # why doesn't this work?
        #if 'SUPORT1' in subcase:
            #suport_id = subcase.get_int_parameter('SUPORT1')[0]
            #suport1 = fem.suport1[suport_id]

        with self.assertRaises(RuntimeError):
            trim1a.verify_trim(suport1_id=0)
            #trim1a.verify_trim(model.suport, model.suport1, model.aestats, model.aeparams,
                               #model.aelinks, model.aesurf, xref=True)
        with self.assertRaises(RuntimeError):
            trim2a.verify_trim(suport1_id=0)
            #trim2a.verify_trim(model.suport, model.suport1, model.aestats, model.aeparams,
                               #model.aelinks, model.aesurf, xref=True)

        trim1b.verify_trim(suport1_id=0)
        trim2b.verify_trim(suport1_id=0)
        #trim1b.verify_trim(model.suport, model.suport1, model.aestats, model.aeparams,
                           #model.aelinks, model.aesurf, xref=True)
        #trim2b.verify_trim(model.suport, model.suport1, model.aestats, model.aeparams,
                           #model.aelinks, model.aesurf, xref=True)
        model.write_bdf('trim.bdf')
        model2 = read_bdf('trim.bdf', debug=None)
        model2._verify_bdf(xref=True)
        #model2.uncross_reference()
        model2._verify_bdf(xref=False)
        model2.cross_reference()
        model2._verify_bdf(xref=True)
        os.remove('trim.bdf')

        #model2.uncross_reference()
        #model2.safe_cross_reference()
        save_load_deck(model)

    def test_gust(self):
        """checks the GUST card"""
        sid = 100
        dload = 200
        wg = 50.
        x0 = 3.
        V = 42.
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        gust = model.add_gust(sid, dload, wg, x0, V=V, comment='gust load')
        model.setup(run_geom_check=True)
        #model.gust.validate()
        model.gust.write()

        gust2 = model.gust.add_card(BDFCard(['GUST', sid, dload, wg, x0, V]), comment='gust load')
        model.setup(run_geom_check=True)
        #model.gust.validate()
        model.gust.write()
        model.setup()
        save_load_deck(model)


    def test_csschd(self):
        """checks the CSSCHD card"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        csschd = model.csschd
        sid = 5
        aesid = 50
        lalpha = 12
        lmach = 15
        lschd = 25

        card = ['CSSCHD', sid, aesid, lalpha, lmach, lschd]
        bdf_card = BDFCard(card, has_none=True)
        csschd.add(sid, aesid, lschd, lalpha='lalpha', lmach=4,
                   comment='')
        with self.assertRaises(ValueError):
            csschd.parse_cards()

        csschd.clear()
        csschd.add(sid, aesid, lschd, lalpha=4, lmach=5.0,
                   comment='')
        csschd.parse_cards()

        #csschd.lalpha[0] = 4
        #csschd.lmach[0] = 5.0
        #with self.assertRaises(TypeError):
            #csschd.validate()
        csschd.lmach[0] = 5
        csschd.validate()

        card = ['CSSCHD', sid, aesid, lalpha, lmach, lschd]
        bdf_card = BDFCard(card, has_none=True)
        csschd.clear()
        csschd.add_card(bdf_card, comment='csschd card')
        csschd.validate()
        csschd.write()

        sid = 6
        csshcdi = model.add_csschd(sid, aesid, lschd, lalpha=lalpha, lmach=lmach,
                                   comment='csschd card')

        label = 'ELEV'
        cid1 = 0
        aelist_id1 = 37
        unused_aesurf = model.add_aesurf(
            aesid, label, cid1, aelist_id1, cid2=None, aelist_id2=None,
            eff=1.0, ldw='LDW', crefc=1.0, crefs=1.0,
            pllim=-np.pi/2., pulim=np.pi/2.,
            hmllim=None, hmulim=None,
            tqllim=None, tqulim=None, comment='aesurf')

        unused_aelist = model.add_aelist(aelist_id1, [1, 2, 3], comment='')

        aefact_sid = aelist_id1
        fractions = [0., 0.5, 1.]
        unused_aefact_elev = model.add_aefact(aefact_sid, fractions, comment='aefact')

        aefact_sid = lalpha
        fractions = [0., 5., 10.]
        unused_aefact_alpha = model.add_aefact(aefact_sid, fractions, comment='aefact')

        aefact_sid = lmach
        fractions = [0., 0.7, 0.8]
        unused_aefact_mach = model.add_aefact(aefact_sid, fractions, comment='aefact')

        aefact_sid = lschd
        fractions = [0., 15., 30., 45.]
        unused_aefact_delta = model.add_aefact(aefact_sid, fractions, comment='aefact')

        #model.cross_reference()
        model.setup()
        csschd.write()
        #csshcd1.write()
        #model.uncross_reference()

        bdf_filename = StringIO()
        model.write_bdf(bdf_filename, close=False)
        #model.safe_cross_reference()

        model.validate()
        save_load_deck(model)

        bdf_filename.seek(0)
        model2 = read_bdf(bdf_filename, punch=True, debug=False)

        bdf_filename2 = StringIO()
        model.write_bdf(bdf_filename2, size=16, close=False)
        model2.write_bdf(bdf_filename2, size=16, close=False)

        #-----------
        if 0:
            csshcd3 = CSSCHD(sid, aesid, lschd, lalpha=None, lmach=None, comment='cssch card')
            csshcd3.write_card()
            with self.assertRaises(RuntimeError):
                csshcd3.validate()
        save_load_deck(model)

    def test_bah_plane_bdf(self):
        """tests the bah_plane"""
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.bdf')
        folder = ''
        run_equivalence = False # buggy
        run_extract_bodies = False # was True in OG
        run_bdf(folder, bdf_filename, debug=False, xref=True, check=True,
                punch=False, mesh_form='combined',
                is_folder=False, print_stats=False,
                encoding=None, sum_load=True, size=8,
                is_double=False, stop=False, nastran='',
                post=-1, dynamic_vars=None, quiet=True,
                dumplines=False, dictsort=False,
                run_extract_bodies=run_extract_bodies,
                run_equivalence=run_equivalence,
                nerrors=0, dev=True,
                crash_cards=None, pickle_obj=True)

    def _test_rotord(self):
        """tests the ROTORD"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log)

        sid = 42
        rstart = 3.14
        rstep = .314
        numstep = 10
        rids = [None]
        rsets = [-31]
        rcords = [10]
        w3s = [13.]
        w4s = [3.]
        rforces = [14]
        brgsets = [17, False]
        rspeeds = [42.1]
        rotord = model.add_rotord(
            sid, rstart, rstep, numstep,
            rids, rsets, rspeeds, rcords, w3s, w4s, rforces, brgsets,
            refsys='ROT', cmout=0.0, runit='RPM', funit='RPM',
            zstein='NO', orbeps=1.e-6, roprt=0, sync=1, etype=1,
            eorder=1.0, threshold=0.02, maxiter=10, comment='rotord')
        rotord.validate()

        sid = 43
        nids = [100, 101, 102]
        rotorg = model.add_rotorg(
            sid, nids, comment='rotorg'
        )
        rotorg.validate()
        save_load_deck(model)

    def _test_build_structure_from_caero(self):
        model = BDF()
        eid = 1
        pid = 1
        igroup = 1
        p1 = [0., 0., 0.]
        p4 = [0., 10., 0.]
        x12 = 1.0
        x43 = 0.5
        model.add_caero1(eid, pid, igroup, p1, x12, p4, x43, cp=0, nspan=10, lspan=0, nchord=10, lchord=0, comment='')
        model.add_paero1(pid)
        caero_bdf_filename = os.path.join(TEST_PATH, 'caero.bdf')
        build_structure_from_caero(model, caero_bdf_filename)

        trim_load_cases = [
                # subcase, name, trim_id, mach, q, label: value
                (10, 'alpha=1', 10, 0.8, 1., {'ANGLEA': 1.,}),
        ]
        build_trim_load_cases(model, trim_load_cases, aeqr=1.0)


    def _test_zona_1(self):
        """zona explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8', log_func=None)  # lots of zona errors
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'f16_ma41.bdf')
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        save_load_deck(model, xref='safe',
                       run_renumber=False, run_convert=False, run_remove_unused=False,
                       run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
                       run_test_bdf=False, run_op2_writer=False)
        with self.assertRaises(NotImplementedError):
            model.zona.convert_to_nastran()

    def _test_zona_2(self):
        """zona explicit test"""
        log = SimpleLogger(level='error', encoding='utf-8', log_func=None)  # lots of zona errors
        bdf_filename = os.path.join(MODEL_PATH, 'aero', 'ztran.bdf')
        model = read_bdf(bdf_filename, xref=False, debug=None, log=log)
        model.safe_cross_reference()
        save_load_deck(model, xref='safe',
                       run_renumber=False, run_convert=False, run_remove_unused=False,
                       run_save_load=False, run_save_load_hdf5=False, run_mass_properties=False,
                       run_test_bdf=False)
        model.zona.convert_to_nastran()

    def _test_zona_3(self):
        """totally fake zona model"""
        bdf_file = get_zona_model()

        model = read_bdf(bdf_filename=bdf_file, validate=True, xref=True, punch=False,
                         skip_cards=None, read_cards=None, encoding=None,
                         log=None, debug=None, mode='zona')
        #with self.assertRaises(AttributeError):

        model.uncross_reference()
        model.write_bdf('zona.bdf')
        model.safe_cross_reference()
        model.write_bdf('zona.bdf')

        bdf_file.seek(0)
        model.clear_attributes()
        model2 = read_bdf('zona.bdf', debug=None)
        os.remove('zona.bdf')

        model2.zona.convert_to_nastran()

def get_zona_model():
    bdf_file = StringIO()
    bdf_file.write(
        '$ pyNastran: version=zona\n'
        'CEND\n'
        'BEGIN BULK\n'
        #'$       acsid, rcsid, cref, bref, sref, symxz, symxy\n'
        #'AEROZ, 10,     0,     1.,   10.,  100., YES\n'
        '$AEROZ  ACSID XZSYM FLIP FMMUNIT FMLUNIT REFC   REFB   REFS\n'
        '$       REFX  REFY  REFZ\n'
        'AEROZ,  0,    YES,  NO,  SLIN,   IN,      22.73,59.394,1175.8\n'
        ',       59.53,0.0,  0.0\n'

        '$       label, type, cid, PANLST, setg, actid\n'
        'AESURFZ,FLAP,  ASYM, 1,   10,       20,   0\n'
        #'AESURFZ,FLAP,  SYM,  1,  10,       20,   0\n'
        'CORD2R, 1,0, 0.,0.,0., 0.,0.,1.,\n'
        ',1.,0.,0.\n'
        '$BODY7,ID,LABEL,IPBODY7, ACOORD, NSEG, IDMESH1\n'
        'BODY7, 1, FUSE,        ,      2,     , 1\n'
        'PANLST3,10, FUSE, \n'
        '$       id,naxial,nradial, \n'
        'SEGMESH,1, 4,     3,       \n'

        # ITYPEi = 1 (Body of Revolution):
        #    Xi, CAMi, YRi
        # ITYPEi = 2 (Elliptical Body):
        #    Xi, YRi, ZRi
        # ITYPEi = 3 (Arbitrary Body):
        #    Xi, IDYi, IDZi
        '$       itype, x1, cam, yr1, zr1, idy1, idz1 \n'
        ',       1,        ,  1.,  1.,    ,\n'
        ',       2,      1.,    ,  1.,  2.,\n'
        ',       3,      2.,    ,    ,    , 13,   14   \n'
        ',       3,      3.,    ,    ,    , 13,   14   \n'

        # y
        'AEFACT,13, 1., 0.,  0.,-1.\n'
        'AEFACT,14, 0., 1., -1., 0.\n'
        '$ MKAEROZ, ID, MACH, METHOD, IDFLT\n'
        'MKAEROZ,   101, 0.8, -1,     -1,  \n'
        '$ TRIM, ID, MKAEROZ, Q,   LABEL1, UX1,    CGX, CGY,\n'
        'TRIM, 100,  101,     42., ALPHA,  5., 0., 0.,  0.,\n'
        '$CGZ, WEIGHT, Ixx, Ixy, Iyy, Ixz, Iyz, Izz\n'
        ',0.,  1e4,    1e3, 1e3, 1e5, 1e3, 1e3, 1e4\n'
        '$TRUE/G, NX,     NY,  NZ,  P,       Q,   R, \n'
        ', TRUE,  FREE, NONE,  32., FREE, NONE, 42., \n'
        '$var, value\n'
        ',17,  1.0,\n'
        '$\n'
        'TRIMVAR,17,VAR\n'
        '$\n'
        '$trimlnk,id,sym, ?,  ?\n'
        'TRIMLNK,10,SYM, -1, 17\n'
        'ACOORD, 2, 0.,0.,0., 1.0,0.\n'
        '$       ID,    MODEL, CP, PANLST, SETG, DZ, EPS\n'
        'SPLINE1,100,        ,   ,    422, 423,\n'
        '$,      NELEM, MELEM\n'
        '$,      10,    10\n'
        'PANLST3,422, FUSE, \n'
        '$       id,naxial,nradial, \n'

        #'$       ID,   MODEL, PANLST, SETG,\n'
        #'SPLINE2,1000,      ,    422,  423,\n'

        '$       ID,   MODEL, CP, PANLST, SETG,\n'
        'SPLINE3,1200,      ,   ,    422,  423,\n'

        'SET1,423,10\n'
        'GRID,10,,0.,0.,0.\n'
        'GRID,11,,1.,0.,0.\n'
        'CONROD,100, 10,11, 101,1.0\n'
        'MAT1,101,3.0e7,,0.3\n'
    )
    bdf_file.seek(0)
    return bdf_file


def build_structure_from_caero(model: BDF,
                               caero_bdf_filename: str,
                               write_panel_xyz: bool=True):
    export_caero_mesh(model, caero_bdf_filename=caero_bdf_filename,
                      is_subpanel_model=True, pid_method='caero',
                      write_panel_xyz=write_panel_xyz)
    model_quads = read_bdf(caero_bdf_filename)
    model.nodes = model_quads.nodes
    model.elements = model_quads.elements
    model.properties = model_quads.properties
    model.materials = model_quads.materials
    os.remove(caero_bdf_filename)
    return model

def _setup_aero_plot(fig_id: Optional[int]=None) -> tuple[Any, Any]:
    """helper for plotting aero panels"""
    fig = None
    ax = None
    if IS_MATPLOTLIB:
        fig = plt.figure(fig_id)
        ax = fig.gca()
        ax.set_ylabel('Y')
        ax.set_xlabel('X')
        ax.grid()
    return fig, ax

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
