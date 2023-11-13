from collections import Counter
import numpy as np
import unittest

from cpylog import SimpleLogger
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.field_writer_8 import print_int_card_blocks
#from pyNastran.bdf.cards.bdf_sets import (
    #SET2, SET3,
    #SEBSET, SEBSET1, SECSET, SECSET1, SEQSET, SEQSET1, #SEUSET, SEUSET1,
#)


class TestSets(unittest.TestCase):

    def test_monpnt(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log)
        monpnt1 = model.monpnt1
        monpnt3 = model.monpnt3

        name = 'test'
        label = 'test2'
        axes = '123'
        comp = 'WING'
        xyz = [0., 0., 0.]
        model.add_monpnt1(name, label, axes, comp, xyz, cp=0,
                          cd=None, comment='monpnt1')
        #monpnt1.raw_fields()
        monpnt1.validate()

        Type = 'CQUAD4'
        table = 'STRESS'
        nddl_item = 'SX1'
        eid = 17
        model.add_monpnt2(name, label, table, Type, nddl_item, eid,
                          comment='monpnt2')
        monpnt2 = model.monpnt2
        #monpnt2.raw_fields()
        #monpnt2.validate()

        grid_set = 43
        elem_set = 44
        model.add_monpnt3(name, label, axes, grid_set, elem_set,
                          xyz, cp=0, cd=None,
                          xflag=None, comment='monpnt3')
        #monpnt3.raw_fields()
        monpnt3.validate()

        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)
        #model.uncross_reference()

        save_load_deck(model)

    def test_set1_01(self):
        model = BDF(debug=False)
        set1 = model.set1

        lines = ['SET1,    1100,    100,     101']
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        card = set1.add_card(card)
        model.setup()
        set1.write(size, 'dummy')
        #card.raw_fields()

        unused_card2 = model.add_set1(1100, [100, 101], is_skin=False, comment='')
        model.setup()
        set1.write(size, 'dummy')
        save_load_deck(model)

    def test_set1_02(self):
        """checks the SET1 card"""
        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3], True], # these are all integers
        ]
        msg = print_int_card_blocks(fields_blocks)
        self.assertEqual('SET1           a      1.       3       1       2       3\n', msg)

        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3, 5, 4], True], # these are all integers
        ]
        msg2 = print_int_card_blocks(fields_blocks)
        #print('%r' % msg2)
        self.assertEqual('SET1           a      1.       3       1       2       3       5       4\n', msg2)

        fields_blocks = [
            'SET1',
            [['a', 1.0, 3], False], # these are not all integers
            [[1, 2, 3, 5, 4, 6], True], # these are all integers
        ]
        msg3 = print_int_card_blocks(fields_blocks)
        #print('%r' % msg3)
        self.assertEqual('SET1           a      1.       3       1       2       3       5       4\n'
                         '               6\n', msg3)

    def test_set1_03(self):
        """checks the SET1 card"""
        model = BDF(debug=False)
        set1 = model.set1
        sid = 10
        ids = [1, 2, 3, 4, 5]
        unused_set1a = model.add_set1(sid, ids, is_skin=False, comment='set1')
        unused_set1b = set1.add_card(BDFCard(['SET1', sid] + ids))
        model.setup()
        #set1a.write_card()
        #set1b.write_card()
        str(set1.write())

    def _test_set2_01(self):
        """checks the SET2 card"""
        bdf = BDF(debug=False)
        lines = ['SET2,     110,      10,   -0.1,    1.1,   -0.1,     1.1']
        card = bdf._process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SET2.add_card(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

        card2 = SET2(110, 10, -0.1, 1.1, -0.1, 1.1, comment='')
        card2.write_card(size, 'dummy')
        card2.raw_fields()

        card3 = bdf.add_set2(110, 10, -0.1, 1.1, -0.1, 1.1)
        card3.write_card(size, 'dummy')
        card3.raw_fields()

    def _test_set2_02(self):
        """checks the SET2 card"""
        bdf = BDF(debug=True)

        set2 = bdf.add_set2(110, 10, -0.1, 1.1, -0.1, 1.1)
        caero = bdf.add_caero4(10, 10, [.0, .0, .0], 1., [.0, 1., .0], 1.)
        spline2 = bdf.add_spline2(10, 10, 10, 11, 110)
        bdf.pop_parse_errors()

        spline2.cross_reference(bdf)

        self.assertEqual(spline2.setg_ref, set2)
        self.assertEqual(set2.macro_ref, caero)

    def test_set3_01(self):
        """checks the SET3 card"""
        model = BDF(debug=False)
        set3 = model.set3
        sid = 10
        ids = [1, 2, 3, 4, 5]
        desc = 'ELEM'
        unused_set3a = model.add_set3(sid, desc, ids, comment='set3')
        #model.set3[sid] = set3a
        set3.add_card(BDFCard(['SET3', sid+1, desc] + ids), comment='set3')
        model.setup(run_geom_check=True)
        unused_set3b = set3.slice_card_by_id(sid+1)
        set3.write(size=8)
        set3.write(size=16)
        #set3a.validate()
        #set3b.validate()
        save_load_deck(model)

    def test_set3_02(self):
        """checks the SET3 card"""
        model = BDF()
        set3 = model.set3

        # list of grid IDs
        grid_list = [1, 2, 3, 4, 5, 6, 7, 13, 15,
                     20, 21, 22, 23, 30, 31, 32, 33]

        # Define the card lines
        card_lines = ['SET3', 1, 'GRID'] + grid_list

        # Add nastran card to BDF object
        model.add_card(card_lines, 'SET3', comment='set3-1', is_list=True)
        model.setup(run_geom_check=True)
        set3b = model.set3.slice_card_by_id(1)
        fields = set3b.write().split()
        # .raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        #str(model.sets[1])
        set3.write(size=8)
        set3.write(size=16)

        set3a = model.add_set3(2, 'GRID', grid_list, comment='set3-2')
        model.setup(run_geom_check=True)
        fields = model.set3.slice_card_by_id(1).write().split()
        # .raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        str(set3a)

    def test_set3_03(self):
        model = BDF()
        card_lines = [
            'SET3, 2, ELEM, 20, THRU, 33, 36, THRU, 44',
            '49, THRU, 62, 91, THRU, 110',
        ]
        model.add_card(card_lines, 'SET3', is_list=False)
        model.setup()
        print(model.set3.write())
        save_load_deck(model)

    def test_radset_02(self):
        """checks the RADSET card"""
        model = BDF()
        radset = model.radset

        # list of grid IDs
        cavity_ids = [1, 2, 3, 4, 5, 6, 7, 13, 15,
                     20, 21, 22, 23, 30, 31, 32, 33]

        # Define the card lines
        card_lines = ['RADSET'] + cavity_ids

        # Add nastran card to BDF object
        model.add_card(card_lines, 'RADSET', comment='radset-1', is_list=True)
        model.setup(run_geom_check=True)
        unused_radsetb = radset.slice_card_by_id(1)
        fields = radset.write().split()
        # .raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        #str(model.sets[1])
        radset.write(size=8)
        radset.write(size=16)

        unused_radseta = model.add_radset(cavity_ids, comment='radset-2')
        model.setup(run_geom_check=True)
        fields = radset.slice_card_by_id(1).write().split()
        # .raw_fields()
        #thru_count = Counter(fields)['THRU']
        #assert thru_count in [0, 1], fields
        str(radset)

    def test_radset_02(self):
        """checks the RADSET card"""
        model = BDF()
        seset = model.seset

        # list of grid IDs
        node_ids = [1, 2, 3, 4, 5, 6, 7, 13, 15,
                    20, 21, 22, 23, 30, 31, 32, 33]

        # Define the card lines
        card_lines = ['SESET', 1] + node_ids

        # Add nastran card to BDF object
        model.add_card(card_lines, 'SESET', comment='seset-1', is_list=True)
        model.setup(run_geom_check=True)
        unused_sesetb = seset.slice_card_by_id(1)
        fields = seset.write().split()
        # .raw_fields()
        thru_count = Counter(fields)['THRU']
        #assert thru_count in [0, 1], fields
        assert thru_count == 3, seset.write()
        #str(model.sets[1])
        seset.write(size=8)
        seset.write(size=16)

        seid = 2
        unused_seseta = model.add_seset(seid, node_ids, comment='seset')
        model.setup(run_geom_check=True)
        fields = seset.slice_card_by_id(1).write().split()
        # .raw_fields()
        #thru_count = Counter(fields)['THRU']
        #assert thru_count in [0, 1], fields
        str(seset)

    def test_aset(self):
        """checks the ASET/ASET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        aset = model.aset
        unused_aset1a = aset.add_set1([1, 'THRU', 10], 4, comment='aset')
        card = BDFCard(['ASET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
        unused_aset1b = aset.add_set1_card(card, comment='aset1')
        aset.write()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_aseta = aset.add_set([1, 2, 3, 4, 5],
                                    [5, 4, 3, 2, 1])
        card = BDFCard(['ASET',
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_asetb = aset.add_set_card(card)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        aset.validate()
        aset.write()
        aset.set_map
        save_load_deck(model)

    def test_omit(self):
        """checks the OMIT/OMIT1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        omit = model.omit
        unused_omit1a = omit.add_set1([1, 'THRU', 10], 4, comment='omit1')
        model.setup()
        assert np.array_equal(omit.component, np.ones(10)*4)
        card = BDFCard(['OMIT1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
        unused_omit1b = omit.add_set1_card(card, comment='omit1')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])

        omit.validate()
        omit.write()
        omit.set_map
        #| OMIT1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        #omita = OMIT([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        #omitb = OMIT.add_card(BDFCard(['OMIT',
                                       #1, 2, 3, 4, 5,
                                       #5, 4, 3, 2, 1]))
        #omita.validate()
        #omitb.validate()
        #omitb.write()
        save_load_deck(model)

    def test_bset(self):
        """checks the BSET/BSET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        bset = model.bset
        unused_bset1a = bset.add_set1([1, 'THRU', 10], 4, comment='bset1')
        card = BDFCard(['BSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
        unused_bset1b = bset.add_set1_card(card, comment='bset1')
        bset.write()
        #| BSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_bseta = bset.add_set([1, 2, 3, 4, 5],
                                    [5, 4, 3, 2, 1], comment='bset')
        card = BDFCard(['BSET',
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_bsetb = bset.add_set_card(card, comment='bset')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        bset.validate()
        bset.write()
        bset.set_map
        save_load_deck(model)

    def test_cset(self):
        """checks the CSET/CSET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        cset = model.cset
        unused_cset1a = cset.add_set1([1, 'THRU', 10], 4, comment='cset')
        card = BDFCard(['CSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])

        unused_cset1b = cset.add_set1_card(card, comment='cset1')
        cset.write()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_cseta = cset.add_set([1, 2, 3, 4, 5],
                             [5, 4, 3, 2, 1], comment='cset')
        card = BDFCard(['CSET',
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_csetb = cset.add_set_card(card, comment='cset')
        model.add_cset([1, 2, 3], '42', comment='cset')
        model.add_cset1([1, 2, 3], [1, 2, 3], comment='cset1')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        cset.validate()
        cset.write()
        cset.set_map
        save_load_deck(model)

    def test_qset(self):
        """checks the QSET/QSET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        qset = model.qset
        unused_qset1a = qset.add_set1([1, 'THRU', 10], 4, comment='qset')
        card = BDFCard(['QSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])

        unused_qset1b = qset.add_set1_card(card, comment='qset1')
        qset.write()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_qseta = qset.add_set([1, 2, 3, 4, 5],
                                    [5, 4, 3, 2, 1], comment='qset')
        card = BDFCard(['QSET',
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_qsetb = qset.add_set_card(card, comment='qset')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        qset.validate()
        qset.write()
        qset.set_map
        save_load_deck(model)

    def test_qset_add_card(self):
        """checks the QSET/QSET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        qset = model.qset
        unused_qset1a = model.add_qset1([1, 'THRU', 10], 4, comment='qset')
        model.setup()

        qset.write()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_qseta = model.add_qset([1, 2, 3, 4, 5],
                                      [5, 4, 3, 2, 1], comment='qset')
        model.setup()

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])

        model.setup()
        qset.validate()
        qset.write()
        qset.set_map
        save_load_deck(model)

    def test_uset(self):
        """checks the USET/USET1 cards"""
        model = BDF(debug=False)
        #add_methods = model._add_methods
        uset = model.uset
        unused_uset1a = model.add_uset1('MYSET1', [1, 'THRU', 10], 4, comment='uset')
        model.setup()
        fields = ['USET1', 'MYSET2', 5,
                  1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
        model.add_card(fields, 'USET1', comment='uset1')
        model.setup()
        #add_methods._add_uset_object(uset1a)
        #add_methods._add_uset_object(uset1b)
        uset.write()

        unused_useta = model.add_uset(
            'MYSET3',
            [1, 2, 3, 4, 5],
            [5, 4, 3, 2, 1], comment='uset')
        model.setup()
        fields = ['USET', 'MYSET4',
                  1, 2, 3, 4, 5,
                  5, 4, 3, 2, 1]
        model.setup()
        model.add_card(fields, 'USET', comment='uset')
        #add_methods._add_uset_object(useta)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        uset.validate()
        model.validate()
        uset.write()
        save_load_deck(model)

    def test_sebset(self):
        """checks the SEBSET/SEBSET1 cards"""
        model = BDF(debug=False)
        sebset = model.sebset
        #add_methods = model._add_methods
        seid = 42
        unused_bset1a = model.add_sebset(seid, [1, 'THRU', 10], 4, comment='bset1')
        card = BDFCard(['SEBSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
        unused_bset1b = model.sebset.add_set1_card(card, comment='sebset1')
        sebset.write()
        #| BSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        seid = 50
        unused_sebseta_id = model.add_sebset(seid,
                                             [1, 2, 3, 4, 5],
                                             [5, 4, 3, 2, 1], comment='sebset')
        card = BDFCard(['SEBSET', seid + 1,
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_sebsetb_id = model.sebset.add_set1_card(card, comment='sebset')
        model.setup()
        sebseta = sebset.slice_card_by_id(seid)
        unused_sebsetb = sebset.slice_card_by_id(seid+1)
        assert len(sebseta.component) == 5, sebseta.component
        #add_methods._add_sebset_object(sebseta)
        #add_methods._add_sebset_object(sebsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        sebset.validate()
        sebset.write()
        save_load_deck(model)

    def test_secset(self):
        """checks the SECSET/SECSET1 cards"""
        model = BDF(debug=False)
        secset = model.secset
        seid = 171
        unused_secset1a = model.add_secset1(seid, [1, 'THRU', 10], 4, comment='cset')
        card = BDFCard(['SECSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9])
        unused_secset1b = secset.add_set1_card(card, comment='secset1')
        secset.write()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        unused_secseta = model.add_secset(seid,
                                          [1, 2, 3, 4, 5],
                                          [5, 4, 3, 2, 1], comment='secset')
        card = BDFCard(['SECSET', seid,
                        1, 2, 3, 4, 5,
                        5, 4, 3, 2, 1])
        unused_secsetb = secset.add_set_card(card, comment='secset')
        model.add_secset(seid, [1, 2, 3], '42', comment='secset')
        model.add_secset1(seid, [1, 2, 3], [1, 2, 3], comment='secset1')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        secset.validate()
        secset.write()
        save_load_deck(model)

    def test_seqset(self):
        """checks the QSET/QSET1 cards"""
        model = BDF(debug=False)
        seqset = model.seqset
        seid = 42
        unused_seqset1a = model.add_seqset1(seid, [1, 'THRU', 10], 4, comment='qset')
        model.add_card(['SEQSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9],
                       'SEQSET1', comment='seqset1')
        seqset.write()
        #| SEQSET1 | SEID |  C  | ID1 | THRU | ID2 |

        unused_seqseta = model.add_seqset(
            seid,
            [1, 2, 3, 4, 5],
            [5, 4, 3, 2, 1], comment='seqset')
        fields = ['SEQSET', seid,
                  1, 2, 3, 4, 5,
                  5, 4, 3, 2, 1]
        model.add_card(fields, 'SEQSET', comment='seqset')
        #add_methods._add_seqset_object(seqsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        seqset.validate()
        seqset.write()
        save_load_deck(model)

    #def test_seuset(self):
        #"""checks the SEUSET/SEUSET1 cards"""
        #model = BDF(debug=False)
        #seuset1a = SEUSET1('MYSET1', [1, 'THRU', 10], 4, comment='seuset')
        #fields = ['SEUSET1', 'MYSET2',
                  #5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
        #model.add_card(fields, 'SEUSET1', comment='seuset1')
        #add_methods._add_seuset_object(seuset1a)
        ##add_methods._add_uset_object(uset1b)
        #seuset1.write()

        #useta = SEUSET('MYSET3', seid, [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='seuset')
        #fields = ['SEUSET', seid, 'MYSET4',
                  #1, 2, 3, 4, 5,
                  #5, 4, 3, 2, 1]
        #model.add_card(fields, 'SEUSET', comment='seuset')
        #add_methods._add_uset_object(useta)

        #nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        #for nid in nids:
            #model.add_grid(nid, [float(nid), 0., 0.])
        #seuseta.validate()
        #model.validate()
        #seuseta.write_card()
        #save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
