from collections import Counter
import unittest

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.field_writer_8 import print_int_card_blocks
from pyNastran.bdf.cards.bdf_sets import (
    SET1, SET2, SET3, ASET, ASET1, OMIT1, BSET, BSET1, CSET, CSET1, QSET, QSET1, USET, USET1,
    SEBSET, SEBSET1, SECSET, SECSET1, SEQSET, SEQSET1, #SEUSET, SEUSET1,
)
from pyNastran.bdf.cards.test.utils import save_load_deck


class TestSets(unittest.TestCase):

    def test_set1_01(self):
        bdf = BDF(debug=False)
        lines = ['SET1,    1100,    100,     101']
        card = bdf._process_card(lines)
        card = BDFCard(card)

        size = 8
        card = SET1.add_card(card)
        card.write_card(size, 'dummy')
        card.raw_fields()

        card2 = SET1(1100, [100, 101], is_skin=False, comment='')
        card2.write_card(size, 'dummy')

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
        sid = 10
        ids = [1, 2, 3, 4, 5]
        set1a = SET1(sid, ids, is_skin=False, comment='set1')
        set1b = SET1.add_card(BDFCard(['SET1', sid] + ids))
        set1a.write_card()

    def test_set2_01(self):
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

    def test_set2_02(self):
        """checks the SET2 card"""
        bdf = BDF(debug=True)

        set2 = bdf.add_set2(110, 10, -0.1, 1.1, -0.1, 1.1)
        caero = bdf.add_caero4(10, 10, [.0, .0, .0], 1., [.0, 1., .0], 1.)
        spline2 = bdf.add_spline2(10, 10, 10, 11, 110)

        spline2.cross_reference(bdf)

        self.assertEqual(spline2.setg_ref, set2)
        self.assertEqual(set2.macro_ref, caero)

    def test_set3_01(self):
        """checks the SET3 card"""
        model = BDF(debug=False)
        sid = 10
        ids = [1, 2, 3, 4, 5]
        desc = 'ELEM'
        set3a = SET3(sid, desc, ids, comment='set3')
        model.sets[sid] = set3a
        model.add_card(BDFCard(['SET3', sid+1, desc] + ids), 'SET3', comment='set3')
        set3b = model.sets[sid]
        set3a.write_card()
        set3a.validate()
        set3b.validate()
        save_load_deck(model)

    def test_set3_02(self):
        """checks the SET3 card"""
        model = BDF()

        # List of grid IDs
        grid_list = [1, 2, 3, 4, 5, 6, 7, 13, 15,
                     20, 21, 22, 23, 30, 31, 32, 33]

        # Define the card lines
        card_lines = ['SET3', 1, 'GRID'] + grid_list

        # Add nastran card to BDF object
        model.add_card(card_lines, 'SET3', comment='set3-1', is_list=True)
        fields = model.sets[1].raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        str(model.sets[1])

        set3a = SET3(2, 'GRID', grid_list, comment='set3-2')
        fields = model.sets[1].raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        str(set3a)

    def test_aset(self):
        """checks the ASET/ASET1 cards"""
        model = BDF(debug=False)
        aset1a = ASET1([1, 'THRU', 10], 4, comment='aset')
        aset1b = ASET1.add_card(BDFCard(['ASET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='aset1')
        aset1a.write_card()
        aset1b.write_card()
        model._add_aset_object(aset1a)
        model._add_aset_object(aset1b)
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        aseta = ASET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        asetb = ASET.add_card(BDFCard(['ASET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]))

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        aseta.validate()
        asetb.validate()
        aseta.write_card()
        asetb.write_card()
        save_load_deck(model)

    def test_omit(self):
        """checks the OMIT/OMIT1 cards"""
        model = BDF(debug=False)
        omit1a = OMIT1([1, 'THRU', 10], 4, comment='omit1')
        self.assertEqual(omit1a.components, 4)
        omit1b = OMIT1.add_card(BDFCard(['OMIT1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='omit1')
        model._add_omit_object(omit1a)
        model._add_omit_object(omit1b)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])

        omit1a.validate()
        omit1b.validate()
        omit1a.write_card()
        omit1b.write_card()
        #| OMIT1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        #omita = OMIT([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        #omitb = OMIT.add_card(BDFCard(['OMIT',
                                       #1, 2, 3, 4, 5,
                                       #5, 4, 3, 2, 1]))
        #omita.validate()
        #omitb.validate()
        #omita.write_card()
        #omitb.write_card()
        save_load_deck(model)

    def test_bset(self):
        """checks the BSET/BSET1 cards"""
        model = BDF(debug=False)
        bset1a = BSET1([1, 'THRU', 10], 4, comment='bset1')
        bset1b = BSET1.add_card(BDFCard(['BSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='bset1')
        bset1a.write_card()
        bset1b.write_card()
        model._add_bset_object(bset1a)
        model._add_bset_object(bset1b)
        #| BSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        bseta = BSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='bset')
        bsetb = BSET.add_card(BDFCard(['BSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]), comment='bset')
        model._add_bset_object(bseta)
        model._add_bset_object(bsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        bseta.validate()
        bsetb.validate()
        bseta.write_card()
        bsetb.write_card()
        save_load_deck(model)

    def test_cset(self):
        """checks the CSET/CSET1 cards"""
        model = BDF(debug=False)
        cset1a = CSET1([1, 'THRU', 10], 4, comment='cset')
        cset1b = CSET1.add_card(BDFCard(['CSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='cset1')
        cset1a.write_card()
        cset1b.write_card()
        model._add_cset_object(cset1a)
        model._add_cset_object(cset1b)
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        cseta = CSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='cset')
        csetb = CSET.add_card(BDFCard(['CSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]), comment='cset')
        model._add_cset_object(cseta)
        model._add_cset_object(csetb)
        model.add_cset([1, 2, 3], '42', comment='cset')
        model.add_cset1([1, 2, 3], [1, 2, 3], comment='cset1')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        cseta.validate()
        csetb.validate()
        cseta.write_card()
        csetb.write_card()
        save_load_deck(model)

    def test_qset(self):
        """checks the QSET/QSET1 cards"""
        model = BDF(debug=False)
        qset1a = QSET1([1, 'THRU', 10], 4, comment='qset')
        qset1b = QSET1.add_card(BDFCard(['QSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='qset1')
        model._add_qset_object(qset1a)
        model._add_qset_object(qset1b)
        qset1a.write_card()
        qset1b.write_card()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        qseta = QSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='qset')
        qsetb = QSET.add_card(BDFCard(['QSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]), comment='qset')
        model._add_qset_object(qseta)
        model._add_qset_object(qsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        qseta.validate()
        qsetb.validate()
        qseta.write_card()
        qsetb.write_card()
        save_load_deck(model)

    def test_uset(self):
        """checks the USET/USET1 cards"""
        model = BDF(debug=False)
        uset1a = USET1('MYSET1', [1, 'THRU', 10], 4, comment='uset')
        fields = ['USET1', 'MYSET2',
                  5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
        model.add_card(fields, 'USET1', comment='uset1')
        model._add_uset_object(uset1a)
        #model._add_uset_object(uset1b)
        uset1a.write_card()
        #uset1b.write_card()

        useta = USET('MYSET3', [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='uset')
        fields = ['USET', 'MYSET4',
                  1, 2, 3, 4, 5,
                  5, 4, 3, 2, 1]
        model.add_card(fields, 'USET', comment='uset')
        model._add_uset_object(useta)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        useta.validate()
        model.validate()
        useta.write_card()
        save_load_deck(model)

    def test_sebset(self):
        """checks the SEBSET/SEBSET1 cards"""
        model = BDF(debug=False)
        seid = 42
        bset1a = SEBSET1(seid, [1, 'THRU', 10], 4, comment='bset1')
        bset1b = SEBSET1.add_card(BDFCard(['SEBSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                comment='sebset1')
        bset1a.write_card()
        bset1b.write_card()
        model._add_sebset_object(bset1a)
        model._add_sebset_object(bset1b)
        #| BSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        sebseta = SEBSET(seid, [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='sebset')
        sebsetb = SEBSET.add_card(BDFCard(['SEBSET', seid,
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]), comment='sebset')
        assert len(sebseta.components) == 5, sebseta.components
        model._add_sebset_object(sebseta)
        model._add_sebset_object(sebsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        sebseta.validate()
        sebsetb.validate()
        sebseta.write_card()
        sebsetb.write_card()
        save_load_deck(model)

    def test_secset(self):
        """checks the SECSET/SECSET1 cards"""
        model = BDF(debug=False)
        seid = 171
        secset1a = SECSET1(seid, [1, 'THRU', 10], 4, comment='cset')
        secset1b = SECSET1.add_card(BDFCard(['SECSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]),
                                    comment='secset1')
        secset1a.write_card()
        secset1b.write_card()
        model._add_secset_object(secset1a)
        model._add_secset_object(secset1b)
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        secseta = SECSET(seid, [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='secset')
        secsetb = SECSET.add_card(BDFCard(['SECSET', seid,
                                           1, 2, 3, 4, 5,
                                           5, 4, 3, 2, 1]), comment='secset')
        model._add_secset_object(secseta)
        model._add_secset_object(secsetb)
        model.add_secset(seid, [1, 2, 3], '42', comment='secset')
        model.add_secset1(seid, [1, 2, 3], [1, 2, 3], comment='secset1')

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        secseta.validate()
        secsetb.validate()
        secseta.write_card()
        secsetb.write_card()
        save_load_deck(model)

    def test_seqset(self):
        """checks the QSET/QSET1 cards"""
        model = BDF(debug=False)
        seid = 42
        seqset1a = SEQSET1(seid, [1, 'THRU', 10], 4, comment='qset')
        model.add_card(['SEQSET1', seid, 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9],
                       'SEQSET1', comment='seqset1')
        model._add_seqset_object(seqset1a)
        seqset1a.write_card()
        #| SEQSET1 | SEID |  C  | ID1 | THRU | ID2 |

        seqseta = SEQSET(seid, [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='seqset')
        fields = ['SEQSET', seid,
                  1, 2, 3, 4, 5,
                  5, 4, 3, 2, 1]
        model.add_card(fields, 'SEQSET', comment='seqset')
        model._add_seqset_object(seqseta)
        #model._add_seqset_object(seqsetb)

        nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for nid in nids:
            model.add_grid(nid, [float(nid), 0., 0.])
        seqseta.validate()
        seqseta.write_card()
        save_load_deck(model)

    #def test_seuset(self):
        #"""checks the SEUSET/SEUSET1 cards"""
        #model = BDF(debug=False)
        #seuset1a = SEUSET1('MYSET1', [1, 'THRU', 10], 4, comment='seuset')
        #fields = ['SEUSET1', 'MYSET2',
                  #5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
        #model.add_card(fields, 'SEUSET1', comment='seuset1')
        #model._add_seuset_object(seuset1a)
        ##model._add_uset_object(uset1b)
        #seuset1a.write_card()
        ##seuset1b.write_card()

        #useta = SEUSET('MYSET3', seid, [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], comment='seuset')
        #fields = ['SEUSET', seid, 'MYSET4',
                  #1, 2, 3, 4, 5,
                  #5, 4, 3, 2, 1]
        #model.add_card(fields, 'SEUSET', comment='seuset')
        #model._add_uset_object(useta)

        #nids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        #for nid in nids:
            #model.add_grid(nid, [float(nid), 0., 0.])
        #seuseta.validate()
        #model.validate()
        #seuseta.write_card()
        #save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
