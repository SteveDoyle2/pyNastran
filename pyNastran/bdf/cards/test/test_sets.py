from __future__ import print_function
from collections import Counter
import unittest

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.field_writer_8 import print_int_card_blocks
from pyNastran.bdf.cards.bdf_sets import (
    SET1, SET3, ASET, ASET1, BSET, BSET1, CSET, CSET1, QSET, QSET1,
)

class TestSets(unittest.TestCase):

    def test_set1_01(self):
        bdf = BDF(debug=False)
        lines = ['SET1,    1100,    100,     101']
        card = bdf.process_card(lines)
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
        set1a = SET1(sid, ids, is_skin=False, comment='')
        set1b = SET1.add_card(BDFCard(['SET1', sid] + ids))
        set1a.write_card()

    def test_set3_01(self):
        """checks the SET3 card"""
        sid = 10
        ids = [1, 2, 3, 4, 5]
        desc = 'ELEM'
        set3a = SET3(sid, desc, ids, comment='')
        set3b = SET3.add_card(BDFCard(['SET3', sid, desc] + ids))
        set3a.write_card()

    def test_set3_02(self):
        """checks the SET3 card"""
        model = BDF()

        # List of grid IDs
        grid_list = [1, 2, 3, 4, 5, 6, 7, 13, 15,
                     20, 21, 22, 23, 30, 31, 32, 33]

        # Define the card lines
        card_lines = ['SET3', 1, 'GRID'] + grid_list

        # Add nastran card to BDF object
        model.add_card(card_lines, 'SET3', comment='', is_list=True)
        fields = model.sets[1].raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        str(model.sets[1])

        set3a = SET3(2, 'GRID', grid_list, comment='')
        fields = model.sets[1].raw_fields()
        thru_count = Counter(fields)['THRU']
        assert thru_count in [0, 1], fields
        str(set3a)


    def test_aset(self):
        """checks the ASET/ASET1 cards"""
        aset1a = ASET1(4, [1, 'THRU', 10])
        aset1b = ASET1.add_card(BDFCard(['ASET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]))
        aset1a.write_card()
        aset1b.write_card()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        aseta = ASET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        asetb = ASET.add_card(BDFCard(['ASET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]))
        aseta.validate()
        asetb.validate()
        aseta.write_card()
        asetb.write_card()


    def test_bset(self):
        """checks the BSET/BSET1 cards"""
        bset1a = BSET1(4, [1, 'THRU', 10])
        bset1b = BSET1.add_card(BDFCard(['BSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]))
        bset1a.write_card()
        bset1b.write_card()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        bseta = BSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        bsetb = BSET.add_card(BDFCard(['BSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]))
        bseta.validate()
        bsetb.validate()
        bseta.write_card()
        bsetb.write_card()

    def _test_cset(self):
        """checks the CSET/CSET1 cards"""
        cset1a = CSET1(4, [1, 'THRU', 10])
        cset1b = CSET1.add_card(BDFCard(['CSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]))
        cset1a.write_card()
        cset1b.write_card()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        cseta = CSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        csetb = CSET.add_card(BDFCard(['CSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]))
        cseta.validate()
        csetb.validate()
        cseta.write_card()
        csetb.write_card()


    def test_qset(self):
        """checks the QSET/QSET1 cards"""
        qset1a = QSET1(4, [1, 'THRU', 10])
        qset1b = QSET1.add_card(BDFCard(['QSET1', 5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9]))
        qset1a.write_card()
        qset1b.write_card()
        #| ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |

        qseta = QSET([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])
        qsetb = QSET.add_card(BDFCard(['QSET',
                                       1, 2, 3, 4, 5,
                                       5, 4, 3, 2, 1]))
        qseta.validate()
        qsetb.validate()
        qseta.write_card()
        qsetb.write_card()


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
