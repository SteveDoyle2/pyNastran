"""
tests:
 - CELAS1, CELAS2, CELAS3, CELAS4
 - PELAS, PELAST

"""
import unittest

from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck


class TestSprings(unittest.TestCase):
    """
    tests:
     - CELAS1, CELAS2, CELAS3, CELAS4
     - PELAS, PELAST

    """
    def test_pelas_01(self):
        """tests PELAS"""
        lines = ['pelas, 201, 1.e+5']
        model = BDF(debug=False)
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        elem_id = model.pelas.add_card(card)
        pelas = model.pelas
        model.setup()
        pelas.write(size, 'dummy')
        #pelas.raw_fields()
        self.assertEqual(pelas.property_id, 201)
        self.assertEqual(pelas.k, 1e5)

    def test_pelas_02(self):
        """tests PELAS"""
        fields = ['pelas', 201, 1.e+5, None, None, 202, 2.e+5]
        model = BDF(debug=False)
        #model.echo = True
        card_name = fields[0]
        model.add_card(fields, card_name, comment='', is_list=True,
                       has_none=True)
        model.setup()
        assert len(model.pelas) == 2, model.pelas

    def test_springs_03(self):
        """tests CELAS1, CELAS2, PELAS, PELAST"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        nids = [3, 4]
        c1 = 1
        c2 = 1
        celas1_id = model.add_celas1(eid, pid, nids, c1, c2, comment='celas1')
        celas1 = model.celas1

        k = 1.0e7
        ge = 0.0
        s = 0.
        pelas_id = model.add_pelas(pid, k, ge, s, comment='pelas')

        model.setup()
        #celas1.raw_fields()
        #pelas.raw_fields()

        tkid = 10
        tgeid = 10
        tknid = 10
        pelast_id = model.add_pelast(pid, tkid, tgeid, tknid, comment='pealst')

        model.setup()
        #pelast.raw_fields()

        eid = 5
        celas2_id = model.add_celas2(eid, k, nids, c1=0, c2=0, ge=0., s=0., comment='celas2')
        #celas2.raw_fields()

        model.add_grid(3, [0., 0., 0.])
        model.add_grid(4, [0., 0., 0.])
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()

        pelas = model.pelas
        celas1 = model.celas1
        celas2 = model.celas2
        pelast = model.pelast

        celas1.write(size=8, is_double=False)
        celas2.write(size=8, is_double=False)
        pelas.write(size=8, is_double=False)
        pelast.write(size=8, is_double=False)

        model._verify_bdf(xref=True)
        save_load_deck(model)

    def test_springs_04(self):
        """tests SPOINT, CELAS3, PELAS, CELAS4"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        s1 = 3
        s2 = 4
        celas3_id = model.add_celas3(eid, pid, [s1, s2], comment='celas3')
        celas3 = model.celas3
        #celas3.raw_fields()
        celas3.write(size=8, is_double=False)

        eid = 2
        k = 1.0e7
        nids = [10, 21]
        celas_id = model.add_celas4(eid, k, nids, comment='celas4')
        celas4 = model.celas4
        #celas4.raw_fields()
        celas4.write(size=8, is_double=False)

        ge = 0.0
        s = 0.
        pelas_id = model.add_pelas(pid, k, ge, s, comment='pelas')
        spoint_id = model.add_spoint([3, 21, 4, 10], comment='spoints')
        spoints = model.spoint
        #spoints.raw_fields()
        spoints.write()

        save_load_deck(model)

    def test_springs_05(self):
        fields = ['CELAS2', '615', '-7.39E6', '613', '1', '613', '5', '0.01']
        model = BDF()
        model.add_card(fields, fields[0])
        model.pop_parse_errors()
        model.validate()
        assert len(model.celas2) == 1, model.celas1


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
