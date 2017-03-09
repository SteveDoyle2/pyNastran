from __future__ import print_function
import unittest
from six.moves import StringIO

from pyNastran.bdf.bdf import BDF, BDFCard, PDAMP, read_bdf#, get_logger2

class TestDampers(unittest.TestCase):
    def test_damper_01(self):
        """tests PDAMP"""
        lines = ['pdamp, 201, 1.e+5']
        model = BDF(debug=False)
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = PDAMP.add_card(card)
        elem.write_card(size, 'dummy')
        elem.raw_fields()
        self.assertEqual(elem.Pid(), 201)
        self.assertEqual(elem.B(), 1e5)

        fields = ['pelas', 201, 1.e+5, None, None, 202, 2.e+5]
        card_name = fields[0]
        model.add_card(fields, card_name, comment='', is_list=True,
                       has_none=True)
        assert len(model.properties) == 2, model.properties

    def test_damper_02(self):
        """tests CELAS1, CELAS2, PDAMP, PDAMPT"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        nids = [3, 4]
        c1 = 1
        c2 = 1
        celas1 = model.add_cdamp1(eid, pid, nids, c1, c2, comment='cdamp1')
        celas1.raw_fields()
        celas1.write_card(size=8, is_double=False)

        b = 1.0e7
        pdamp = model.add_pdamp(pid, b, comment='pdamp')
        pdamp.raw_fields()
        pdamp.write_card(size=8, is_double=False)

        tbid = 10
        pdampt = model.add_pdampt(pid, tbid, comment='pdampt')
        pdampt.raw_fields()
        pdampt.write_card(size=8, is_double=False)

        eid = 5
        cdamp2 = model.add_cdamp2(eid, b, nids, comment='cdamp2')
        cdamp2.raw_fields()
        cdamp2.write_card(size=8, is_double=False)

        model.add_grid(3, xyz=[0., 0., 0.])
        model.add_grid(4, xyz=[0., 0., 0.])
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)
        model2 = read_bdf(bdf_file, punch=True, debug=False)

    def test_damper_03(self):
        model = BDF(debug=False)
        eid = 1
        pid = 2
        s1 = 3
        s2 = 4
        cdamp3 = model.add_cdamp3(eid, pid, [s1, s2], comment='cdamp3')
        cdamp3.raw_fields()
        cdamp3.write_card(size=8, is_double=False)

        b = 1.0e7
        pelas = model.add_pdamp(pid, b, comment='pdamp')
        spoints = model.add_spoint([3, 4], comment='spoints')
        spoints.raw_fields()
        spoints.write_card()

        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)
        model2 = read_bdf(bdf_file, punch=True, debug=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
