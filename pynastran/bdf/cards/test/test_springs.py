import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, PELAS

bdf = BDF(debug=False)
class TestSprings(unittest.TestCase):
    def test_pelas_01(self):
        lines = ['pelas, 201, 1.e+5']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = PELAS.add_card(card)
        elem.write_card(size, 'dummy')
        elem.raw_fields()
        self.assertEqual(elem.Pid(), 201)
        self.assertEqual(elem.K(), 1e5)

    def test_pelas_02(self):
        fields = ['pelas', 201, 1.e+5, None, None, 202, 2.e+5]
        model = BDF(debug=False)
        #model.echo = True
        card_name = fields[0]
        model.add_card(fields, card_name, comment='', is_list=True,
                      has_none=True)
        assert len(model.properties) == 2, model.properties


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
