from io import StringIO
import unittest
from itertools import count

from pyNastran.dev.bdf_vectorized.bdf import BDF, BDFCard # GRID,

class TestNodes(unittest.TestCase):
    def test_grid_01(self):
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        card_count = {'GRID': 1,}

        model = BDF(debug=False)
        model.allocate(card_count)
        data1 = BDFCard(['GRID', nid, cp, 0., 0., 0., cd, ps, seid])

        nodes = model.grid
        nodes.add(data1)

        #print(n1)
        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        nodes.write_card(f, size=16, write_header=False)
        nodes.write_card(f, size=16, is_double=True, write_header=False)

        # small field
        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        msg = f.getvalue()
        card = 'GRID           1       2      0.      0.      0.\n'
        self.assertCardEqual(msg, card)

        # large field
        f = StringIO()
        nodes.write_card(f, size=16, write_header=False)
        card = ('GRID*                  1               2              0.              0.\n'
                '*                     0.\n')
        msg = f.getvalue()
        self.assertCardEqual(msg, card)

    def test_grid_02(self):
        model = BDF(debug=False)
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        card_count = {'GRID': 2,}

        model = BDF(debug=False)
        model.allocate(card_count)

        nodes = model.grid
        data1 = BDFCard(['GRID', nid, cp, 0., 0., 0., cd, ps, seid])
        data2 = BDFCard(['GRID', nid+1, cp, 0., 0., 0., cd, ps, seid])
        nodes.add(data1)
        nodes.add(data2)
        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        #print(f.getvalue())

    def test_grid_03(self):
        model = BDF(debug=False)
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        card_count = {'GRID': 2,}

        model = BDF(debug=False)
        model.allocate(card_count)

        nodes = model.grid
        data1 = BDFCard(['GRID', nid, cp, 0., 0., 0., cd, ps, seid])
        data2 = BDFCard(['GRID', nid+1, cp, 0., 0., 0., cd, ps, seid])
        data3 = BDFCard(['GRID', nid+2, cp, 0., 0., 0., cd, ps, seid])
        nodes.add(data1)
        nodes.add(data2)
        nodes.resize(3, refcheck=False)
        #print('nodes.node_id = %s' % nodes.node_id)
        nodes.add(data3)
        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        #print(f.getvalue())

    def test_grid_04(self):
        model = BDF(debug=False)
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        card_count = {'GRID': 2,}

        model = BDF(debug=False)
        model.allocate(card_count)

        nodes = model.grid
        data1 = BDFCard(['GRID', nid, cp, 0., 0., 0., cd, ps, seid])
        data2 = BDFCard(['GRID', nid+1, cp, 0., 0., 0., cd, ps, seid])
        data3 = BDFCard(['GRID', nid+2, cp, 0., 0., 0., cd, ps, seid])
        data4 = BDFCard(['GRID', nid+3, cp, 0., 0., 0., cd, ps, seid])
        nodes.add(data1)
        nodes.add(data2)
        nodes.resize(4, refcheck=False)
        #print('nodes.node_id = %s' % nodes.node_id)
        nodes.add(data3)
        self.assertEqual(len(nodes.node_id), 4)
        self.assertEqual(nodes.n, 4)
        self.assertEqual(nodes.i, 3)
        nodes.shrink(refcheck=False)
        #print('nodes.node_id = %s' % nodes.node_id)
        self.assertEqual(len(nodes.node_id), 3)
        self.assertEqual(nodes.n, 3)
        self.assertEqual(nodes.i, 3)

        nodes.resize(4, refcheck=False)
        nodes.add(data4)
        self.assertEqual(len(nodes.node_id), 4)
        self.assertEqual(nodes.n, 4)
        self.assertEqual(nodes.i, 4)

        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        #print(f.getvalue())

        nodes.resize(2, refcheck=False)
        self.assertEqual(len(nodes.node_id), 2)
        self.assertEqual(nodes.n, 2)
        self.assertEqual(nodes.i, 2)

        f = StringIO()
        nodes.write_card(f, size=8, write_header=False)
        #print(f.getvalue())

    def test_spoint_01(self):
        model = BDF(debug=False)
        card_count = {'SPOINT': 2,}

        model = BDF(debug=False)
        model.allocate(card_count)

        spoints = model.spoint
        data1 = BDFCard(['SPOINT', 1, 2])
        data2 = BDFCard(['SPOINT', 4, 5])
        data3 = BDFCard(['SPOINT', 10, 20])
        data4 = BDFCard(['SPOINT', 30, 'THRU', 40])
        spoints.add(data1)
        spoints.add(data2)

        f = StringIO()
        spoints.write_card(f, size=8)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(len(spoints), 4)
        self.assertEqual(spoints.n, 4)

        spoints.add(data3)
        f = StringIO()
        spoints.write_card(f, size=8)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(len(spoints), 6)
        self.assertEqual(spoints.n, 6)

        spoints.add(data4)
        f = StringIO()
        spoints.write_card(f, size=16)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(len(spoints), 17)
        self.assertEqual(spoints.n, 17)

        spoints.add(data4)
        f = StringIO()
        spoints.write_card(f, size=16)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(spoints.max(), 40)
        self.assertEqual(spoints.min(), 1)
        self.assertEqual(spoints.n, 17)

        if 1 in spoints:
            pass
        else:
            raise RuntimeError('invalid catch')

        if 29 in spoints:
            raise RuntimeError('invalid catch')

        data5 = BDFCard(['SPOINT', 50, 'THRU', 55, 59, '61', 'THRU', '64'])
        spoints.add(data5)
        f = StringIO()
        spoints.write_card(f, size=16)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(spoints.max(), 64)
        self.assertEqual(spoints.min(), 1)
        self.assertEqual(spoints.n, 28)

        spoints.remove([1, 2, 64])
        #self.assertRaises(KeyError, spoints.remove([1, 2]))
        f = StringIO()
        spoints.write_card(f, size=16)
        #print('spoints %s' % f.getvalue())
        self.assertEqual(spoints.n, 25)

        #self.assertRaises(KeyError, spoints.remove([64]))


    def assertCardEqual(self, msg, card):
        if card != msg:
            ref = '\nERROR\n'
            scard = card.split('\n')
            smsg = msg.split('\n')
            #print(scard)
            #print(smsg)
            for i, sc, sm in zip(count(), scard, smsg):
                if sc != sm:
                    ref += 'i=%s\n  card=%r\n  msg =%r\n' % (i, sc, sm)
            print(ref)
            self.assertEqual(msg, card), ref


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
