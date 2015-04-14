from __future__ import unicode_literals
import unittest

from pyNastran.bdf.bdf import GRID, SPOINTs as SPOINT, BDFCard

class TestNodes(unittest.TestCase):
    def test_grid_01(self):
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        datai = [nid, cp, 0., 0., 0., cd, ps, seid]
        n1 = GRID(data=datai)
        #print n1

        msg = n1.write_card(size=8)
        #print(msg)
        msg = n1.write_card(size=16)
        #print(msg)
        msg = n1.write_card(size=16, is_double=True)
        #print(msg)
        if 0:
            msg = n1.write_card(size=8)
            #print('%r' % msg)
            # small field
            self.assertEqual(msg, 'GRID           1       2      0.      0.      0.                        \n')
            msg = n1.write_card(size=16)

            # large field
            card = ('GRID*                  1               2             .-0             .-0\n'
                    '*                    .-0                                                \n')
            print('%r' % msg)
            ref = 'ERROR\n'
            if card != msg:
                scard = card.split('\n')
                smsg  = msg.split('\n')
                i = 0
                print(scard)
                print(smsg)
                for sc, sm in zip(scard, smsg):
                    if sc != sm:
                        ref += 'i=%s\ncard=%r\nmsg =%r\n' % (i, sc, sm)
                    i += 1
            print(ref)
            self.assertEqual(msg, card), ref

    def test_spoint_01(self):
        #      12345678 2345678 2345678 2345678 2345678 2345678
        msg = 'SPOINT         1       3       5\n'
        card = BDFCard(['SPOINT', 1, 3, 5])
        s1 = SPOINT(card)
        assert list(s1.spoints) == [1, 3, 5], '\n%s' % list(s1.spoints)
        assert s1.write_card() == msg, '\n%s---\n%s' % (s1.write_card(), msg)

        #      12345678 2345678 2345678 2345678 2345678 2345678
        msg = 'SPOINT         1    THRU       5\n'
        card = BDFCard(['SPOINT', 1, 'THRU', 5])
        s2 = SPOINT(card)
        assert list(s2.spoints) == [1, 2, 3, 4, 5], '\n%s' % list(s2.spoints)
        #assert s2.write_card() == msg, '\n%s---\n%s' % (s2.write_card(), msg)

        #       12345678 2345678 2345678 2345678 2345678 2345678
        msg  = 'SPOINT         7\n'
        msg += 'SPOINT         1    THRU       5\n'
        card = BDFCard(['SPOINT', 1, 2, 3, 4, 5, 7])
        s3 = SPOINT(card)
        assert list(s3.spoints) == [1, 2, 3, 4, 5, 7], '\n%s' % list(s3.spoints)
        #assert s3.write_card() == msg, '\n%s---\n%s' % (s3.write_card(), msg)

        #       12345678 2345678 2345678 2345678 2345678 2345678
        msg  = 'SPOINT         7\n'
        msg += 'SPOINT         1    THRU       5\n'
        card = BDFCard(['SPOINT', 1, 'THRU', 5, 7])
        s4 = SPOINT(card)
        assert list(s4.spoints) == [1, 2, 3, 4, 5, 7], '\n%s' % list(s4.spoints)
        #assert s4.write_card() == msg, '\n%s---\n%s' % (s4.write_card(), msg)


        #       12345678 2345678 2345678 2345678 2345678 2345678
        msg  = 'SPOINT         7\n'
        msg += 'SPOINT         1    THRU       5\n'
        card = BDFCard(['SPOINT', 1, 'THRU', 5, 7])
        s5 = SPOINT(card)
        assert list(s5.spoints) == [1, 2, 3, 4, 5, 7], '\n%s' % list(s5.spoints)
        assert s5.write_card() == msg, '\n%s---\n%s' % (s5.write_card(), msg)




if __name__ == '__main__':  # pragma: no cover
    unittest.main()
