from __future__ import unicode_literals
import unittest

from pyNastran.bdf.bdf import GRID

class TestNodes(unittest.TestCase):
    def test_nodes_01(self):
        nid = 1
        cp = 2
        cd = 0
        ps = ''
        seid = 0
        datai = [nid, cp, 0., 0., 0., cd, ps, seid]
        n1 = GRID(data=datai)
        #print n1

        msg = n1.write_bdf2(size=8)
        print msg
        msg = n1.write_bdf2(size=16)
        print msg
        msg = n1.write_bdf2(size=16, double=True)
        print msg
        if 0:
            msg = n1.write_bdf2(size=8)
            #print '%r' % msg
            # small field
            self.assertEqual(msg, 'GRID           1       2      0.      0.      0.                        \n')
            msg = n1.write_bdf2(size=16)

            # large field
            card = ('GRID*                  1               2             .-0             .-0\n'
                    '*                    .-0                                                \n')
            print '%r' % msg
            ref = 'ERROR\n'
            if card != msg:
                scard = card.split('\n')
                smsg  = msg.split('\n')
                i = 0
                print scard
                print smsg
                for sc, sm in zip(scard, smsg):
                    if sc != sm:
                        ref += 'i=%s\ncard=%r\nmsg =%r\n' % (i, sc, sm)
                    i += 1
            print ref
            self.assertEqual(msg, card), ref


if __name__ == '__main__':
    unittest.main()
