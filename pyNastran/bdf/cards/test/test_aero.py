import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF, CORD2R
from pyNastran.bdf.cards.aero import (
    FLFACT, AEFACT, AEPARM, AERO, AEROS, CAERO1, CAERO2, PAERO1, PAERO2)

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

comment_bad = 'this is a bad comment'
comment_good = 'this is a good comment\n'
model = BDF(debug=None)
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
    def test_aefact_1(self):
        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', '98', '.3', '0.7', '1.0']
        model.add_card(data, data[0], comment_good, is_list=True)

        msg = 'this is a bad commentAEFACT        97      .3      .7      1.\n'
        aefact97 = model.aefacts[97]
        aefact98 = model.aefacts[98]
        self.assertTrue(all(aefact97.Di == [.3, .7, 1.0]))
        self.assertTrue(all(aefact98.Di == [.3, .7, 1.0]))

        out = aefact97.write_card(8, None)
        self.assertEqual(msg, out)

        msg = 'this is a good comment\nAEFACT        98      .3      .7      1.\n'
        out = aefact98.write_card(8, None)
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

        Di = [1., 2., 3.]
        aefact = AEFACT(200, Di, comment='')

   # def test_aelink_1(self):
    def test_aelist_1(self):
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1201, 1202]
        aelist75 = model.aelists[75]
        #print(aelist.elements)
        #print(elements)
        self.assertTrue(elements == aelist75.elements)

        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1108, 1202]
        data = ['AELIST', 76, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1108, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        aelist76 = model.aelists[76]
        #print(aelist76 .elements)
        #print(elements)
        self.assertFalse(elements == aelist76.elements)

        elements = list(set(elements))
        elements.sort()
        self.assertTrue(elements == aelist76.elements)

    def test_aeparm_1(self):
        aeparm = AEPARM(100, 'THRUST', 'lb')
        aeparm.write_card()

   # def test_aestat_1(self):
   # def test_aesurf_1(self):
   # def test_aesurfs_1(self):

    def test_aero_1(self):
        acsid = 0.
        velocity = None
        cref = 1.0
        rho_ref = 1.0
        aero = AERO(acsid, velocity, cref, rho_ref, sym_xz=0., sym_xy=0,
                    comment='')
        aero.write_card()

    def test_aeros_1(self):
        acsid = 0.
        velocity = None
        cref = 1.0
        bref = 2.0
        sref = 100.
        aeros = AEROS(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0,
                      comment='')
        aeros.write_card()


    def test_caero1_1(self):
        eid = 1
        pid = 10
        cp = 4
        nspan = None
        lspan = 3
        nchord = None
        lchord = 4
        igid = None
        p1 = [0., 0., 0.]
        x12 = 5.
        p4 = [2., 3., 4.]
        x43 = 1.
        caero = CAERO1(eid, pid, cp, nspan, lspan, nchord, lchord, igid, p1,
                       x12, p4, x43)
        caero.validate()
        caero.write_card()

    def test_caero2_1(self):
        model = BDF()
        eid = 1
        pid = 10
        cp = 4
        nsb = 0
        nint = 0
        lsb = 3
        lint = 6
        igid = None
        p1 = [0., 1., 2.]
        x12 = 10.
        caero = CAERO2(eid, pid, cp, nsb, nint, lsb, lint, igid, p1, x12,
                       comment='')
        caero.validate()
        aefact = AEFACT(lint, [0., 1., 2., 3., 4., 5.])
        aefact.validate()
        model.aefacts[lint] = aefact

        orient = 'Z'
        width = 10.
        AR = 2.
        lrsb = 0
        lrib = 3
        lth1 = 0
        lth2 = 0
        thi = 0
        thn = 0
        paero = PAERO2(pid, orient, width, AR, lrsb, lrib, lth1, lth2, thi, thn)
        paero.validate()
        model.paeros[pid] = paero

        coord = CORD2R(cp, rid=0, origin=None, zaxis=None, xzplane=None,
                       comment='')
        coord.validate()
        model.coords[cp] = coord

        aefact = AEFACT(lrib, [0., 1., 2., 3., 4., 5.])
        aefact.validate()
        model.aefacts[lrib] = aefact

        acsid = 0
        velocity = None
        cref = 1.0
        rho_ref = 1.0

        aero = AERO(acsid, velocity, cref, rho_ref,
                          comment='')
        aero.validate()
        model.aero = aero

        paero.cross_reference(model)
        caero.cross_reference(model)
        xyz, elems = caero.get_points_elements_3d()

   # def test_caero3_1(self):
   # def test_caero4_1(self):
   # def test_caero5_1(self):

   # def test_paero1_1(self):
   # def test_paero2_1(self):
   # def test_paero3_1(self):
   # def test_paero4_1(self):
   # def test_paero5_1(self):

   # def test_spline1_1(self):
   # def test_spline2_1(self):
   # def test_spline3_1(self):
   # def test_spline4_1(self):
   # def test_spline5_1(self):


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
