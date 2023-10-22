import copy
import unittest

from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck


class TestContact(unittest.TestCase):

    def test_bsurf(self):
        model = BDF(debug=False)
        bsurf = model.bsurf

        lines = ['BSURF,    1100,    100,     101']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurf.add_card(card)

        lines = ['BSURF,    1100,    11,THRU,15']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurf.add_card(card)

        lines = ['BSURF,    1101,    11,THRU,15,',
                 ',1,2']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurf.add_card(card)

        model.setup()
        size = 8
        bsurf.write(size, 'dummy')
        save_load_deck(model)

    def test_bsurfs(self):
        model = BDF(debug=False)
        bsurfs = model.bsurfs

        lines = ['BSURFS,    1100,    100,     101']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurfs.add_card(card)

        lines = ['BSURFS,    1100,    11,THRU,15']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurfs.add_card(card)

        lines = ['BSURFS,    1101,    11,THRU,15,',
                 ',1,2']
        card = model._process_card(lines)
        card = BDFCard(card)
        card = bsurfs.add_card(card)

        model.setup()
        size = 8
        bsurfs.write(size, 'dummy')
        save_load_deck(model)

    def test_contact_01(self):
        """checks the BSURF cards"""
        model = BDF(debug=False)

        lines = [
            'BSURF          3       1       2       3       4       5       6       7',
            '               8       9      10      11      12      13      14      15',
            '              16      17      18      19      20      21      22      23',
        ]
        unused_card = model.add_card(copy.deepcopy(lines), 'BSURF', is_list=False)
        model.setup()
        out = model.bsurf.write(8, None)
        out_expected = 'BSURF          3       1    THRU      23\n'
        msg = f'out     ={out!r}\nexpected={out_expected!r}'
        assert out == out_expected, msg
        #lines2 = out.split('\n')
        #for line, line2 in zip(lines, lines2):
            #self.assertEqual(line, line2)

    def test_contact_2(self):
        sid = 42
        eids = [1, 2, 3]
        model = BDF(debug=False)
        bsurf = model.bsurf
        sid = 42
        bsurfi = model.add_bsurf(sid, eids, comment='bsurf')
        #bsurfi.raw_fields()

        sid = 43
        g1s = [10, 11, 12]
        g2s = [20, 21, 22]
        g3s = [30, 31, 32]
        #bsurfsi = model.add_bsurfs(sid, eids, g1s, g2s, g3s, comment='bsurfs')
        bsurfsi = model.add_bsurfs(sid, eids, comment='bsurfs')
        #bsurfs.raw_fields()

        contact_set_id = 44
        source_ids = [37, 38]
        target_ids = [47, 48]
        frictions = [0.11, 0.22]
        min_distances = [0.001, 0.001]
        max_distances = [0.1, 0.2]
        HAS_CONTACT = False
        if HAS_CONTACT:
            bctset = model.add_bctset(contact_set_id, source_ids, target_ids, frictions,
                                      min_distances, max_distances,
                                      comment='bctset')
            bctset.raw_fields()

            contract_region = 100
            surface = 'BOT'
            contact_type = 'RIGID'
            offset = .1012
            master_grid_point = 101
            bcrpara = model.add_bcrpara(contract_region, surface, offset, contact_type,
                                        master_grid_point, comment='bcrpara')
            bcrpara.raw_fields()
        model.setup()
        model.validate()

        if HAS_CONTACT:
            contact_region = 102
            params = {'cat' : 111, 'dog' : 222, 'frog' : 0.}
            bctpara = model.add_bctpara(contact_region, params, comment='bctpara')
            bctpara.raw_fields()
            str(bctpara)

            contact_region = 300
            contact_sets = [301, 302]
            bctadd = model.add_bctadd(contact_region, contact_sets, comment='bctadd')
            bctadd.raw_fields()
        save_load_deck(model)

    def _test_contact_3(self):
        """
        tests:
         - BLSEG
         - BCONP -> BFRIC
        """
        model = BDF(debug=False, log=None, mode='msc')
        nodes = [2, 3]
        contact_id = 4
        master = 5
        slave = 6
        sfac = 1.2
        friction_id = 7
        ptype = 8
        cid = 9

        line_id = master
        model.add_blseg(line_id, nodes, comment='blseg_master')

        line_id = slave
        blseg = model.add_blseg(line_id, nodes, comment='blseg_slave')
        bconp = model.add_bconp(contact_id, slave, master, sfac, friction_id, ptype, cid,
                                comment='bconp')
        mu1 = 0.2
        bfric = model.add_bfric(friction_id, mu1, fstiff=None, comment='bfric')
        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [0., 0., 0.])
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1]
        xzplane = [1., 0., 0.]
        model.add_cord2r(9, origin, zaxis, xzplane, rid=0, setup=True, comment='')
        blseg.raw_fields()
        bconp.raw_fields()
        save_load_deck(model)

    def test_contact_bgset(self):
        """
        |   1   |  2   |  3   |   4  |    5    | 6  |  7   |   8  |  9 |
        | BGSET | GSID | SID1 | TID1 | SDIST1  |    | EXT1 |      |    |
        |       |      | SID2 | TID2 | SDIST2  |    | EXT2 |      |    |
        """
        model = BDF(debug=True, log=None, mode='msc')
        gsid = 1
        sids = [1, 2, 3]
        tids = [10, 20, 30]
        sdists = [100., 200., 300.]
        exts = [0.01, 0.02, 0.03]
        card = [
            'BGSET',
            gsid, sids[0], tids[0], sdists[0], None, exts[0], None, None,
            None, sids[1], tids[1], sdists[1], None, exts[1], None, None,
            None, sids[2], tids[2], sdists[2], None, exts[2], None, None,
        ]
        card = model.add_card(card, 'BGSET', comment='comment', ifile=None, is_list=True, has_none=True)
        model.setup()
        save_load_deck(model)
        #bgset = model.add_bgset

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
