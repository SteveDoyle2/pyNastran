import copy
import unittest

from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck

#from pyNastran.bdf.field_writer_8 import print_card

class TestContact(unittest.TestCase):

    def test_contact_01(self):
        """checks the BSURF cards"""
        model = BDF(debug=False)

        lines = [
            'BSURF          3       1       2       3       4       5       6       7',
            '               8       9      10      11      12      13      14      15',
            '              16      17      18      19      20      21      22      23',
        ]
        unused_card = model.add_card(copy.deepcopy(lines), 'BSURF', is_list=False)
        out = model.bsurf[3].write_card(size=8, is_double=False)
        lines2 = out.split('\n')
        for line, line2 in zip(lines, lines2):
            self.assertEqual(line, line2)

    def test_contact_2(self):
        eids = [1, 2, 3]
        model = BDF(debug=False)
        sid = 42
        bsurf = model.add_bsurf(sid, eids, comment='bsurf')
        bsurf.raw_fields()

        sid = 43
        g1s = [10, 11, 12]
        g2s = [20, 21, 22]
        g3s = [30, 31, 32]
        bsurfs = model.add_bsurfs(sid, eids, g1s, g2s, g3s, comment='bsurfs')
        bsurfs.raw_fields()

        contact_set_id = 44
        source_ids = [37, 38]
        target_ids = [47, 48]
        frictions = [0.11, 0.22]
        min_distances = [0.001, 0.001]
        max_distances = [0.1, 0.2]
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
        model.validate()

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

    def test_contact_3(self):
        """
        tests:
         - BLSEG
         - BCONP -> BFRIC
        """
        log = SimpleLogger(level='warning')
        model = BDF(log=log, mode='msc')
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
        bconp = model.add_bconp(
            contact_id, slave, master, sfac, friction_id, ptype, cid,
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
        log = SimpleLogger(level='warning')
        model = BDF(log=log, mode='msc')
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
        save_load_deck(model)
        #bgset = model.add_bgset

    def test_contact_bcbody(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log, mode='msc')
        contact_id = 42
        bsid = 3
        bcbody = model.add_bcbody(
            contact_id, bsid, {},
            dim='3D', behav='DEFORM',
            istype=0, fric=0,
            idispl=0, comment='bcbody')
        bcbody.cross_reference(model)
        bcbody.raw_fields()
        bcbody.write_card(size=8)
        bcbody.write_card(size=16)
        contact_id += 1
        grow = [0.1, 0.2, 0.3, 1, 2, 3]
        bcbody_grow = model.add_bcbody(
            contact_id, bsid, {'GROW': grow},
            dim='3D', behav='DEFORM',
            istype=0, fric=0,
             idispl=0, comment='bcbody')
        contact_id += 1

        bcbody_rigid = model.add_bcbody(
            contact_id, bsid, {},
            dim='3D', behav='RIGID',
            istype=0, fric=0,
             idispl=0, comment='bcbody')

        heat = [0.1] * 17
        bcbody_heat = model.add_bcbody(
            contact_id, bsid, {'HEAT': heat},
            dim='3D', behav='SYMM',
            istype=0, fric=0,
            idispl=0, comment='bcbody')

        csid = 10
        params = {
            'BIAS': 20.0,
            'DDULMT': 1.0,
            'ERRBAS': 1,
            'ERROR': 3.0,
            'FKIND': 1,
            'FNTOL': 4.0,
            'FTYPE': 1,
            'IBSEP': 1,
            'ISPLIT': 1,
            'ICSEP': 1,
            'LINQUAD': 1,
            'LINCNT': -1,
            'METHOD': 'NODESURF',
            'MAXENT': 3,
            'MAXNOD': 42,
            'NLGLUE': 0,
            'NBODIES': 2,
            'SLDLMT': 1.0,
            'SEGSYM': 0,
            'THKOFF': 1,
        }
        bcpara = model.add_bcpara(csid, params, comment='bcpara')
        bcpara.raw_fields()
        bcpara.write_card(size=8)
        bcpara.write_card(size=16)

        bct_params1 = {
            'TYPE': 1,
            'NSIDE': 2,
            'TBIRTH': 3.0,
            'TDEATH': 4.0,
            'INIPENE': 2,
            'PDEPTH': 6.0,
        }
        bct_params2 = {
            'SEGNORM': -1,
            'OFFTYPE': 1,
            'OFFSET': 9.0,
            'TZPENE': 10.0,
            'CSTIFF': 1,
        }
        bct_params3 = {
            'TIED': 1,
            'TIEDTOL': 13.0,
            'EXTFAC': 0.09,
        }
        bctpara1 = model.add_bctpara(csid, bct_params1, comment='bctpara')
        bctpara2 = model.add_bctpara(csid+1, bct_params2, comment='bctpara')
        bctpara3 = model.add_bctpara(csid+2, bct_params3, comment='bctpara')
        bctpara1.raw_fields()
        bctpara1.write_card(size=8)
        bctpara1.write_card(size=16)
        save_load_deck(model)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
