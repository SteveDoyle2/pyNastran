"""various solid element tests"""
import copy
import unittest

import numpy as np

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.cards.elements.solid import (
    #CTETRA4, CHEXA8, CPENTA6,
    #CTETRA10, CHEXA20,
    CPENTA15
)
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  #mass_properties_breakdown
from pyNastran.bdf.mesh_utils.loads import sum_forces_moments, sum_forces_moments_elements


class TestSolids(unittest.TestCase):
    """various solid element tests"""
    def test_cpenta_01(self):
        """tests a cpenta15"""
        lines = [
            'CPENTA,85,22,201,202,203,205,206,207,+PN2',
            '+PN2,209,210,217,  ,  ,  ,213,214,',
            ',218'
        ]
        bdf = BDF(debug=False)
        card = bdf._process_card(lines)
        card = BDFCard(card)

        solid = CPENTA15.add_card(card, comment='cpenta15')
        solid.write_card(size=8, is_double=False)
        solid.write_card(size=16, is_double=False)
        solid.raw_fields()

        node_ids = solid.node_ids
        assert node_ids == [201, 202, 203, 205, 206, 207,
                            209, 210, 217, None, None, None, 213, 214, 218], node_ids
        nids = [201, 202, 203, 205, 206, 207,
                209, 210, 217, None, None, None, 213, 214, 218]
        CPENTA15.add_card(card, comment='spike')
        eid = 85
        pid = 22
        bdf.add_cpenta(eid, pid, nids, comment='spike')

    #def test_cpenta_01b(self):
        #pass
          # ..todo:: this is wrong...
        #lines = [  # this will fail!
        #    'CPENTA,85,22,201,202,203,205,206,207,+PN2',
        #    '+PN2,209,210,217,  ,  ,  ,213,214,218'
        #]

    #def test_solid_02(self):
        #"""checks nonlinear static solid material"""
        #pass

    def test_solid_01(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],
            #['CPENTA',19, pid+1, 21, 22, 23, 24, 25, 26],

            # static
            #['PIHEX', pid+1, mid],
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho]
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # CTETRA
        eid = 8
        mid = 2
        pid = 4
        nsm = 0.
        volume = 1. / 3.
        rho = 0.1
        self.check_solid(model, eid, 'CTETRA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, volume)

        eid = 9
        volume = 1.0
        self.check_solid(model, eid, 'CPENTA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, volume)

        eid = 7
        volume = 2.0
        self.check_solid(model, eid, 'CHEXA', pid, 'PSOLID', mid, 'MAT1', nsm, rho, volume)

    def test_solid_02(self):
        """tests CHEXA, CTETRA, CPENTA, PLSOLID, MATHP"""
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # hyperelastic
            ['PLSOLID', pid, mid, 'GRID'],
            ['MATHP', mid, None, None, None, rho],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        # CTETRA
        eid = 8
        nsm = 0.
        volume = 1. / 3.
        self.check_solid(model, eid, 'CTETRA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, volume)

        eid = 9
        volume = 1.0
        self.check_solid(model, eid, 'CPENTA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, volume)

        eid = 7
        volume = 2.0
        self.check_solid(model, eid, 'CHEXA', pid, 'PLSOLID', mid, 'MATHP', nsm, rho, volume)

    def test_solid_03(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CIHEX1', 17, pid+1, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CIHEX2', 18, pid+1, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['PIHEX', pid+1, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho],
            ['MATS1', mid, None, 'PLASTIC', 0.0, 1, 1, 100000., ],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()
        save_load_deck(model)

    def test_solid_04(self):
        """checks linear static solid material"""
        mid = 2
        pid = 4
        rho = 0.1
        table_id = 42
        cards = [
            #$ Solid Nodes
            ['GRID', 11, 0, 0., 0., 0., 0],
            ['GRID', 12, 0, 1., 0., 0., 0],
            ['GRID', 13, 0, 1., 1., 0., 0],
            ['GRID', 14, 0, 0., 1., 0., 0],

            ['GRID', 15, 0, 0., 0., 2., 0],
            ['GRID', 16, 0, 1., 0., 2., 0],
            ['GRID', 17, 0, 1., 1., 2., 0],
            ['GRID', 18, 0, 0., 1., 2., 0],

            # Solids
            ['CHEXA', 7, pid, 11, 12, 13, 14, 15, 16, 17, 18],
            ['CTETRA', 8, pid, 11, 12, 13, 15],

            # Solid Nodes
            ['GRID', 21, 0, 0., 0., 0., 0,],
            ['GRID', 22, 0, 1., 0., 0., 0,],
            ['GRID', 23, 0, 1., 1., 0., 0,],
            ['GRID', 24, 0, 0., 0., 2., 0,],
            ['GRID', 25, 0, 1., 0., 2., 0,],
            ['GRID', 26, 0, 1., 1., 2., 0,],
            ['CPENTA', 9, pid, 21, 22, 23, 24, 25, 26],

            # static
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho],
            ['MATS1', mid, table_id, 'PLASTIC', 0.0, 1, 1, 100000., ],
            #['TABLEST'],
            ['TABLES1', table_id, 1, None, None, None, None, None, None,
             1.0, 10.0, 2.0, 10.0, 'ENDT'],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        model.cross_reference()

        mat = model.Material(mid)
        mat.E()

        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=False)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        save_load_deck(model)

    def test_solids_ctetra4(self):
        """tests a CTETRA4"""
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model = BDF(debug=False)
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(15, [0., 0., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=0.1)
        nids = [11, 12, 13, 15]
        model.add_ctetra(eid, pid, nids, comment='ctetra')

        eid2 = eid + 1
        pid2 = pid + 1
        model.add_ctetra(eid2, pid2, nids, comment='ctetra')

        global_ply_ids = [1, 2, 3]
        nlayers = len(global_ply_ids)
        thicknesses = [0.1] * nlayers
        mids = [mid] * nlayers
        thetas = [0.] * nlayers
        pcomps = model.add_pcomps(pid2, global_ply_ids, mids, thicknesses, thetas,
                                  cordm=0, psdir=13, sb=None, nb=None, tref=0.0, ge=0.0,
                                  failure_theories=None, interlaminar_failure_theories=None,
                                  souts=None, comment='pcomps')
        pcomps.raw_fields()
        end_checks(model)

        model.cross_reference()
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=False)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        save_load_deck(model)

    def test_solids_ctetra4_mat9(self):
        """tests a CTETRA4"""
        eid = 10
        pid = 20
        mid = 30
        #E = 3.e7
        #G = None
        #nu = 0.3
        model = BDF(debug=False)
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(15, [0., 0., 2.])
        model.add_psolid(pid, mid)
        mat9 = model.add_mat9(mid, G11=1000., rho=0.2)
        mat9.raw_fields()

        nids = [11, 12, 13, 15]
        model.add_ctetra(eid, pid, nids, comment='ctetra')

        end_checks(model)

        model.cross_reference()
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=False)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

    def test_solids_ctetra10(self):
        """tests a CTETRA10"""
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model = BDF(debug=False)
        g110 = model.add_grid(110, [0., 0., 0.])
        g120 = model.add_grid(120, [1., 0., 0.])
        g130 = model.add_grid(130, [1., 1., 0.])
        g140 = model.add_grid(140, [0., 0., 2.])

        model.add_grid(111, g110.xyz+g120.xyz)
        model.add_grid(112, g120.xyz+g130.xyz)
        model.add_grid(113, g130.xyz+g110.xyz)

        model.add_grid(121, g110.xyz+g140.xyz)
        model.add_grid(122, g120.xyz+g140.xyz)
        model.add_grid(123, g130.xyz+g140.xyz)

        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [
            110, 120, 130, 140,
            111, 112, 113,
            121, 122, 123
        ]
        model.add_ctetra(eid, pid, nids, comment='ctetra10')
        end_checks(model)

        with self.assertRaises(RuntimeError):
            model.get_length_breakdown(property_ids=None, stop_if_no_length=True)
        with self.assertRaises(RuntimeError):
            model.get_area_breakdown(property_ids=None, stop_if_no_area=True)

        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=False)
        model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True, detailed=True)
        model.get_volume_breakdown(property_ids=None, stop_if_no_volume=True)

        save_load_deck(model)

    def test_solids_cpyram5(self):
        """tests a CPYRAM5"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [1., 0., 0.])
        model.add_grid(30, [1., 1., 0.])
        model.add_grid(40, [0., 1., 0.])
        model.add_grid(50, [0., 0., 1.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [10, 20, 30, 40, 50]
        model.add_cpyram(eid, pid, nids, comment='cpyram')

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_cpyram(eid, pid, nids2, comment='cpyram13')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)

        # quad face; g1/g3 are opposite
        # tri face: g1/g3 are adjacent
        #model.add_pload4(sid, eids, pressures, g1=None, g34=None, cid=0, nvector=None,
                         #surf_or_line='SURF', line_load_dir='NORM', comment='')

        # --------------------------------------------------------------------------------------------
        eid = 10
        pressures = 10.
        # -------------------
        # quad faces
        sid = 10
        g1 = 10
        g3 = 30
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 11
        g1 = 30
        g3 = 10
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        # -------------------
        # tri faces
        sid = 21
        g1 = 10
        g3 = 20
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 22
        g1 = 20
        g3 = 30
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 23
        g1 = 30
        g3 = 40
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 24
        g1 = 40
        g3 = 10
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        # -------------------
        # tri faces top
        g3 = 50

        sid = 31
        g1 = 10
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 32
        g1 = 20
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 33
        g1 = 30
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        sid = 34
        g1 = 40
        model.add_pload4(sid, eid, pressures, g1=g1, g34=g3, cid=0, nvector=None,
                         surf_or_line='SURF', line_load_dir='NORM', comment='')

        # --------------------------------------------------------------------------------------------
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        eids = None
        nids = None
        p0 = [0., 0., 0.]
        # --------------------------------------------------------------------------------------------
        def sum_loads(loadcase_id, force):
            force1, moment1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False, xyz_cid0=None)
            force2, moment2 =  sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                           include_grav=False, xyz_cid0=None)
            assert np.allclose(force1, force), 'sid=%s force1=%s force2=%s' % (loadcase_id, force1, force2)
            assert np.allclose(force2, force), 'sid=%s force1=%s force2=%s' % (loadcase_id, force1, force2)

        sum_loads(10, [0., 0., 10.])
        sum_loads(11, [0., 0., 10.])

        sum_loads(21, [0., 5., 0.])
        sum_loads(22, [-5., 0., -5.])
        sum_loads(23, [0., -5., -5.])
        sum_loads(24, [-5., 0., 0.])

        with self.assertRaises(RuntimeError):
            sum_loads(31, [0., 5., 0.])
        with self.assertRaises(RuntimeError):
            sum_loads(32, [-5., 0., -5.])
        with self.assertRaises(RuntimeError):
            sum_loads(33, [0., -5., -5.])
        with self.assertRaises(RuntimeError):
            sum_loads(34, [-5., 0., 0.])

        # -------------------


        end_checks(model)
        save_load_deck(model)

    def test_solids_cpenta(self):
        """tests a CPENTA6"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(21, [0., 0., 0.])
        model.add_grid(22, [1., 0., 0.])
        model.add_grid(23, [1., 1., 0.])
        model.add_grid(24, [0., 0., 2.])
        model.add_grid(25, [1., 0., 2.])
        model.add_grid(26, [1., 1., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [21, 22, 23, 24, 25, 26]
        model.add_cpenta(eid, pid, nids, comment='cpenta')

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_cpenta(eid, pid, nids2, comment='cpenta15')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)
        end_checks(model)
        save_load_deck(model)

    def test_solids_chexa(self):
        """tests a CHEXA8"""
        model = BDF(debug=False)
        eid = 10
        pid = 20
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])

        model.add_grid(15, [0., 0., 2.])
        model.add_grid(16, [1., 0., 2.])
        model.add_grid(17, [1., 1., 2.])
        model.add_grid(18, [0., 1., 2.])
        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [11, 12, 13, 14, 15, 16, 17, 18]
        elem = model.add_chexa(eid, pid, nids, comment='chexa')
        elem.write_card(size=8)
        elem.write_card(size=16)
        elem.write_card_16(is_double=False)

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_chexa(eid, pid, nids2, comment='chexa16')
        elem2.write_card(size=8)
        elem2.write_card(size=16)
        elem2.write_card_16(is_double=False)

        end_checks(model)
        elem = model.elements[eid]
        assert elem.Mass() > 0, elem.Mass()
        save_load_deck(model)

    def check_solid(self, model, eid, etype, pid, ptype, mid, mtype, nsm, rho, volume):
        """checks that various solid methods work"""
        mass = rho * volume
        element = model.elements[eid]
        assert max(element.node_ids) > 0
        assert pid in model.properties, 'pid is missing for\n%s' % str(element)
        self.assertEqual(element.type, etype)
        self.assertEqual(element.eid, eid)
        self.assertEqual(element.pid_ref.type, ptype)
        self.assertEqual(element.Pid(), pid)
        self.assertEqual(element.pid_ref.mid_ref.type, mtype)
        self.assertEqual(element.Mid(), mid)
        self.assertEqual(element.Volume(), volume)
        self.assertEqual(element.Mass(), mass)
        mass_mp = mass_properties(model, element_ids=eid)[0]
        mass_mp_nsm = mass_properties_nsm(model, element_ids=eid)[0]
        unused_centroid, unused_xe, unused_ye, unused_ze = element.material_coordinate_system()
        assert np.allclose(mass, mass_mp)
        assert np.allclose(mass, mass_mp_nsm)


def end_checks(model):
    """various checks"""
    model.validate()
    model._verify_bdf(xref=False)
    model.cross_reference()
    model._verify_bdf(xref=True)
    model.uncross_reference()
    model.cross_reference()
    model.pop_xref_errors()
    mass, cg, inertia = mass_properties(model)
    assert mass > 0, 'mass=%s, cg=%s, inertia=%s' % (mass, cg, inertia)

    #bdf_filename = 'solid_test.bdf'
    #model.write_bdf(bdf_filename)
    #read_bdf(bdf_filename, debug=False)
    #os.remove(bdf_filename)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
