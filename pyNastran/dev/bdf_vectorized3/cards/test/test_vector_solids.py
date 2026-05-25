"""various solid element tests"""
import copy
import unittest

import numpy as np

from pyNastran.bdf.bdf import BDFCard # BDF,
from pyNastran.dev.bdf_vectorized3.bdf import BDF

from pyNastran.bdf.cards.elements.solid import (
    #CTETRA4, CHEXA8, CPENTA6,
    #CTETRA10, CHEXA20,
    CPENTA15
)
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
#from pyNastran.bdf.mesh_utils.mass_properties import (
    #mass_properties, mass_properties_nsm)  #mass_properties_breakdown
#from pyNastran.bdf.mesh_utils.loads import sum_forces_moments, sum_forces_moments_elements


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
            ['PSOLID', pid, mid, 0],
            ['MAT1', mid, 1.0, 2.0, 3.0, rho]
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        #model.cross_reference()
        model.setup()

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
        #model.cross_reference()
        model.setup()

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
            ['MATS1', mid, None, 'PLASTIC', 0.0, 1, 1, 100000., ],
        ]
        model = BDF(debug=False)
        for fields in cards:
            model.add_card(fields, fields[0], is_list=True)
        #model.cross_reference()
        model.setup()
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
        #model.cross_reference()
        model.setup()
        mat = model.Material(mid)[0]
        mat.E

        model.get_mass_breakdown() # property_ids=None, stop_if_no_mass=True
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
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
        #pcomps.raw_fields()
        model.setup()
        end_checks(model)

        #model.cross_reference()
        model.get_mass_breakdown(stop_if_no_mass=True)
        #model.get_mass_breakdown(property_ids=None, stop_if_no_mass=True)
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
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
        model.add_mat9(mid, G11=1000., rho=0.2)
        #mat9.raw_fields()

        nids = [11, 12, 13, 15]
        model.add_ctetra(eid, pid, nids, comment='ctetra')
        model.setup()
        end_checks(model)

        #model.cross_reference()
        model.get_mass_breakdown() # property_ids=None,, stop_if_no_mass=True
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
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

        g110_xyz = np.array([0., 0., 0.])
        g120_xyz = np.array([1., 0., 0.])
        g130_xyz = np.array([1., 1., 0.])
        g140_xyz = np.array([0., 0., 2.])

        g110 = model.add_grid(110, g110_xyz)
        g120 = model.add_grid(120, g120_xyz)
        g130 = model.add_grid(130, g130_xyz)
        g140 = model.add_grid(140, g140_xyz)

        model.add_grid(111, g110_xyz+g120_xyz)
        model.add_grid(112, g120_xyz+g130_xyz)
        model.add_grid(113, g130_xyz+g110_xyz)

        model.add_grid(121, g110_xyz+g140_xyz)
        model.add_grid(122, g120_xyz+g140_xyz)
        model.add_grid(123, g130_xyz+g140_xyz)

        model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [
            110, 120, 130, 140,
            111, 112, 113,
            121, 122, 123
        ]
        model.add_ctetra(eid, pid, nids, comment='ctetra10')
        model.setup()
        end_checks(model)

        with self.assertRaises(RuntimeError):
            model.get_length_breakdown(property_ids=None, stop_if_no_length=True)
        with self.assertRaises(RuntimeError):
            model.get_area_breakdown(property_ids=None, stop_if_no_area=True)

        model.get_mass_breakdown(stop_if_no_mass=True) # property_ids=None,
        #model.get_mass_breakdown_detailed(property_ids=None, stop_if_no_mass=True)
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
        model.add_cpyram(eid, pid, nids2, comment='cpyram13')
        model.setup()
        elem2 = model.Element(eid)[0]
        elem2.write(size=8)
        elem2.write(size=16)
        #elem2.write_card_16(is_double=False)

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
        model.setup()
        #model.pop_parse_errors()
        #model.cross_reference()

        eids = None
        nids = None
        p0 = [0., 0., 0.]
        # --------------------------------------------------------------------------------------------
        def sum_loads(loadcase_id, force):
            force1, moment1 = model.sum_forces_moments() # , p0, loadcase_id, include_grav=False
            #force2, moment2 =  model.sum_forces_moments(p0, loadcase_id, eids, nids,
                                                           #include_grav=False, xyz_cid0=None)
            force2 = force1
            #force1, moment1 = sum_forces_moments(model, p0, loadcase_id, include_grav=False, xyz_cid0=None)
            #force2, moment2 =  sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                                           #include_grav=False, xyz_cid0=None)
            assert np.allclose(force1, force), 'sid=%s force1=%s force2=%s' % (loadcase_id, force1, force2)
            #assert np.allclose(force2, force), 'sid=%s force1=%s force2=%s' % (loadcase_id, force1, force2)

        if 0:
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
        model.setup()
        elem2 = model.Element(eid)[0]
        elem2.write(size=8)
        elem2.write(size=16)
        elem2.write_16(is_double=False)
        end_checks(model)
        save_load_deck(model)

    def test_solids_cpentcz(self):
        """tests a CPENTCZ"""
        model = BDF(debug=False)
        cpentcz = model.cpentcz
        psolcz = model.psolcz
        mcid = 42
        thickness = 0.1
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
        psolcz.add(pid, mid, thickness, mcid=mcid)
        model.add_cord2r(mcid, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [21, 22, 23, 24, 25, 26]
        cpentcz.add(eid, pid, nids, comment='cpenta')

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = cpentcz.add(eid, pid, nids2, comment='cpenta15')
        model.setup()
        elem2 = model.Element(eid)[0]
        elem2.write(size=8)
        elem2.write(size=16)
        elem2.write_16(is_double=False)
        end_checks(model)
        save_load_deck(model)

    def test_solids_chexcz(self):
        """tests a CHEXCZ, PSOLCZ"""
        model = BDF(debug=False)
        chexcz = model.chexcz
        psolcz = model.psolcz
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
        psolcz.add(pid, mid, 0.1, mcid=0, comment='')
        model.add_mat1(mid, E, G, nu, rho=1.0)
        nids = [11, 12, 13, 14, 15, 16, 17, 18]
        elem = chexcz.add(eid, pid, nids, comment='chexa')
        model.setup(run_geom_check=False)
        elem = model.Element(eid)[0]
        elem.write(size=8)
        #elem.write_card(size=16)
        #elem.write_card_16(is_double=False)

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = chexcz.add(eid, pid, nids2, comment='chexa16')
        model.setup(run_geom_check=False)
        elem2 = model.Element(eid)[0]
        elem2.write(size=8)
        #elem2.write_card(size=16)
        #elem2.write_card_16(is_double=False)

        eid10 = 1
        pid10 = 10
        #mid10 = 10
        #bulk = 1.
        #c = 1.
        #rho = None
        #model.add_mat10(mid10, bulk, rho, c, ge=0.0, gamma=None,
                        #table_bulk=None, table_rho=None, table_ge=None, table_gamma=None,
                        #comment='')
        thickness = 0.1
        mcid = 42
        model.psolcz.add(pid10, mid, thickness, mcid=mcid, comment='comment')
        #model.add_plsolid(pid10, mid10, stress_strain='GRID', ge=0., comment='')
        elem3 = chexcz.add(eid10, pid10, nids, comment='chexa')
        model.setup()

        chexcz.write(size=8)
        chexcz.write(size=16)
        chexcz.write_16(is_double=False)
        massi = chexcz.mass()

        end_checks(model)
        elem = model.Element(eid)[0]
        assert elem.mass() > 0, elem.mass()
        save_load_deck(model, run_remove_unused=False)

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
        model.setup()
        elem = model.Element(eid)[0]
        elem.write(size=8)
        #elem.write_card(size=16)
        #elem.write_card_16(is_double=False)

        nids2 = copy.deepcopy(nids)
        nids2.append(None)
        eid += 1
        elem2 = model.add_chexa(eid, pid, nids2, comment='chexa16')
        model.setup()
        elem2 = model.Element(eid)[0]
        elem2.write(size=8)
        #elem2.write_card(size=16)
        #elem2.write_card_16(is_double=False)

        eid10 = 1
        pid10 = 10
        mid10 = 10
        bulk = 1.
        c = 1.
        rho = None
        model.add_mat10(mid10, bulk, rho, c, ge=0.0, gamma=None,
                        table_bulk=None, table_rho=None, table_ge=None, table_gamma=None,
                        comment='')
        model.add_psolid(pid10, mid10, fctn='FLUID')
        #model.add_plsolid(pid10, mid10, stress_strain='GRID', ge=0., comment='')
        elem3 = model.add_chexa(eid10, pid10, nids, comment='chexa')
        model.cross_reference()

        chexa = model.chexa
        chexa.write(size=8)
        chexa.write(size=16)
        chexa.write_16(is_double=False)
        chexa.mass()

        end_checks(model)
        elem = model.Element(eid)[0]
        assert elem.mass() > 0, elem.mass()
        save_load_deck(model)

    def test_solids_chexa_pcompls(self):
        """tests a CHEXA8/"""
        model = BDF(debug=False)
        chexa = model.chexa
        pcompls = model.pcompls
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
        nids = [11, 12, 13, 14, 15, 16, 17, 18]
        elem = model.add_chexa(eid, pid, nids, comment='chexa')

        #model.add_psolid(pid, mid)
        model.add_mat1(mid, E, G, nu, rho=1.0)

        #| PCOMPLS |  PID | DIRECT |  CORDM |   SB   |  ANAL  |
        #|         |  C8  |  BEH8  |  INT8  | BEH8H  | INT8H  |
        #|         |  C20 |  BEH20 |  INT20 | BEH20H | INT20H |
        #|         |  ID1 |  MID1  |   T1   | THETA1 |        |
        #|         |  ID2 |  MID2  |   T2   | THETA2 |        |
        #| PCOMPLS | 782  | 1      |        |        |        |
        #|         | 1001 |   171  |   .3   |  12.3  |        |
        #|         | 100  |   175  |   .7   |  77.7  |        |
        card_lines = [
            f'PCOMPLS,{pid},1',
            ',1001,171,0.3,12.23'
        ]
        model.add_card(card_lines, 'PCOMPLS', comment='', ifile=None, is_list=False, has_none=True)
        model.setup()
        mass = chexa.mass()
        rho = chexa.rho()
        assert np.isnan(mass)
        assert np.isnan(rho)
        save_load_deck(model, run_remove_unused=False, run_mass_properties=False, )

    def test_pcomps_slice_card_by_property_id(self):
        """Tests PCOMPS.slice_card_by_property_id with multiple properties."""
        model = BDF(debug=False)
        mid = 1
        E = 3.e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.1)

        # PCOMPS with 2 layers
        model.add_pcomps(
            pid=10,
            global_ply_ids=[101, 102],
            mids=[mid, mid],
            thicknesses=[0.1, 0.2],
            thetas=[0., 45.],
        )
        # PCOMPS with 3 layers
        model.add_pcomps(
            pid=20,
            global_ply_ids=[201, 202, 203],
            mids=[mid, mid, mid],
            thicknesses=[0.05, 0.1, 0.05],
            thetas=[0., 90., 0.],
        )
        model.setup()

        pcomps = model.pcomps
        assert len(pcomps.property_id) == 2

        # Slice single property
        p10 = pcomps.slice_card_by_property_id(np.array([10]))
        assert len(p10.property_id) == 1
        assert p10.property_id[0] == 10
        assert p10.nlayer[0] == 2
        assert len(p10.thickness) == 2
        assert np.allclose(p10.thickness, [0.1, 0.2])

        p20 = pcomps.slice_card_by_property_id(np.array([20]))
        assert len(p20.property_id) == 1
        assert p20.nlayer[0] == 3
        assert len(p20.thickness) == 3
        assert np.allclose(p20.thickness, [0.05, 0.1, 0.05])

        # Slice both
        p_both = pcomps.slice_card_by_property_id(np.array([10, 20]))
        assert len(p_both.property_id) == 2
        assert p_both.nlayer.sum() == 5

    def test_pcompls_slice_card_by_property_id(self):
        """Tests PCOMPLS.slice_card_by_property_id with multiple properties."""
        model = BDF(debug=False)
        mid = 30
        E = 3.e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=1.0)

        # Add two PCOMPLS via card lines (no add_pcompls helper exists)
        card_lines1 = [
            'PCOMPLS,10,1',
            ',1001,30,0.3,12.0',
            ',1002,30,0.5,45.0',
        ]
        model.add_card(card_lines1, 'PCOMPLS', comment='', ifile=None, is_list=False, has_none=True)

        card_lines2 = [
            'PCOMPLS,20,1',
            ',2001,30,0.1,0.0',
        ]
        model.add_card(card_lines2, 'PCOMPLS', comment='', ifile=None, is_list=False, has_none=True)
        model.setup()

        pcompls = model.pcompls
        assert len(pcompls.property_id) == 2

        # Slice single
        p10 = pcompls.slice_card_by_property_id(np.array([10]))
        assert len(p10.property_id) == 1
        assert p10.property_id[0] == 10
        assert p10.nply[0] == 2
        assert len(p10.thickness) == 2
        assert np.allclose(p10.thickness, [0.3, 0.5])

        p20 = pcompls.slice_card_by_property_id(np.array([20]))
        assert len(p20.property_id) == 1
        assert p20.nply[0] == 1
        assert len(p20.thickness) == 1

        # Slice both
        p_both = pcompls.slice_card_by_property_id(np.array([10, 20]))
        assert len(p_both.property_id) == 2
        assert p_both.nply.sum() == 3


    def check_solid(self, model: BDF, eid, etype, pid, ptype, mid, mtype, nsm, rho, volume):
        """checks that various solid methods work"""
        mass = rho * volume
        element = model.Element(eid)[0]
        assert element.nodes.max() > 0
        if ptype == 'PSOLID':
            assert pid in model.psolid.property_id, 'pid is missing for\n%s' % str(element)
        elif ptype == 'PLSOLID':
            assert pid in model.plsolid.property_id, 'pid is missing for\n%s' % str(element)
        else:
            raise NotImplementedError(ptype)
        self.assertEqual(element.type, etype)
        self.assertEqual(element.element_id, eid)
        #self.assertEqual(element.pid_ref.type, ptype)
        self.assertEqual(element.property_id, pid)
        #self.assertEqual(element.pid_ref.mid_ref.type, mtype)
        #self.assertEqual(element.Mid(), mid)
        self.assertEqual(element.volume(), volume)
        self.assertEqual(element.mass(), mass)
        mass_mp = model.mass_sum(element_id=eid)
        #mass_mp = mass_properties(model, element_ids=eid)[0]
        #mass_mp_nsm = mass_properties_nsm(model, element_ids=eid)[0]
        #unused_centroid, unused_xe, unused_ye, unused_ze = element.material_coordinate_system()
        assert np.allclose(mass, mass_mp)
        #assert np.allclose(mass, mass_mp_nsm)


def end_checks(model):
    """various checks"""
    #model.setup()
    model.validate()
    model._verify_bdf(xref=False)
    #model.cross_reference()
    model._verify_bdf(xref=True)
    #model.uncross_reference()
    #model.cross_reference()
    eids, mass, cg, inertia = model.inertia()
    assert mass.sum() > 0, 'mass=%s, cg=%s, inertia=%s' % (mass, cg, inertia)

    #bdf_filename = 'solid_test.bdf'
    #model.write_bdf(bdf_filename)
    #read_bdf(bdf_filename, debug=False)
    #os.remove(bdf_filename)

class TestSolidVolumeHigherOrder(unittest.TestCase):
    """Tests for higher-order solid element volume via Gauss quadrature.

    Verifies:
    - Straight-sided higher-order elements match linear volume (tol=1e-10)
    - Curved elements (midside nodes displaced) produce different volumes
    - Vectorized computation works on batches
    """

    def test_volume_ctetra10_straight(self):
        """CTETRA10 with midside nodes at edge midpoints matches CTETRA4 volume."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import (
            volume_ctetra10, volume_ctetra,
        )
        corners = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], dtype=float)
        mids = np.array([(corners[i]+corners[j])/2
                         for i, j in [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]])
        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]

        vol_quad = volume_ctetra10(all_nodes)[0]
        n1, n2, n3, n4 = [corners[i:i+1] for i in range(4)]
        vol_lin = volume_ctetra(n1, n2, n3, n4)[0]

        # tol=1e-10: quadrature exact for linear geometry
        assert np.isclose(vol_quad, vol_lin, atol=1e-10), \
            f'CTETRA10 straight: {vol_quad} != {vol_lin}'
        assert vol_quad > 0, f'volume should be positive: {vol_quad}'

    def test_volume_ctetra10_curved(self):
        """CTETRA10 with midside nodes pushed outward gives larger volume."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import (
            volume_ctetra10, volume_ctetra,
        )
        corners = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], dtype=float)
        centroid = corners.mean(axis=0)
        edges = [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edges])
        # push midside nodes away from centroid
        for k in range(6):
            d = mids[k] - centroid
            mids[k] += 0.1 * d / np.linalg.norm(d)

        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]
        vol_curved = volume_ctetra10(all_nodes)[0]

        n1, n2, n3, n4 = [corners[i:i+1] for i in range(4)]
        vol_lin = volume_ctetra(n1, n2, n3, n4)[0]

        assert vol_curved > vol_lin, \
            f'curved tet volume {vol_curved} should exceed linear {vol_lin}'

    def test_volume_cpenta15_straight(self):
        """CPENTA15 with midside nodes at edge midpoints matches CPENTA6 volume."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import (
            volume_cpenta15, volume_cpenta,
        )
        corners = np.array([
            [0,0,-1], [1,0,-1], [0,1,-1],
            [0,0,+1], [1,0,+1], [0,1,+1],
        ], dtype=float)
        # Midside: bottom edges (0-1,1-2,2-0), top edges (3-4,4-5,5-3), vertical (0-3,1-4,2-5)
        edge_pairs = [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]

        vol_quad = volume_cpenta15(all_nodes)[0]
        n1,n2,n3,n4,n5,n6 = [corners[i:i+1] for i in range(6)]
        vol_lin = volume_cpenta(n1, n2, n3, n4, n5, n6)[0]

        # tol=1e-10: exact for linear
        assert np.isclose(vol_quad, vol_lin, atol=1e-10), \
            f'CPENTA15 straight: {vol_quad} != {vol_lin}'
        # exact volume = triangle area * height = 0.5 * 2 = 1.0
        assert np.isclose(vol_quad, 1.0, atol=1e-10), f'expected 1.0, got {vol_quad}'

    def test_volume_chexa20_straight(self):
        """CHEXA20 with midside nodes at edge midpoints matches CHEXA8 volume."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import (
            volume_chexa20, volume_chexa,
        )
        corners = np.array([
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1],
        ], dtype=float)
        # 12 edges in Nastran order: 1-2,2-3,3-4,4-1, 5-6,6-7,7-8,8-5, 1-5,2-6,3-7,4-8
        edge_pairs = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]

        vol_quad = volume_chexa20(all_nodes)[0]
        n1,n2,n3,n4,n5,n6,n7,n8 = [corners[i:i+1] for i in range(8)]
        vol_lin = volume_chexa(n1, n2, n3, n4, n5, n6, n7, n8)[0]

        # tol=1e-10: exact for linear
        assert np.isclose(vol_quad, vol_lin, atol=1e-10), \
            f'CHEXA20 straight: {vol_quad} != {vol_lin}'
        assert np.isclose(vol_quad, 1.0, atol=1e-10), f'expected 1.0, got {vol_quad}'

    def test_volume_chexa20_barrel(self):
        """CHEXA20 barrel shape (midside nodes pushed out) gives volume > 1.0."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import volume_chexa20
        corners = np.array([
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1],
        ], dtype=float)
        edge_pairs = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        center = np.array([0.5, 0.5, 0.5])
        for k in range(12):
            d = mids[k] - center
            mids[k] += 0.1 * d / np.linalg.norm(d)
        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]
        vol = volume_chexa20(all_nodes)[0]
        assert vol > 1.0, f'barrel hex volume {vol} should exceed 1.0'

    def test_volume_cpyram13_straight(self):
        """CPYRAM13 with midside nodes at edge midpoints matches CPYRAM5 volume."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import (
            volume_cpyram13, volume_cpyram,
        )
        corners = np.array([
            [-1,-1,0], [1,-1,0], [1,1,0], [-1,1,0], [0,0,1],
        ], dtype=float)
        # Base edges: 0-1, 1-2, 2-3, 3-0; Lateral: 0-4, 1-4, 2-4, 3-4
        edge_pairs = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        all_nodes = np.vstack([corners, mids])[np.newaxis, :, :]

        vol_quad = volume_cpyram13(all_nodes)[0]
        n1,n2,n3,n4,n5 = [corners[i:i+1] for i in range(5)]
        vol_lin = volume_cpyram(n1, n2, n3, n4, n5)[0]

        # tol=1e-6: finite difference derivatives limit precision
        assert np.isclose(vol_quad, vol_lin, atol=1e-6), \
            f'CPYRAM13 straight: {vol_quad} != {vol_lin}'
        # exact: (2*2*1)/3 = 4/3
        assert np.isclose(vol_quad, 4./3., atol=1e-6), f'expected 1.333333, got {vol_quad}'

    def test_volume_batch(self):
        """Vectorized volume works on multiple elements simultaneously."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_volume import volume_chexa20
        corners = np.array([
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1],
        ], dtype=float)
        edge_pairs = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        single = np.vstack([corners, mids])[np.newaxis, :, :]
        batch = np.concatenate([single]*5, axis=0)
        vols = volume_chexa20(batch)
        assert vols.shape == (5,)
        assert np.allclose(vols, 1.0, atol=1e-10)

    def test_volume_ctetra10_via_model(self):
        """CTETRA10 volume() method dispatches to quadrature when midside nodes present."""
        model = BDF(debug=False)
        # Unit tet: nodes at origin, (1,0,0), (0,1,0), (0,0,1)
        corners_xyz = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        for i, xyz in enumerate(corners_xyz, start=1):
            model.add_grid(i, xyz)
        # Midside nodes at edge midpoints
        edges = [(1,2),(2,3),(1,3),(1,4),(2,4),(3,4)]
        mid_nids = []
        for k, (i, j) in enumerate(edges, start=5):
            xyz = (np.array(corners_xyz[i-1]) + np.array(corners_xyz[j-1])) / 2.
            model.add_grid(k, xyz.tolist())
            mid_nids.append(k)
        nids = [1, 2, 3, 4] + mid_nids
        model.add_ctetra(1, 1, nids)
        model.add_psolid(1, 1)
        model.add_mat1(1, 3e7, None, 0.3, rho=1.0)
        model.setup()

        elem = model.ctetra
        vol = elem.volume()[0]
        # tol=1e-10: straight-sided tet, volume = 1/6
        assert np.isclose(vol, 1./6., atol=1e-10), f'expected 1/6, got {vol}'

    def test_volume_chexa20_via_model(self):
        """CHEXA20 volume() method dispatches to quadrature when midside nodes present."""
        model = BDF(debug=False)
        corners_xyz = [
            [0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1],
        ]
        for i, xyz in enumerate(corners_xyz, start=1):
            model.add_grid(i, xyz)
        edge_pairs = [(1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8)]
        mid_nids = []
        for k, (i, j) in enumerate(edge_pairs, start=9):
            xyz = (np.array(corners_xyz[i-1]) + np.array(corners_xyz[j-1])) / 2.
            model.add_grid(k, xyz.tolist())
            mid_nids.append(k)
        nids = list(range(1, 9)) + mid_nids
        model.add_chexa(1, 1, nids)
        model.add_psolid(1, 1)
        model.add_mat1(1, 3e7, None, 0.3, rho=1.0)
        model.setup()

        elem = model.chexa
        vol = elem.volume()[0]
        assert np.isclose(vol, 1.0, atol=1e-10), f'expected 1.0, got {vol}'


class TestConsistentMassMatrix(unittest.TestCase):
    """Tests for consistent mass matrix of solid elements.

    Verifies:
    - Total mass conservation: sum of scalar mass block = rho * volume (tol=1e-10)
    - Symmetry: M = M^T
    - Positive semi-definiteness: all eigenvalues >= 0
    - Higher-order element mass matches linear for straight-sided geometry
    """

    def _check_mass_matrix(self, M: np.ndarray, nnodes: int, expected_mass: float,
                           atol: float = 1e-10):
        """Common checks for a single-element mass matrix."""
        ndof = 3 * nnodes
        assert M.shape == (1, ndof, ndof), f'shape {M.shape} != (1, {ndof}, {ndof})'
        # Symmetry
        assert np.allclose(M[0], M[0].T, atol=1e-14), 'mass matrix not symmetric'
        # Total mass: sum of NxN scalar block (x-direction)
        M_scalar = np.zeros((nnodes, nnodes))
        for i in range(nnodes):
            for j in range(nnodes):
                M_scalar[i, j] = M[0, 3*i, 3*j]
        total_mass = M_scalar.sum()
        assert np.isclose(total_mass, expected_mass, atol=atol), \
            f'total mass {total_mass} != expected {expected_mass}'
        # Positive semi-definite
        eigvals = np.linalg.eigvalsh(M[0])
        assert eigvals.min() >= -1e-14, f'negative eigenvalue: {eigvals.min()}'

    def test_ctetra4_mass_matrix(self):
        """CTETRA4 consistent mass matrix: total mass = rho * V = 2.0 * 1/6."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_mass import consistent_mass_ctetra4
        corners = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
        all_nodes = corners[np.newaxis, :, :]
        rho = np.array([2.0])
        M = consistent_mass_ctetra4(all_nodes, rho)
        self._check_mass_matrix(M, 4, 2.0 / 6.)

    def test_chexa8_mass_matrix(self):
        """CHEXA8 consistent mass matrix: total mass = rho * V = 3.0 * 1.0."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_mass import consistent_mass_chexa8
        corners = np.array([
            [0,0,0],[1,0,0],[1,1,0],[0,1,0],
            [0,0,1],[1,0,1],[1,1,1],[0,1,1],
        ], dtype=float)
        all_nodes = corners[np.newaxis, :, :]
        rho = np.array([3.0])
        M = consistent_mass_chexa8(all_nodes, rho)
        self._check_mass_matrix(M, 8, 3.0)

    def test_cpenta6_mass_matrix(self):
        """CPENTA6 consistent mass matrix: total mass = rho * V = 1.5 * 1.0."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_mass import consistent_mass_cpenta6
        corners = np.array([
            [0,0,-1],[1,0,-1],[0,1,-1],
            [0,0,+1],[1,0,+1],[0,1,+1],
        ], dtype=float)
        all_nodes = corners[np.newaxis, :, :]
        rho = np.array([1.5])
        M = consistent_mass_cpenta6(all_nodes, rho)
        # volume = 0.5 * 2 = 1.0
        self._check_mass_matrix(M, 6, 1.5 * 1.0)

    def test_ctetra10_mass_matches_linear(self):
        """CTETRA10 straight-sided mass matches CTETRA4 total mass."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_mass import (
            consistent_mass_ctetra4, consistent_mass_ctetra10,
        )
        corners = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
        edges = [(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edges])
        rho = np.array([5.0])

        M4 = consistent_mass_ctetra4(corners[np.newaxis], rho)
        M10 = consistent_mass_ctetra10(np.vstack([corners, mids])[np.newaxis], rho)

        # Total mass should be the same
        mass4 = sum(M4[0, 3*i, 3*j] for i in range(4) for j in range(4))
        mass10 = sum(M10[0, 3*i, 3*j] for i in range(10) for j in range(10))
        assert np.isclose(mass4, mass10, atol=1e-10), f'{mass4} != {mass10}'
        self._check_mass_matrix(M10, 10, 5.0 / 6.)

    def test_chexa20_mass_matches_linear(self):
        """CHEXA20 straight-sided mass matches CHEXA8 total mass."""
        from pyNastran.dev.bdf_vectorized3.cards.elements.solid_mass import (
            consistent_mass_chexa8, consistent_mass_chexa20,
        )
        corners = np.array([
            [0,0,0],[1,0,0],[1,1,0],[0,1,0],
            [0,0,1],[1,0,1],[1,1,1],[0,1,1],
        ], dtype=float)
        edge_pairs = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        mids = np.array([(corners[i]+corners[j])/2 for i, j in edge_pairs])
        rho = np.array([7.8])

        M8 = consistent_mass_chexa8(corners[np.newaxis], rho)
        M20 = consistent_mass_chexa20(np.vstack([corners, mids])[np.newaxis], rho)

        mass8 = sum(M8[0, 3*i, 3*j] for i in range(8) for j in range(8))
        mass20 = sum(M20[0, 3*i, 3*j] for i in range(20) for j in range(20))
        # tol=1e-8: 3x3x3 quadrature may have small integration error for linear element
        assert np.isclose(mass8, mass20, atol=1e-8), f'{mass8} != {mass20}'
        self._check_mass_matrix(M20, 20, 7.8, atol=1e-8)

    def test_mass_matrix_via_model(self):
        """consistent_mass_matrix() works through the element class interface."""
        model = BDF(debug=False)
        model.add_grid(1, [0,0,0])
        model.add_grid(2, [1,0,0])
        model.add_grid(3, [0,1,0])
        model.add_grid(4, [0,0,1])
        model.add_ctetra(1, 1, [1,2,3,4])
        model.add_psolid(1, 1)
        model.add_mat1(1, 3e7, None, 0.3, rho=2.7)
        model.setup()

        M = model.ctetra.consistent_mass_matrix()
        expected_mass = 2.7 / 6.
        M_scalar = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                M_scalar[i, j] = M[0, 3*i, 3*j]
        total = M_scalar.sum()
        assert np.isclose(total, expected_mass, atol=1e-10), f'{total} != {expected_mass}'


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
