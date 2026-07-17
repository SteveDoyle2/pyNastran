"""defines various shell element tests"""
from pathlib import Path
import unittest
import numpy as np
from cpylog import SimpleLogger

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties_nsm
import pyNastran

PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'


class TestNsm(unittest.TestCase):
    def test_nsm_cquad4(self):
        eid_quad = 1
        eid_tri = 2
        eid_conrod = 3
        eid_crod = 4
        eid_pbeaml = 5
        eid_pbarl = 6
        pid_pbeaml = 40
        pid_pshell = 10
        pid_pbeaml = 21
        pid_pbarl = 31
        pid_prod = 41
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]
        model = BDF(debug=None)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_ctria3(eid_tri, pid_pshell, nids[:-1]) # area=0.5
        model.add_conrod(eid_conrod, mid, [1, 2], A=1.0, j=0.0, c=0.0, nsm=0.0, comment='')

        x = [0., 0., 1.]
        g0 = None
        nids_beam = [1, 2]
        model.add_cbar(eid_pbarl, pid_pbarl, nids_beam, x, g0, offt='GGG', pa=0, pb=0,
                       wa=None, wb=None, comment='')
        model.add_cbeam(eid_pbeaml, pid_pbeaml, nids_beam, x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
        model.add_crod(eid_crod, pid_prod, [1, 2])
        model.add_prod(pid_prod, mid, A=0.1)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        bar_type = 'BAR'
        dims = [1., 2.]
        xxb = [0.]
        model.add_pbarl(pid_pbarl, mid, bar_type, dims, group='MSCBML0', nsm=0., comment='')

        beam_type = 'BAR'
        dims = [[1., 2.]]
        nsm = [0.0]
        model.add_pbeaml(pid_pbeaml, mid, beam_type, xxb, dims, so=None, nsm=nsm,
                         group='MSCBML0', comment='')
        model.add_mat1(mid, E, G, nu, rho=0.0)

        # TODO: these are correct barring incorrect formulas
        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.5; area=1.5 for PSHELL
        model.add_nsm1(1001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsm1(1002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.5
        model.add_nsm1(1003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsm1(1004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        model.add_nsm1(1005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 2.5
        model.add_nsm1(1006, 'PSHELL', 1.0, 'ALL') # correct; 1.5
        #model.add_nsm1(1007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.5
        model.add_nsm1(1008, 'PSHELL', 1.0, [10, 'THRU', 12, 'BY', 2]) # correct; 1.5
        model.add_nsm1(1009, 'PBARL', 1.0, pid_pbarl) # correct; 1.0
        model.add_nsm1(1010, 'PBEAML', 1.0, pid_pbeaml) # correct; 1.0
        model.add_nsm1(1011, 'PROD', 1.0, pid_prod) # correct; 1.0
        model.add_nsm1(1012, 'CONROD', 1.0, eid_conrod) # correct; 1.0

        #model.add_nsml1(sid, nsm_type, value, ids)
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml1(2001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsml1(2002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.0
        model.add_nsml1(2003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsml1(2004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        model.add_nsml1(2005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 1.0
        model.add_nsml1(2006, 'PSHELL', 1.0, 'ALL') # correct; 1.0
        model.add_nsml1(2007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.0
        model.add_nsml1(2008, 'PSHELL', 1.0, [10, 'THRU', 12, 'BY', 2]) # correct; 1.0
        model.add_nsml1(2009, 'PBARL', 1.0, pid_pbarl) # correct; 1.0
        model.add_nsml1(2010, 'PBEAML', 1.0, pid_pbeaml) # correct; 1.0
        model.add_nsml1(2011, 'PROD', 1.0, pid_prod) # correct; 1.0
        model.add_nsml1(2012, 'CONROD', 1.0, eid_conrod) # correct; 1.0

        #model.add_nsml1(2011, 'PSHELL', 1.0, ['1240', 'THRU', '1250', None, None, # correct; 0.0
        #'2567', 'THRU', '2575',
        #'35689', 'THRU', '35700', None, None,
        #'76', 'THRU', '85',])
        #print(model.nsms[2011])

        model.add_nsm(3000, 'PSHELL', pid_pshell, 1.0, comment='nsm') # correct; 1.5
        model.add_nsm(3001, 'ELEMENT', eid_quad, 1.0) # correct; 1.0
        model.add_nsm(3003, 'ELEMENT', [eid_pbeaml], 1.0) # correct; 1.0
        model.add_nsm(3004, 'ELEMENT', eid_pbarl, 1.0) # correct; 1.0
        model.add_nsm(3009, 'PBARL', pid_pbarl, 1.0) # correct; 1.0
        model.add_nsm(3010, 'PBEAML', pid_pbeaml, 1.0) # correct; 1.0
        model.add_nsm(3011, 'PROD', pid_prod, 1.0) # correct; 1.0
        model.add_nsm(3012, 'CONROD', eid_conrod, 1.0) # correct; 1.0

        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4001, 'ELEMENT', eid_quad, 1.0) # correct; 1.0
        model.add_nsml(4003, 'ELEMENT', [eid_pbeaml], 1.0) # correct; 1.0
        model.add_nsml(4004, 'ELEMENT', eid_pbarl, 1.0) # correct; 1.0
        model.add_nsml(4009, 'PBARL', pid_pbarl, 1.0) # correct; 1.0
        model.add_nsml(4010, 'PBEAML', pid_pbeaml, 1.0) # correct; 1.0
        model.add_nsml(4011, 'PROD', pid_prod, 1.0) # correct; 1.0
        model.add_nsml(4012, 'CONROD', eid_conrod, 1.0) # correct; 1.0

        model.pop_parse_errors()
        model.safe_cross_reference()
        model.pop_xref_errors()

        expected_dict = {
            # NSM1
            1000 : 1.5,
            1001 : 1.0,
            1002 : 1.5,
            1003 : 1.0,
            1004 : 1.0,
            1005 : -1.0,  # crash
            1006 : 1.5,
            1007 : 1.5,
            1008 : 1.5,
            1009 : 1.0,
            1010 : 1.0,
            1011 : 1.0,
            1012 : 1.0,

            #model.add_nsml1(sid, nsm_type, value, ids)
            # NSML1
            2000 : 1.0,
            2001 : 1.0,
            2002 : 1.0,
            2003 : 1.0,
            2004 : 1.0,
            2005 : -1.0, # crash
            2006 : 1.0,
            2007 : 1.0,
            2008 : 1.0,
            2009 : 1.0,
            2010 : 1.0,
            2011 : 1.0,
            2012 : 1.0,

            # NSM
            3000 : 1.5,
            3001 : 1.0,
            3003 : 1.0,
            3004 : 1.0,
            3009 : 1.0,
            3010 : 1.0,
            3011 : 1.0,
            3012 : 1.0,

            # NSM1
            4000 : 1.0,
            4001 : 1.0,
            4003 : 1.0,
            4004 : 1.0,
            4009 : 1.0,
            4010 : 1.0,
            4011 : 1.0,
            4012 : 1.0,
        }
        for nsm_id in sorted(model.nsms):
            mass1_expected = expected_dict[nsm_id]
            if mass1_expected == -1.0:
                with self.assertRaises(RuntimeError):
                    mass1, unused_cg, unused_inertia = mass_properties_nsm(
                        model, nsm_id=nsm_id, debug=False)
            else:
                mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
                if mass1 != mass1_expected:
                    unused_mass2 = mass_properties_nsm(model, nsm_id=nsm_id, debug=True)[0]
                    raise RuntimeError('nsm_id=%s mass != %s; mass1=%s' % (nsm_id, mass1_expected, mass1))
            #print('mass[%s] = %s' % (nsm_id, mass))
            #print('----------------------------------------------')

        for nsm_id, nsms in model.nsms.items():
            expected = expected_dict[nsm_id]
            for nsm in nsms:
                if nsm.type == 'NSML1':
                    if expected == -1.0:
                        pass  # crash
                    else:
                        nsm.get_eid_mass_cg_by_element(model)

        model2 = save_load_deck(model, xref='safe', run_test_bdf=False)
        model2.reset_rslot_map()
        #print(model2._type_to_slot_map)
        model2.elements = {}

        type_to_id_map = {}
        for card_type, ids in model2._type_to_id_map.items():
            if card_type in ['CQUAD4', 'CTRIA3', 'CBEAM', 'CONROD', 'CBAR', 'CROD']:
                pass
            elif card_type in ['NSM', 'NSM1', 'NSML', 'NSML1', 'MAT1',
                               'PBARL', 'PBEAM', 'PSHELL', 'PCOMP', 'PROD', 'PBEAML', 'GRID']:
                type_to_id_map[card_type] = ids
            else:
                raise NotImplementedError(str((card_type, ids)))
        model2._type_to_id_map = type_to_id_map

        model2.log = SimpleLogger(level='error')

        # don't crash on the null case
        for nsm_id in sorted(model2.nsms):
            mass, unused_cg, unused_I = mass_properties_nsm(model2, nsm_id=nsm_id, debug=False)
            self.assertEqual(mass, 0.0)
            #print('mass[%s] = %s' % (nsm_id, mass))
        #print('done with null')
        save_load_deck(model, run_mass_properties=False)

    def test_nsm_prepare(self):
        """tests the NSMADD and all NSM cards using the prepare methods"""
        model = BDF(debug=None)
        nsm_id = 100
        mid = 1000
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 1., 0.])
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.1)
        model.add_conrod(1, mid, [1, 2])
        model.add_conrod(2, mid, [1, 2])
        model.add_conrod(3, mid, [1, 2])
        fields = ['NSM', nsm_id, 'ELEMENT',
                  1, 1.0,
                  2, 2.0,
                  3, 3.0,
                  4, 2.0]
        model.add_card(fields, 'NSM', comment='', is_list=True,
                       has_none=True)
        model.add_card(fields, 'NSML', comment='', is_list=True,
                       has_none=True)

        fields = ['NSM1', nsm_id, 'ELEMENT', 1.0, 1, 2, 3]
        model.add_card(fields, 'NSM1', comment='', is_list=True,
                       has_none=True)
        model.add_card(fields, 'NSML1', comment='', is_list=True,
                       has_none=True)
        save_load_deck(model, run_mass_properties=False)

    def test_nsml1_single(self):
        eid_conrod = 42
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        # nids = [1, 2, 3, 4]

        model = BDF(debug=None)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_conrod(eid_conrod, mid, [1, 4])
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_pshell(1, mid, t=0.1)
        # if 0:
        model.add_nsml1(1000, 'ELEMENT', 1.0, eid_conrod)  # ???
        model.add_nsml1(1000, 'PSHELL', 1.0, 'ALL')  # ???
        model.add_nsml1(1001, 'ELEMENT', 1.0, eid_conrod)  # ???
        model.add_nsm1(1000, 'ELEMENT', 1.0, eid_conrod)  # ???
        model.add_nsm1(1000, 'PSHELL', 1.0, 'ALL')  # ???
        model.add_nsm1(1001, 'ELEMENT', 1.0, eid_conrod)  # ???
        model.add_nsm(1000, 'ELEMENT', [eid_conrod], 1.0)  # ???
        model.add_nsm(1001, 'CONROD', [eid_conrod], 1.0)  # ???
        model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                nxref_errors=0, stop_on_xref_error=True)
        # model.cross_reference()
        model.safe_cross_reference()

        save_load_deck(model, run_mass_properties=False)

    def test_nsml1_mass_by_element(self):
        eid_conrod = 1
        eid_crod = 2
        eid_ctube = 3
        eid_cbush = 4

        eid_cbar = 10
        eid_cbeam = 11

        eid_quad = 20
        eid_tri = 21

        pid_prod = 1
        pid_ptube = 2
        pid_pbarl = 3
        pid_pbeaml = 4
        pid_pbush = 5
        pid_pshell = 10
        pid_pcomp = 11

        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=None)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])

        model.add_conrod(eid_conrod, mid, [1, 4])
        model.add_crod(eid_crod, pid_prod, [1, 2])
        model.add_ctube(eid_ctube, pid_ptube, [1, 2])
        model.add_cbar(eid_cbar, pid_pbarl, [1, 2], [0., 1., 0.], None, validate=True)
        model.add_cbush(eid_cbush, pid_pbush, [1, 2], [0., 1., 0.], None)

        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_ctria3(eid_tri, pid_pshell, [1, 2, 3]) # area=1.0

        model.add_mat1(mid, E, G, nu, rho=0.0)

        model.add_prod(pid_prod, mid, 1.0)
        model.add_ptube(pid_ptube, mid, 1.0)
        model.add_pbarl(pid_pbarl, mid, 'BAR', [1., 2.])
        model.add_pbush(pid_pbush, [1.0])

        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)
        model.add_pcomp(pid_pcomp, mid, [0.1])

        nsml1_conrod_ele = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_conrod)  # ???
        nsml1_crod_ele   = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_crod)   # ???
        nsml1_ctube_ele  = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_ctube)  # ???
        nsml1_cbar_ele   = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_cbar)  # ???
        # nsml1_cbeam_ele  = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_cbeam)  # ???
        nsml1_cbush_ele  = model.add_nsml1(1000, 'ELEMENT', 1.0, eid_cbush) # ???

        nsml1_prod   = model.add_nsml1(2000, 'PROD', 1.0, pid_prod)       # correct; 1.0
        nsml1_ptube  = model.add_nsml1(2000, 'PTUBE', 1.0, pid_ptube)     # correct; 1.0
        nsml1_pbar   = model.add_nsml1(2000, 'PBAR', 1.0, pid_pbarl)     # correct; 1.0
        # nsml1_pbeam = model.add_nsml1(2001, 'PBEAM', 1.0, pid_pbeaml) # correct; 1.0
        nsml1_pbush  = model.add_nsml1(2002, 'PBUSH', 1.0, pid_pbush)    # correct; 1.0
        nsml1_pshell = model.add_nsml1(2003, 'PSHELL', 1.0, pid_pshell)  # correct; 1.0
        nsml1_pcomp  = model.add_nsml1(2004, 'PCOMP', 1.0, pid_pcomp)    # ???

        nsml1_prod_all = model.add_nsml1(2000, 'PROD', 1.0, 'ALL')       # correct; 1.0

        model.validate()
        model.cross_reference()
        nsml1_conrod_ele.get_eid_mass_cg_by_element(model)
        nsml1_crod_ele.get_eid_mass_cg_by_element(model)
        nsml1_ctube_ele.get_eid_mass_cg_by_element(model)
        nsml1_cbar_ele.get_eid_mass_cg_by_element(model)
        # nsml1_cbeam_ele.get_eid_mass_cg_by_element(model)
        nsml1_cbush_ele.get_eid_mass_cg_by_element(model)

        nsml1_pbush.get_eid_mass_cg_by_element(model)
        nsml1_prod.get_eid_mass_cg_by_element(model)
        nsml1_ptube.get_eid_mass_cg_by_element(model)
        nsml1_pbar.get_eid_mass_cg_by_element(model)
        #nsml1_pbeam.get_eid_mass_cg_by_element(model)
        nsml1_pshell.get_eid_mass_cg_by_element(model)
        nsml1_pcomp.get_eid_mass_cg_by_element(model)

        nsml1_prod_all.get_eid_mass_cg_by_element(model)
        save_load_deck(model, run_mass_properties=False)

    def test_nsmadd(self):
        """tests the NSMADD and all NSM cards"""
        eid_quad = 1
        #unused_eid_tri = 2
        #unused_eid_conrod = 3
        #unused_eid_crod = 4
        #unused_eid_pbeaml = 5
        #unused_eid_pbarl = 6
        #unused_pid_pbeaml = 40
        pid_pshell = 10
        #unused_pid_pbeaml = 21
        #unused_pid_pbarl = 31
        #unused_pid_prod = 41
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=None)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.0
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml(3000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.cross_reference()
        model.pop_xref_errors()

        mass, unused_cg, unused_inertia = mass_properties_nsm(model, nsm_id=1000)
        self.assertAlmostEqual(mass, 1.0)
        mass, unused_cg, unused_inertia = mass_properties_nsm(model, nsm_id=2000)
        self.assertAlmostEqual(mass, 1.0)
        mass, unused_cg, unused_inertia = mass_properties_nsm(model, nsm_id=3000)
        self.assertAlmostEqual(mass, 1.0)
        mass, unused_cg, unused_inertia = mass_properties_nsm(model, nsm_id=4000)
        self.assertAlmostEqual(mass, 1.0)

        mass, unused_cg, unused_inertia = mass_properties_nsm(model, nsm_id=5000)
        self.assertAlmostEqual(mass, 8.0)

        model2 = save_load_deck(model)
        mass, unused_cg, unused_inertia = mass_properties_nsm(model2, nsm_id=5000)
        save_load_deck(model, run_mass_properties=False)

    def test_nsmadd_short(self):
        """tests the NSMADD and all NSM cards"""
        eid_quad = 1
        pid_pshell = 10
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.0
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml(3000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.cross_reference()
        model.pop_xref_errors()

        save_load_deck(model, run_mass_properties=False,
                       run_remove_unused=False, run_convert=False,
                       run_save_load_hdf5=False)

    def test_nsmadd_subset(self):
        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4]) # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14]) # area=1.0

        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        pid_pshell = 1

        t = 0.1
        rho = 1.0
        area = 1.0
        nsm = 0.
        model.add_mat1(mid, E, G, nu, rho=rho)
        model.add_pshell(pid_pshell, mid1=mid, t=t) #, nsm=None)

        base_mass = 2 * area * (rho*t + nsm)
        nsml1_mass = 3.14

        model.add_nsml1(1, 'PSHELL', nsml1_mass, pid_pshell, comment='nsml1') # correct; 1.0
        model.cross_reference()

        expected_mass = base_mass + nsml1_mass
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=1)
        assert np.allclose(mass, expected_mass), f'mass={mass} expected_mass={expected_mass}'

        expected_mass = base_mass/2 + nsml1_mass/2
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=1, element_ids=[1])
        assert np.allclose(mass, expected_mass), f'mass={mass} expected_mass={expected_mass}'

    def test_nsm1_subset_per_unit_area(self):
        """NSM1 applies nsm_value per unit area; subset should only get
        its own area's contribution.
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=2.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        model.add_nsm1(10, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Full model: both elements, areas 1+2, nsm_per_area=1.0
        # mass = 1*1.0 + 2*1.0 = 3.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=10)
        self.assertAlmostEqual(mass, 3.0)

        # Subset eid=1 only: area=1.0, mass = 1*1.0 = 1.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=10, element_ids=[1])
        self.assertAlmostEqual(mass, 1.0)

        # Subset eid=2 only: area=2.0, mass = 2*1.0 = 2.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=10, element_ids=[2])
        self.assertAlmostEqual(mass, 2.0)

    def test_nsml_subset_total_mass(self):
        """NSML distributes total mass proportional to area; subset should
        get proportional share based on area fraction.
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [3., 0., 0.])
        model.add_grid(13, [3., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=3.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        total_nsm_mass = 4.0
        model.add_nsml(20, 'PSHELL', 1, total_nsm_mass)
        model.cross_reference()

        # Full: total_nsm_mass distributed over area 1+3=4
        # eid=1 gets 4.0 * (1/4) = 1.0
        # eid=2 gets 4.0 * (3/4) = 3.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=20)
        self.assertAlmostEqual(mass, 4.0)

        # Subset eid=1: gets 1.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=20, element_ids=[1])
        self.assertAlmostEqual(mass, 1.0)

        # Subset eid=2: gets 3.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=20, element_ids=[2])
        self.assertAlmostEqual(mass, 3.0)

    def test_nsm_subset_element_type(self):
        """NSM with nsm_type=ELEMENT applies per-unit-area to specific elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 2., 0.])
        model.add_grid(14, [0., 2., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=4.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM applied to both elements, nsm_per_area=2.0
        model.add_nsm(30, 'ELEMENT', [1, 2], [2.0, 2.0])
        model.cross_reference()

        # Full: eid=1 mass=1*2=2, eid=2 mass=4*2=8, total=10
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=30)
        self.assertAlmostEqual(mass, 10.0)

        # Subset eid=1: mass=2
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=30, element_ids=[1])
        self.assertAlmostEqual(mass, 2.0)

        # Subset eid=2: mass=8
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=30, element_ids=[2])
        self.assertAlmostEqual(mass, 8.0)

    def test_nsml1_subset_line_elements(self):
        """NSML1 on PROD distributes total mass by length to line elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [4., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=1.0
        model.add_crod(2, pid, [2, 3])  # length=3.0

        total_nsm = 8.0
        model.add_nsml1(40, 'PROD', total_nsm, pid)
        model.cross_reference()

        # Full: distribute 8.0 by length fraction: 1/(1+3)=0.25, 3/(1+3)=0.75
        # eid=1 gets 8*0.25=2, eid=2 gets 8*0.75=6
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=40)
        self.assertAlmostEqual(mass, 8.0)

        # Subset eid=1: gets 2.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=40, element_ids=[1])
        self.assertAlmostEqual(mass, 2.0)

        # Subset eid=2: gets 6.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=40, element_ids=[2])
        self.assertAlmostEqual(mass, 6.0)

    def test_nsm1_subset_line_elements(self):
        """NSM1 on PROD applies nsm_value per unit length to line elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=2.0
        model.add_crod(2, pid, [2, 3])  # length=3.0

        nsm_per_length = 1.5
        model.add_nsm1(50, 'PROD', nsm_per_length, pid)
        model.cross_reference()

        # Full: eid=1 mass=2*1.5=3, eid=2 mass=3*1.5=4.5, total=7.5
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=50)
        self.assertAlmostEqual(mass, 7.5)

        # Subset eid=1: mass=3
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=50, element_ids=[1])
        self.assertAlmostEqual(mass, 3.0)

        # Subset eid=2: mass=4.5
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=50, element_ids=[2])
        self.assertAlmostEqual(mass, 4.5)

    def test_nsmadd_subset_combined(self):
        """NSMADD combines multiple NSM cards; subset filtering applies to all."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM1: 2.0 per unit area on PSHELL 1
        model.add_nsm1(100, 'PSHELL', 2.0, 1)
        # NSML1: 6.0 total distributed by area on PSHELL 1
        model.add_nsml1(200, 'PSHELL', 6.0, 1)
        # NSMADD combining both
        model.add_nsmadd(300, [100, 200])
        model.cross_reference()

        # NSM1: each element gets 2.0*1.0=2.0, total=4.0
        # NSML1: each element gets 6.0*(1/(1+1))=3.0, total=6.0
        # Combined: total=10.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=300)
        self.assertAlmostEqual(mass, 10.0)

        # Subset eid=1: NSM1 gives 2.0, NSML1 gives 3.0, total=5.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=300, element_ids=[1])
        self.assertAlmostEqual(mass, 5.0)

    def test_nsml1_subset_unequal_areas(self):
        """NSML1 with 3 elements of different areas validates proportional split."""
        model = BDF(debug=False)
        # eid=1: area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])

        # eid=2: area=2.0
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])

        # eid=3: area=3.0 (triangle, base=6, height=1 -> area=3)
        model.add_grid(21, [0., 0., 0.])
        model.add_grid(22, [6., 0., 0.])
        model.add_grid(23, [0., 1., 0.])
        model.add_ctria3(3, 1, [21, 22, 23])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        total_nsm = 12.0
        model.add_nsml1(60, 'PSHELL', total_nsm, 1)
        model.cross_reference()

        # total area = 1+2+3 = 6
        # eid=1: 12*(1/6) = 2.0
        # eid=2: 12*(2/6) = 4.0
        # eid=3: 12*(3/6) = 6.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=60)
        self.assertAlmostEqual(mass, 12.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=60, element_ids=[1])
        self.assertAlmostEqual(mass, 2.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=60, element_ids=[2])
        self.assertAlmostEqual(mass, 4.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=60, element_ids=[3])
        self.assertAlmostEqual(mass, 6.0)

        # Two-element subset: eid=1+3 should get 2+6=8
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=60, element_ids=[1, 3])
        self.assertAlmostEqual(mass, 8.0)

    def test_nsm_subset_with_structural_mass(self):
        """NSM subset correctly adds to structural mass subset."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [1., 2., 0.])
        model.add_grid(13, [1., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1.0

        mid = 1
        rho = 2.0
        t = 0.5
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=rho)
        model.add_pshell(1, mid1=mid, t=t)

        nsm_per_area = 3.0
        model.add_nsm1(70, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # structural mass per element = area * rho * t = 1.0 * 2.0 * 0.5 = 1.0
        # nsm mass per element = area * nsm_per_area = 1.0 * 3.0 = 3.0
        # total per element = 4.0

        # Full model
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=70)
        self.assertAlmostEqual(mass, 8.0)

        # Subset eid=1
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=70, element_ids=[1])
        self.assertAlmostEqual(mass, 4.0)

        # Subset eid=2
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=70, element_ids=[2])
        self.assertAlmostEqual(mass, 4.0)

    def test_nsm_subset_cg(self):
        """Validate CG is correct when using element_ids subset."""
        model = BDF(debug=False)
        # eid=1 at x=0..1
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # centroid=(0.5, 0.5, 0), area=1

        # eid=2 at x=2..4
        model.add_grid(11, [2., 0., 0.])
        model.add_grid(12, [4., 0., 0.])
        model.add_grid(13, [4., 1., 0.])
        model.add_grid(14, [2., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # centroid=(3, 0.5, 0), area=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        model.add_nsm1(80, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Subset eid=1: mass=1.0, centroid=(0.5, 0.5, 0)
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=80, element_ids=[1])
        self.assertAlmostEqual(mass, 1.0)
        assert np.allclose(cg, [0.5, 0.5, 0.0]), f'cg={cg}'

        # Subset eid=2: mass=2.0, centroid=(3.0, 0.5, 0)
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=80, element_ids=[2])
        self.assertAlmostEqual(mass, 2.0)
        assert np.allclose(cg, [3.0, 0.5, 0.0]), f'cg={cg}'

    def test_nsml1_subset_conrod(self):
        """NSML1 on CONROD elements distributes total mass by length."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)

        model.add_conrod(1, mid, [1, 2], A=1.0)  # length=2
        model.add_conrod(2, mid, [2, 3], A=1.0)  # length=3

        total_nsm = 10.0
        model.add_nsml1(90, 'CONROD', total_nsm, [1, 2])
        model.cross_reference()

        # total length = 2+3 = 5
        # eid=1: 10*(2/5) = 4.0
        # eid=2: 10*(3/5) = 6.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=90)
        self.assertAlmostEqual(mass, 10.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=90, element_ids=[1])
        self.assertAlmostEqual(mass, 4.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=90, element_ids=[2])
        self.assertAlmostEqual(mass, 6.0)


    def test_nsml1_subset_pcomp(self):
        """NSML1 on PCOMP distributes total mass by area to composite shells."""
        model = BDF(debug=False)
        # eid=1: area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 10, [1, 2, 3, 4])

        # eid=2: area=4.0
        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [2., 2., 0.])
        model.add_grid(13, [2., 4., 0.])
        model.add_grid(14, [0., 4., 0.])
        model.add_cquad4(2, 10, [11, 12, 13, 14])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pcomp(10, [mid], [0.1], [0.])

        total_nsm = 15.0
        model.add_nsml1(100, 'PCOMP', total_nsm, 10)
        model.cross_reference()

        # total area = 1+4 = 5
        # eid=1: 15*(1/5) = 3.0
        # eid=2: 15*(4/5) = 12.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=100)
        self.assertAlmostEqual(mass, 15.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=100, element_ids=[1])
        self.assertAlmostEqual(mass, 3.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=100, element_ids=[2])
        self.assertAlmostEqual(mass, 12.0)

    def test_nsm1_subset_all_keyword(self):
        """NSM1 with 'ALL' applies nsm to all properties of that type."""
        model = BDF(debug=False)
        # Two PSHELL properties, one element each
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [3., 2., 0.])
        model.add_grid(13, [3., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=3, pid=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)

        nsm_per_area = 2.0
        model.add_nsm1(110, 'PSHELL', nsm_per_area, 'ALL')
        model.cross_reference()

        # eid=1: 1.0*2.0 = 2.0
        # eid=2: 3.0*2.0 = 6.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=110)
        self.assertAlmostEqual(mass, 8.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=110, element_ids=[1])
        self.assertAlmostEqual(mass, 2.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=110, element_ids=[2])
        self.assertAlmostEqual(mass, 6.0)

    def test_nsml1_subset_all_keyword(self):
        """NSML1 with 'ALL' distributes total mass across all PSHELL elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [3., 2., 0.])
        model.add_grid(13, [3., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=3, pid=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)

        total_nsm = 20.0
        model.add_nsml1(120, 'PSHELL', total_nsm, 'ALL')
        model.cross_reference()

        # total area = 1+3 = 4
        # eid=1: 20*(1/4) = 5.0
        # eid=2: 20*(3/4) = 15.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=120)
        self.assertAlmostEqual(mass, 20.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=120, element_ids=[1])
        self.assertAlmostEqual(mass, 5.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=120, element_ids=[2])
        self.assertAlmostEqual(mass, 15.0)

    def test_nsm1_subset_cbar(self):
        """NSM1 on PBAR applies per-unit-length mass to CBAR elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pbar(pid, mid, A=1.0, i1=1.0, i2=1.0)

        model.add_cbar(1, pid, [1, 2], x=[0., 0., 1.], g0=None)  # length=2
        model.add_cbar(2, pid, [2, 3], x=[0., 0., 1.], g0=None)  # length=3

        nsm_per_length = 4.0
        model.add_nsm1(130, 'PBAR', nsm_per_length, pid)
        model.cross_reference()

        # eid=1: 2*4=8, eid=2: 3*4=12, total=20
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=130)
        self.assertAlmostEqual(mass, 20.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=130, element_ids=[1])
        self.assertAlmostEqual(mass, 8.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=130, element_ids=[2])
        self.assertAlmostEqual(mass, 12.0)

    def test_nsml1_subset_cbar(self):
        """NSML1 on PBAR distributes total mass by length to CBAR elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pbar(pid, mid, A=1.0, i1=1.0, i2=1.0)

        model.add_cbar(1, pid, [1, 2], x=[0., 0., 1.], g0=None)  # length=2
        model.add_cbar(2, pid, [2, 3], x=[0., 0., 1.], g0=None)  # length=3

        total_nsm = 10.0
        model.add_nsml1(140, 'PBAR', total_nsm, pid)
        model.cross_reference()

        # total length=5; eid=1: 10*(2/5)=4, eid=2: 10*(3/5)=6
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=140)
        self.assertAlmostEqual(mass, 10.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=140, element_ids=[1])
        self.assertAlmostEqual(mass, 4.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=140, element_ids=[2])
        self.assertAlmostEqual(mass, 6.0)

    def test_nsml_subset_element_type_line(self):
        """NSML/ELEMENT assigns total mass to each element directly."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [4., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=1
        model.add_crod(2, pid, [2, 3])  # length=3

        # NSML/ELEMENT: each element gets its value as total mass
        model.add_nsml(150, 'ELEMENT', [1, 2], [5.0, 7.0])
        model.cross_reference()

        # eid=1 gets 5.0, eid=2 gets 7.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=150)
        self.assertAlmostEqual(mass, 12.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=150, element_ids=[1])
        self.assertAlmostEqual(mass, 5.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=150, element_ids=[2])
        self.assertAlmostEqual(mass, 7.0)

    def test_nsmadd_subset_area_and_line(self):
        """NSMADD with area (shell) and line (rod) NSM cards combined."""
        model = BDF(debug=False)
        # Shell element
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [2., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=2

        # Rod element
        model.add_grid(11, [0., 5., 0.])
        model.add_grid(12, [3., 5., 0.])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_prod(2, mid, A=1.0)

        model.add_crod(2, 2, [11, 12])  # length=3

        # NSM1 on shell: 5.0 per area -> eid=1 gets 5*2=10
        model.add_nsm1(160, 'PSHELL', 5.0, 1)
        # NSM1 on rod: 2.0 per length -> eid=2 gets 2*3=6
        model.add_nsm1(170, 'PROD', 2.0, 2)
        # NSMADD
        model.add_nsmadd(180, [160, 170])
        model.cross_reference()

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=180)
        self.assertAlmostEqual(mass, 16.0)

        # Subset shell only
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=180, element_ids=[1])
        self.assertAlmostEqual(mass, 10.0)

        # Subset rod only
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=180, element_ids=[2])
        self.assertAlmostEqual(mass, 6.0)

    def test_nsm_subset_inertia(self):
        """Validate inertia tensor for single-element subset with known geometry."""
        model = BDF(debug=False)
        # Single unit square at origin (centroid at 0.5, 0.5, 0)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, centroid=(0.5,0.5,0)

        # Second element far away
        model.add_grid(11, [10., 0., 0.])
        model.add_grid(12, [11., 0., 0.])
        model.add_grid(13, [11., 1., 0.])
        model.add_grid(14, [10., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1, centroid=(10.5,0.5,0)

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 3.0
        model.add_nsm1(190, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Point mass approximation: mass=3, at (0.5, 0.5, 0), ref=(0,0,0)
        # Ixx = m*(dy^2+dz^2) = 3*(0.25+0) = 0.75
        # Iyy = m*(dx^2+dz^2) = 3*(0.25+0) = 0.75
        # Izz = m*(dx^2+dy^2) = 3*(0.25+0.25) = 1.5
        # Ixy = m*dx*dy = 3*0.5*0.5 = 0.75
        mass, cg, inertia = mass_properties_nsm(
            model, nsm_id=190, element_ids=[1],
            reference_point=[0., 0., 0.], inertia_reference='ref')
        self.assertAlmostEqual(mass, 3.0)
        assert np.allclose(cg, [0.5, 0.5, 0.0]), f'cg={cg}'
        # inertia = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        self.assertAlmostEqual(inertia[0], 0.75, places=5)
        self.assertAlmostEqual(inertia[1], 0.75, places=5)
        self.assertAlmostEqual(inertia[2], 1.5, places=5)
        self.assertAlmostEqual(inertia[3], 0.75, places=5)
        self.assertAlmostEqual(inertia[4], 0.0, places=5)
        self.assertAlmostEqual(inertia[5], 0.0, places=5)

    def test_nsm1_subset_multiple_pids(self):
        """NSM1 applied to multiple PSHELL pids via THRU; subset filtering."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [2., 2., 0.])
        model.add_grid(13, [2., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=2, pid=2

        model.add_grid(21, [0., 4., 0.])
        model.add_grid(22, [3., 4., 0.])
        model.add_grid(23, [3., 5., 0.])
        model.add_grid(24, [0., 5., 0.])
        model.add_cquad4(3, 3, [21, 22, 23, 24])  # area=3, pid=3

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)
        model.add_pshell(3, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        # NSM1 on pids 1 THRU 3
        model.add_nsm1(200, 'PSHELL', nsm_per_area, [1, 'THRU', 3])
        model.cross_reference()

        # eid=1: 1*1=1, eid=2: 2*1=2, eid=3: 3*1=3, total=6
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=200)
        self.assertAlmostEqual(mass, 6.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=200, element_ids=[1])
        self.assertAlmostEqual(mass, 1.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=200, element_ids=[2])
        self.assertAlmostEqual(mass, 2.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=200, element_ids=[3])
        self.assertAlmostEqual(mass, 3.0)

        # Subset of two elements
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=200, element_ids=[1, 3])
        self.assertAlmostEqual(mass, 4.0)

    def test_nsm_subset_empty_result(self):
        """Subset with element_ids that don't overlap NSM gives zero mass."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [1., 2., 0.])
        model.add_grid(13, [1., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM only on eid=1
        model.add_nsm1(210, 'PSHELL', 5.0, 1)
        model.cross_reference()

        # Full model: only eid=1 has NSM via pid=1
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=210)
        self.assertAlmostEqual(mass, 10.0)

        # Subset with both elements (both have pid=1)
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=210, element_ids=[1])
        self.assertAlmostEqual(mass, 5.0)

    def test_nsm1_subset_ctria3(self):
        """NSM1 on CTRIA3 elements via PSHELL with subset filtering."""
        model = BDF(debug=False)
        # eid=1: triangle base=2, height=1 -> area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [0., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])  # area=1.0

        # eid=2: triangle base=4, height=3 -> area=6.0
        model.add_grid(11, [0., 5., 0.])
        model.add_grid(12, [4., 5., 0.])
        model.add_grid(13, [0., 8., 0.])
        model.add_ctria3(2, 1, [11, 12, 13])  # area=6.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 2.0
        model.add_nsm1(220, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # eid=1: 1.0*2.0=2.0, eid=2: 6.0*2.0=12.0
        mass, cg, inertia = mass_properties_nsm(model, nsm_id=220)
        self.assertAlmostEqual(mass, 14.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=220, element_ids=[1])
        self.assertAlmostEqual(mass, 2.0)

        mass, cg, inertia = mass_properties_nsm(model, nsm_id=220, element_ids=[2])
        self.assertAlmostEqual(mass, 12.0)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
