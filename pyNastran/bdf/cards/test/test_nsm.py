"""defines various shell element tests"""
import os
import unittest
from cpylog import SimpleLogger

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties_nsm
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


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
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_ctria3(eid_tri, pid_pshell, nids[:-1]) # area=0.5
        model.add_conrod(eid_conrod, mid, [1, 2], A=1.0, j=0.0, c=0.0, nsm=0.0, comment='')

        x = [0., 0., 1.]
        g0 = None
        model.add_cbar(eid_pbarl, pid_pbarl, [1, 2], x, g0, offt='GGG', pa=0, pb=0,
                       wa=None, wb=None, comment='')
        model.add_cbeam(eid_pbeaml, pid_pbeaml, nids, x, g0, offt='GGG', bit=None,
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
        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.5
        model.add_nsm1(1001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsm1(1002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.5
        model.add_nsm1(1003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsm1(1004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        model.add_nsm1(1005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 2.5
        model.add_nsm1(1006, 'PSHELL', 1.0, 'ALL') # correct; 1.5
        model.add_nsm1(1007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.5
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
        model.cross_reference()
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
                    mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
            else:
                mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
                if mass1 != mass1_expected:
                    unused_mass2 = mass_properties_nsm(model, nsm_id=nsm_id, debug=True)[0]
                    raise RuntimeError('nsm_id=%s mass != %s; mass1=%s' % (nsm_id, mass1_expected, mass1))
            #print('mass[%s] = %s' % (nsm_id, mass))
            #print('----------------------------------------------')

        model2 = save_load_deck(model, run_test_bdf=False)
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

    def test_nsm_prepare(self):
        """tests the NSMADD and all NSM cards using the prepare methods"""
        model = BDF()
        nsm_id = 100
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


    def test_nsmadd(self):
        """tests the NSMADD and all NSM cards"""
        eid_quad = 1
        unused_eid_tri = 2
        unused_eid_conrod = 3
        unused_eid_crod = 4
        unused_eid_pbeaml = 5
        unused_eid_pbarl = 6
        unused_pid_pbeaml = 40
        pid_pshell = 10
        unused_pid_pbeaml = 21
        unused_pid_pbarl = 31
        unused_pid_prod = 41
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=False)
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

        mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=5000)
        self.assertAlmostEqual(mass, 8.0)
        model2 = save_load_deck(model)
        mass, unused_cg, unused_I = mass_properties_nsm(model2, nsm_id=5000)

    #def test_nsm(self):
        #"""tests a complete nsm example"""
        #bdf_filename = os.path.join(MODEL_PATH, 'nsm', 'nsm.bdf')
        #bdf_filename = os.path.join(MODEL_PATH, 'nsm', 'TEST_NSM_SOL101.bdf')
        #model = read_bdf(bdf_filename)
        #print('    %6s %-9s %s' % ('nsm_id', 'mass', 'nsm'))
        #mass0 = mass_properties_nsm(model, debug=False)[0]
        #for nsm_id in sorted(chain(model.nsms, model.nsmadds)):
            #mass, cg, I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
            #print('    %-6s %-9.4g %.4g' % (nsm_id, mass, mass-mass0))

        #area_breakdown = model.get_area_breakdown()
        #for pid in [20000, 20010]:
            #print('pid=%s area=%.3f' % (pid, area_breakdown[pid]))


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
