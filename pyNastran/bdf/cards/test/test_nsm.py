"""defines various shell element tests"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import unittest
from six.moves import StringIO
import numpy as np
from numpy import array

from pyNastran.bdf.bdf import PCOMP, MAT1, BDF
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.mass_properties import _mass_properties_new

class TestNsm(unittest.TestCase):
    def test_nsm_cquad4(self):
        eid_quad = 1
        eid_tri = 2
        eid_conrod = 3
        eid_pbeaml = 4
        eid_pbarl = 5
        pid_pbeaml = 40
        pid_pshell = 10
        pid_pbeaml = 21
        pid_pbarl = 31
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]
        model = BDF()
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
        #model.add_nsml1(sid, nsm_type, value, ids)

        # TODO: these are correct barring incorrect formulas
        model.add_nsml1(1000, 'PSHELL', 1.0, pid_pshell) # correct; 1.0
        model.add_nsml1(1001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsml1(1002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.5
        model.add_nsml1(1003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsml1(1004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        #model.add_nsml1(1005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 2.5
        model.add_nsml1(1006, 'PSHELL', 1.0, 'ALL') # correct; 1.0
        model.add_nsml1(1007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.5
        model.add_nsml1(1008, 'PSHELL', 1.0, [10, 'THRU', 12, 'BY', 2]) # correct; 1.5
        model.add_nsml1(1009, 'PBARL', 1.0, pid_pbarl) # correct; 1.0
        model.add_nsml1(1010, 'PBEAML', 1.0, pid_pbeaml) # correct; 1.0
        #model.add_nsml1(1011, 'PSHELL', 1.0, ['1240', 'THRU', '1250', None, None, # correct; 0.0
                                              #'2567', 'THRU', '2575',
                                              #'35689', 'THRU', '35700', None, None,
                                              #'76', 'THRU', '85',])
        #print(model.nsms[1011])
        model.pop_parse_errors()
        model.cross_reference()
        model.pop_xref_errors()

        for nsm_id in sorted(model.nsms):
            mass, cg, I = _mass_properties_new(model, nsm_id=nsm_id, dev=True)
            print('mass[%s] = %s' % (nsm_id, mass))
            #print('----------------------------------------------')

        model2 = save_load_deck(model)
        #for nsm_id in model2.nsms:  # TODO: totally wrong
            #mass, cg, I = _mass_properties_new(model2, nsm_id=nsm_id, dev=True)
            #print('mass2[%s] = %s' % (nsm_id, mass))

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
