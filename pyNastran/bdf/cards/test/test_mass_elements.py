# coding: utf-8
from __future__ import print_function
import os
import unittest
import numpy as np

from six import StringIO

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck
#from pyNastran.op2.op2 import OP2, read_op2

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')


class TestMassElements(unittest.TestCase):
    """
    The cards tested are:
     * CMASS1
    """
    def test_conm1(self):
        """tests a CONM1"""
        model = BDF(debug=False)
        mass_matrix = np.zeros((6,6), dtype='float32')

        nid = 10
        eid = 10
        mass = 42.
        mass_matrix[0,0] = mass_matrix[1,1] = mass_matrix[2,2] = mass

        conm1 = model.add_conm1(eid, nid, mass_matrix, cid=0, comment='conm1')
        conm1.write_card(size=8)
        conm1.write_card(size=16)
        conm1.write_card(size=16, is_double=True)
        conm1.raw_fields()

        model.add_grid(nid, [0., 0., 0.])
        model.validate()
        save_load_deck(model, run_convert=False)

    def test_conm2(self):
        """tests a conm2"""
        model = BDF(debug=False)
        nid = 10
        eid = 20
        massi = 42.
        model.add_conm2(eid, nid, massi, cid=0, X=None, I=None, comment='conm2')
        model.add_grid(nid, [0., 0., 0.])
        model.validate()
        model.cross_reference()
        model = save_load_deck(model, run_convert=False)
        pids_to_mass, mass_type_to_mass = model.get_mass_breakdown(property_ids=None)
        assert len(pids_to_mass) == 0, pids_to_mass
        assert mass_type_to_mass['CONM2'] == 42., mass_type_to_mass

        mass, cg, I = model.mass_properties(element_ids=None, mass_ids=None, reference_point=None,
                                            sym_axis=None, scale=None)
        assert np.allclose(mass, massi), 'massi=%s mass=%s' % (massi, mass)
        assert np.array_equal(cg, np.zeros(3))
        assert np.array_equal(I, np.zeros(6))

        mass, cg, I = model.mass_properties(element_ids=None, mass_ids=[20], reference_point=None,
                                            sym_axis=None, scale=None)
        assert np.allclose(mass, massi), 'massi=%s mass=%s' % (massi, mass)
        assert np.array_equal(cg, np.zeros(3))
        assert np.array_equal(I, np.zeros(6))

        mass, cg, I = model.mass_properties(element_ids=None, mass_ids=[42], reference_point=None,
                                            sym_axis=None, scale=None)
        assert np.allclose(mass, 0.), 'massi=%s mass=%s' % (massi, mass) ## TODO: is this reasonable behavior

    def test_cmass1(self):
        """tests a CMASS1, PMASS, CMASS2, DDVAL"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        g1 = 1
        c1 = 3

        g2 = 2
        c2 = 4
        cmass1 = model.add_cmass1(eid, pid, [g1, g2], c1, c2, comment='cmass1')
        cmass1.write_card(size=8)
        cmass1.write_card(size=16)
        cmass1.write_card(size=16, is_double=True)
        cmass1.raw_fields()

        mass = 142.
        pmass = model.add_pmass(pid, mass, comment='pmass')
        pmass.write_card(size=8)
        pmass.write_card(size=16)
        pmass.write_card(size=16, is_double=True)
        pmass.raw_fields()

        eid = 10
        cmass2 = model.add_cmass2(eid, mass, [g1, g2], c1, c2, comment='cmass2')
        cmass2.write_card(size=8)
        cmass2.write_card(size=16)
        cmass2.write_card(size=16, is_double=True)
        cmass2.raw_fields()

        model.add_grid(g1, [0., 0., 0.])
        model.add_grid(g2, [0., 0., 0.])

        oid = 3
        ddvals = 1. # promoted to a list
        ddval = model.add_ddval(oid, ddvals, comment='ddval')
        ddval.write_card(size=8)
        ddval.write_card(size=16)
        ddval.write_card(size=16, is_double=True)
        ddval.raw_fields()
        model.validate()
        model.cross_reference()

        cmass1.write_card(size=8)
        cmass2.write_card(size=8)
        pmass.write_card(size=8)
        ddval.write_card(size=8)
        save_load_deck(model, run_convert=False)

    def test_mass_3_4(self):
        """tests a CMASS3, PMASS, CMASS4"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        s1 = 1
        s2 = 2
        cmass3 = model.add_cmass3(eid, pid, [s1, s2], comment='cmass3')
        cmass3.write_card(size=8)
        cmass3.write_card(size=16)
        cmass3.write_card(size=16, is_double=True)
        cmass3.raw_fields()

        mass = 142.
        pmass = model.add_pmass(pid, mass, comment='pmass')
        pmass.write_card(size=8)
        pmass.write_card(size=16)
        pmass.write_card(size=16, is_double=True)
        pmass.raw_fields()

        eid = 10
        cmass4 = model.add_cmass4(eid, mass, [s1, s2], comment='cmass4')
        cmass4.write_card(size=8)
        cmass4.write_card(size=16)
        cmass4.write_card(size=16, is_double=True)
        cmass4.raw_fields()

        model.validate()
        model.pop_parse_errors()

        cmass3.write_card(size=8)
        cmass4.write_card(size=8)
        pmass.write_card(size=8)
        #save_load_deck(model)

    def test_nsm(self):
        """tests the NSM card"""
        model = BDF(debug=False)
        sid = 1
        nsmi = 0.1
        pid = 10
        mid = 100
        model.add_nsm(sid, 'PSHELL', pid, nsmi, comment='nsm')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_ctria3(1, pid, [1, 2, 3]) # A=0.5
        model.add_cquad4(2, pid, [1, 2, 3, 4]) # A=1.0
        model.add_pshell(pid, mid, t=0.1)
        model.add_mat1(mid, 3.0e7, None, 0.3)
        save_load_deck(model, run_convert=False)

    def test_nsm1(self):
        """tests the NSM1 card"""
        model = BDF(debug=False)
        sid = 1
        Type = 'PSHELL'
        nsmi = 0.1
        pid = 10
        mid = 100
        model.add_nsm1(sid, 'PSHELL', nsmi, [pid])
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_ctria3(1, pid, [1, 2, 3]) # A=0.5
        model.add_cquad4(2, pid, [1, 2, 3, 4]) # A=1.0
        model.add_pshell(pid, mid, t=0.1)
        model.add_mat1(mid, 3.0e7, None, 0.3)
        save_load_deck(model, run_convert=False)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()

