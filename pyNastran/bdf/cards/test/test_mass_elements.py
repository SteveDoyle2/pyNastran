# coding: utf-8
from __future__ import print_function
import os
import unittest
import numpy as np

from six import iteritems, StringIO

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
#from pyNastran.op2.op2 import OP2, read_op2
#from pyNastran.f06.test.f06_unit_tests import run_model

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
        model.validate()
        read_write(model)

    def test_cmass1(self):
        """tests a CMASS1, PMASS, CMASS2, DDVAL"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        g1 = 1
        c1 = 3

        g2 = 2
        c2 = 4
        cmass1 = model.add_cmass1(eid, pid, g1, c1, g2, c2, comment='cmass1')
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
        cmass2 = model.add_cmass2(eid, mass, g1, c1, g2, c2, comment='cmass2')
        cmass2.write_card(size=8)
        cmass2.write_card(size=16)
        cmass2.write_card(size=16, is_double=True)
        cmass2.raw_fields()

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
        read_write(model)

    def test_mass_3_4(self):
        """tests a CMASS3, PMASS, CMASS4"""
        model = BDF(debug=False)
        eid = 1
        pid = 2
        s1 = 1
        s2 = 2
        cmass3 = model.add_cmass3(eid, pid, s1, s2, comment='cmass3')
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
        cmass4 = model.add_cmass4(eid, mass, s1, s2, comment='cmass4')
        cmass4.write_card(size=8)
        cmass4.write_card(size=16)
        cmass4.write_card(size=16, is_double=True)
        cmass4.raw_fields()

        model.validate()
        model.pop_parse_errors()

        cmass3.write_card(size=8)
        cmass4.write_card(size=8)
        pmass.write_card(size=8)
        read_write(model)

def read_write(model):
    """tests the add_card methods"""
    bdf_file = StringIO()
    model.write_bdf(bdf_file, close=False)
    bdf_file.seek(0)
    model2 = read_bdf(bdf_file, xref=False, punch=True)
    return model2
    #msg = bdf_file.getvalue()
    #bdf_file.close()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

