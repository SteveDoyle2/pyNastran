# coding: utf-8
from __future__ import print_function
import os
import unittest

from six import iteritems, StringIO

import pyNastran
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.op2.op2 import OP2, read_op2
from pyNastran.f06.test.f06_unit_tests import run_model

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')


class TestMassElements(unittest.TestCase):
    """
    The cards tested are:
     * CMASS1
    """
    def test_conm1(self):
        model = BDF(debug=False)
        mass_matrix = np.zeros((6,6), dtype='float32')
        mass_matrix[0,0] = mass_matrix[1,1] = mass_matrix[2,2] = mass
        nid = 42
        eid = 10
        conm1 = model.add_conm1(eid, nid, mass_matrix, cid=0, comment='conm1')
        conm1.write_card(size=8)
        conm1.write_card(size=16)
        conm1.write_card(size=16, is_double=True)
        conm1.raw_fields()

        model.validate()

    def test_cmass1(self):
        """tests a CMASS1, DDVAL"""
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

        oid = 3
        ddvals = [1]
        ddval = model.add_ddval(oid, ddvals, comment='ddval')
        ddval.write_card(size=8)
        ddval.write_card(size=16)
        ddval.write_card(size=16, is_double=True)
        ddval.raw_fields()
        model.validate()
        model.cross_reference()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()

