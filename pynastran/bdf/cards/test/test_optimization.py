# coding: utf-8
from __future__ import print_function
import unittest

from six import iteritems

import os
#import StringIO
#import cStringIO
import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.f06.test.f06_unit_tests import run_model

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')
model_path = os.path.join(pyNastran.__path__[0], '..', 'models')

class TestOpt(unittest.TestCase):
    """
    The cards tested are:
     * DEQATN
    """
    @unittest.expectedFailure
    def test_opt_1(self):
        bdfname = os.path.join(model_path, 'sol200', 'model_200.bdf')
        f06name = None
        op2name = os.path.join(model_path, 'sol200', 'model_200_nx.op2')
        bdf, op2 = run_model(bdfname, op2name,
                             f06_has_weight=False, vectorized=True,
                             encoding='utf-8')

        #op2 = OP2()
        subcase_ids = op2.subcase_key.keys()
        #print('subcase_ids = ', subcase_ids)
        for subcase_id in subcase_ids:
            assert isinstance(subcase_id, int), subcase_id
            for key, dresp in sorted(iteritems(bdf.dresps)):
                dresp.calculate(op2, subcase_id)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
