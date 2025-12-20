"""tests the pyNastran solver"""
import os
import unittest
from pathlib import Path
from cpylog import SimpleLogger

import pyNastran
from pyNastran.dev.bdf_vectorized.solver.solver import Solver
from pyNastran.utils import print_bad_path


PKG_PATH = Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'dev' / 'bdf_vectorized' / 'solver' / 'test'
log = SimpleLogger('warning', encoding='utf8')

class TestSolverSpring(unittest.TestCase):
    """tests the pyNastran solver"""

    def test_celas1(self):
        """runs a 1 element CELAS1 problem"""
        fargs = {
            '--k' : 1.0, '--f' : 1.0, '--m' : 1.0,
            '--debug' : False,
            'BDFNAME' : os.path.join(TEST_PATH, 'celas1.bdf'),
            'BDFBASE' : 'celas1',
        }
        #C:\NASA\m4\formats\git\pyNastran\pyNastran\dev\bdf_vectorized\solver\test
        assert os.path.exists(fargs['BDFNAME']), print_bad_path(fargs['BDFNAME'])
        solver = Solver(fargs, log=log)
        solver.run_solver()

    def test_celas2(self):
        """runs a 1 element CELAS2 problem"""
        fargs = {
            '--k' : 1.0, '--f' : 1.0, '--m' : 1.0,
            '--debug' : False,
            'BDFNAME' : os.path.join(TEST_PATH, 'celas2.bdf'),
            'BDFBASE' : 'celas2',
        }
        log = SimpleLogger('warning', encoding='utf8')
        solver = Solver(fargs, log=log)
        solver.run_solver()

    #def test_crod(self):
        #"""runs a 1 element CROD problem"""
        #fargs = {
            #'--k' : 1.0, '--f' : 1.0, '--m' : 1.0,
            #'--debug' : False,
            #'BDFNAME' : os.path.join(test_path, 'crod.bdf'),
            #'BDFBASE' : 'crod',
        #}
        #log = SimpleLogger('warning', encoding='utf8')
        #solver = Solver(fargs, log=log)
        #solver.run_solver()

    def test_conrod(self):
        """runs a 1 element CONROD problem"""
        fargs = {
            '--k' : 1.0, '--f' : 1.0, '--m' : 1.0,
            '--debug' : False,
            'BDFNAME' : os.path.join(TEST_PATH, 'conrod.bdf'),
            'BDFBASE' : 'conrod',
        }
        log = SimpleLogger('warning', encoding='utf8')
        solver = Solver(fargs, log=log)
        solver.run_solver()

    def test_cquad4_1(self):
        """runs a 1 element CQUAD4 problem"""
        fargs = {
            '--k' : 1.0, '--f' : 1.0, '--m' : 1.0,
            '--debug' : False,
            'BDFNAME' : os.path.join(TEST_PATH, 'cquad4.bdf'),
            'BDFBASE' : 'cquad4',
        }
        log = SimpleLogger('warning', encoding='utf8')
        solver = Solver(fargs, log=log)
        solver.run_solver()

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
