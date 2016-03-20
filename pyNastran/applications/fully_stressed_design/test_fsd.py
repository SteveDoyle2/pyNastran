import os
import pyNastran
test_path = pyNastran.__path__[0]
from numpy import dot, array_equal
import unittest

from pyNastran.bdf.bdf import read_bdf, DRESP1
from pyNastran.applications.fully_stressed_design.fully_stressed_design import fully_stressed_design


class TestFSD(unittest.TestCase):
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_fsd_01(self):
        bdf_filename = os.path.abspath(os.path.join(test_path, '..', 'models', 'sol_101_elements',
                                                    'static_solid_shell_bar.bdf'))
        model = read_bdf(bdf_filename)

        oid = 1
        label = 'STRESS'
        response_type = 'STRESS'
        property_type = 'PSHELL'
        region = 'STRESS'
        atta = None
        attb = None
        atti = []
        model.dresps[dresp.oid] = DRESP1(oid, label, response_type,
                                         property_type,
                                         region,
                                         atta, attb, atti)
        keywords = {
            'scr' : 'yes',
            'bat' : 'no',
            'old' : 'no',
        }
        max_stress = 1.033589E+04
        target_stress = 5.0E+3
        regions = {
            4 : [0.0001, 0.7, -target_stress, target_stress],
        }
        regions2 = fully_stressed_design(model, keywords=keywords)
        self.assertTrue(np.array_equal(i, res), 'A i=%s res=%s' % (i, res))

    def atest_fsd_02(self):
        bdf_filename = os.path.abspath(os.path.join(test_path, '..', 'models', 'bwb',
                                                    'BWB_saero.bdf'))
        keywords = {
            'scr' : 'yes',
            'bat' : 'no',
            'old' : 'no',
            'parallel'  : '4',
        }
        #max_stress = 1.033589E+04
        target_stress = 30000.0 # 30
        regions = {
            4 : [0.0001, 0.7, -target_stress, target_stress],
        }
        regions2 = fully_stressed_design(bdf_filename, keywords=keywords)
        self.assertTrue(np.array_equal(i, res), 'A i=%s res=%s' % (i, res))


if __name__ == '__main__':
    unittest.main()
