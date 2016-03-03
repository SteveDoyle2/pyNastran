import pyNastran
test_path = pyNastran.__path__[0]
from numpy import dot, array_equal

from pyNastran.applications.fully_stressed_design.fully_stressed_design import fully_stressed_design


class TestFSD(Tester):
    #def _spike(self):
        #op2 = OP2()
        #op2.set_results('solidStress.oxx')
        #op2.read_op2(op2_filename, vectorized=False)

    def test_fsd_01(self):
        fully_stressed_design
        a = np.array([1., 2., 0.1])
        i = filter1d(a, zero_tol=0.5)
        res = np.array([0, 1])
        self.assertTrue(np.array_equal(i, res), 'A i=%s res=%s' % (i, res))