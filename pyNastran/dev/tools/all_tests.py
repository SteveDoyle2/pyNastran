import unittest
from pyNastran.dev.tools.inrel.test_inrel import TestInrel
from pyNastran.dev.tools.test_huth import TestHuth
from pyNastran.dev.tools.pressure_map.test_pressure_map import (
    TestPressureMap,
    TestPanelModelBalance, TestPanelModelLoads,
    TestAeroCombinedArrays, TestAeroXyzScale,
    TestFullModelForce, TestFullModelBalance,
    TestLCP,)
from pyNastran.dev.tools.pressure_map.test_structural_lcp import (
    TestConsistentPressure, TestNodalLoads, TestStructuralLCP,
)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
