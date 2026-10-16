import unittest
from pyNastran.dev.tools.inrel.test_inrel import TestInrel
from pyNastran.dev.tools.test_huth import TestHuth
from pyNastran.dev.tools.pressure_map.test_pressure_map import (
    TestPressureMap,
    TestPanelModelBalance, TestPanelModelLoads,
    TestAeroCombinedArrays, TestAeroXyzScale,
    TestFullModelForce, TestFullModelBalance)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
