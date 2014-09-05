import os
import unittest
import pkg_resources

from distutils.spawn import find_executable

from openmdao.main.api import SimulationRoot
from nastranwrapper.test.nastran_models.bar3_static_nastran import Bar3Static

ORIG_DIR = os.getcwd()
DIRECTORY = pkg_resources.resource_filename('nastranwrapper', 'test')

class TestBar3Truss(unittest.TestCase):

    def setUp(self):
        SimulationRoot.chroot(DIRECTORY)

    def tearDown(self):
        SimulationRoot.chroot(ORIG_DIR)

    def test_one_iteration(self):

        static = Bar3Static()
        static.delete_tmp_files = False
        static.stdout = os.devnull
        static.stderr = os.devnull
        static.nastran_filename = "bdf_files/bar3_op2.bdf"
        static.nastran_command = "python"
        static.nastran_command_args = ["fake_nastran.py", 
                                       "bdf_files/test_bar3truss_correct_input.bdf",
                                       "test_bar3truss_correct_output.out",
                                       "test_bar3truss_correct_output.op2",
                                       ]

        # set some variables.
        static.bar1_area = 18.88
        static.bar2_area = 45
        static.bar3_area = 2902.55950333333333
        static.load_x_dir = 50000
        static.load_y_dir = 100000
        static.loadmag = 20
        static.Youngs_Modulus = 70000000
        static.weight_density = 0.289

        static.run()

        # these values were gotten by running it with real nastran
        self.assertAlmostEqual(static.bar1_stress, 2098.76171875 )
        self.assertAlmostEqual(static.bar2_stress, 2088.0517578125 )
        self.assertAlmostEqual(static.bar3_stress, 10.709878921508789 )
        self.assertAlmostEqual(static.displacement_x_dir, 0.0070315715856850147 )
        self.assertAlmostEqual(static.displacement_y_dir, -0.0069601726718246937 )
        self.assertAlmostEqual(static.weight, 118613.7265625 )


if __name__ == "__main__":
    unittest.main()
