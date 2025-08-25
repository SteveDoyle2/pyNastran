import os
import unittest
from pyNastran.bdf.mesh_utils.utils import cmd_line, cmd_line_create_flutter
TEST_DIR = os.path.dirname(__file__)


class TestFlutter(unittest.TestCase):
    """test for bdf flutter"""
    def test_bdf_flutter(self):
        """tests a flutter sweep"""
        #UNITS eas EAS1 EAS2 SWEEP_UNIT N CONST_TYPE CONST_VAL
        # [-o OUT_BDF_FILENAME] [--size SIZE | --clean]
        bdf_filename_out = os.path.join(TEST_DIR, 'test_flutter.bdf')

        #bdf flutter english_in mach .05 0.5       101 alt 2500
        args = [
            'bdf', 'flutter', 'english_in',
            'mach', '0.05', '0.5', '21',
            'alt', '2500', 'ft',
        ]
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(args, quiet=True)
        args = [
            'bdf', 'flutter', 'english_in',
            'alt', -10_000, 100_000, 'ft', 41,
            'mach', 0.8, 'NA',
            '-o', bdf_filename_out,
            '--clean',
            '--eas_limit', 1000, 'knots',
        ]
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(args, quiet=True)

        args = [
            'bdf', 'flutter', 'si_mm',
            'tas', 50, 1000, 'ft/s', 11,
            'alt', 0, 'm',
            '-o', bdf_filename_out,
            '--clean',
            '--eas_limit', 1000, 'knots',
        ]
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(args, quiet=True)

        args = [
            'bdf', 'flutter', 'si_mm',
            'eas', 50, 1000, 'ft/s', 11,
            'mach', 0.8, 'none',
            '-o', bdf_filename_out,
        ]
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(args, quiet=True)

    def test_flutter_sweep_eas_const_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0.', 'm']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0', 'ft']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0', 'cm']
        with self.assertRaises(AssertionError):
            cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_eas_const_mach(self):
        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'mach', '0.4', 'none']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_mach_const_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'mach', '0.1', '0.9', '10', 'alt', '0.', 'm']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_tas_const_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'tas', '0.1', '100.1', 'ft/s', '10', 'alt', '0.', 'm']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_alt_const_tas(self):
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '40000', 'ft', '10', 'tas', '100.', 'm/s']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_alt_const_mach(self):
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_ft', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'si',         'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)

    def test_flutter_sweep_alt_const_mach_size(self):
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--size', '8']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--size', '16']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--clean']
        cmd_line_create_flutter(args, quiet=True)
        cmd_line(argv=args, quiet=True)


if __name__ == '__main__':
    unittest.main()
