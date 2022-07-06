import unittest
from pyNastran.bdf.mesh_utils.utils import cmd_line

class TestFlutter(unittest.TestCase):
    """test for bdf flutter"""
    def test_flutter_eas_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0.', 'm']
        cmd_line(argv=args, quiet=True)

        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0', 'ft']
        cmd_line(argv=args, quiet=True)

        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'alt', '0', 'cm']
        with self.assertRaises(AssertionError):
            cmd_line(argv=args, quiet=True)

    def test_flutter_eas_mach(self):
        args = ['bdf', 'flutter', 'english_in', 'eas', '0.1', '100.1', 'knots', '10', 'mach', '0.4', 'none']
        cmd_line(argv=args, quiet=True)

    def test_flutter_mach_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'mach', '0.1', '0.9', '10', 'alt', '0.', 'm']
        cmd_line(argv=args, quiet=True)

    def test_flutter_tas_alt(self):
        args = ['bdf', 'flutter', 'english_in', 'tas', '0.1', '100.1', 'ft/s', '10', 'alt', '0.', 'm']
        cmd_line(argv=args, quiet=True)

    def test_flutter_alt_mach(self):
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_ft', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'si',         'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none']
        cmd_line(argv=args, quiet=True)

    def test_flutter_alt_mach_size(self):
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--size', '8']
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--size', '16']
        cmd_line(argv=args, quiet=True)
        args = ['bdf', 'flutter', 'english_in', 'alt', '0', '10000', 'ft', '10', 'mach', '0.5', 'none', '--clean']
        cmd_line(argv=args, quiet=True)

if __name__ == '__main__':
    unittest.main()
