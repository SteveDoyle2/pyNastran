"""defines various GUI unit tests"""
from __future__ import print_function
import unittest

from six import PY2
#from pyNastran.gui.arg_handling import get_inputs
from pyNastran.bdf.test.test_bdf import run_argparse
if PY2:
    FileNotFoundError = IOError


class TestBDFParsing(unittest.TestCase):
    """tests parsing of the pyNastranGUI command line"""
    def test_parse_1(self):
        """tests parsing of the pyNastranGUI command line"""
        #keys_to_remove = ['noupdate', 'log', 'test', 'geomscript', 'postscript', 'qt',
                          #'plugin', 'is_groups', 'user_geom', 'user_points', 'debug']
        with open('fem.bdf', 'w') as unused_bdf_file:
            pass

        args = ['test_bdf', 'fem.bdf']
        out = run_argparse(argv=args)
        print(out)
        #print(out)

        #args = ['test_bdf', 'fem.bdf', 'fem2.bdf']
        #out = run_argparse(argv=args)
        #assert isinstance(out, string_types), out # error

        args = ['test_bdf', 'fem.bdf', '-x']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--xref']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--safe']
        out = run_argparse(argv=args)

        #args = ['test_bdf', 'fem.bdf', '--safe', '--xref']
        #out = run_argparse(argv=args)
        #assert isinstance(out, string_types), out # error

        args = ['test_bdf', 'fem.bdf', '-xc']
        out = run_argparse(argv=args)
        #print(out)

        args = ['test_bdf', 'fem.bdf', '-xcldL']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--encoding', 'latin1']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '-e', '100']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--dumplines']
        out = run_argparse(argv=args)
        args = ['test_bdf', 'fem.bdf', '--dictsort']
        out = run_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--profile']
        out = run_argparse(argv=args)
        str(out)

        #args = ['test_bdf', '-h']
        #out = run_argparse(argv=args)
        #print(out)


#def remove_args(dicti, *keys_to_remove):
    #"""removes keys from a dictionary to declutter the comparison"""
    #for key in keys_to_remove:
        #if key in dicti:
            #del dicti[key]

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
