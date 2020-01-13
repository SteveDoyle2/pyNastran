"""defines various GUI unit tests"""
import os
import unittest

#from pyNastran.gui.arg_handling import get_inputs
from pyNastran.bdf.test.test_bdf import test_bdf_argparse


class TestBDFParsing(unittest.TestCase):
    """tests parsing of the pyNastranGUI command line"""
    def test_parse_1(self):
        """tests parsing of the pyNastranGUI command line"""
        #keys_to_remove = ['noupdate', 'log', 'test', 'geomscript', 'postscript', 'qt',
                          #'plugin', 'is_groups', 'user_geom', 'user_points', 'debug']
        with open('fem.bdf', 'w') as unused_bdf_file:
            pass

        args = ['test_bdf', 'fem.bdf']
        out = test_bdf_argparse(argv=args)
        #print(out)

        #args = ['test_bdf', 'fem.bdf', 'fem2.bdf']
        #out = test_bdf_argparse(argv=args)
        #assert isinstance(out, str), out # error

        args = ['test_bdf', 'fem.bdf', '-x']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--xref']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--safe']
        out = test_bdf_argparse(argv=args)

        #args = ['test_bdf', 'fem.bdf', '--safe', '--xref']
        #out = test_bdf_argparse(argv=args)
        #assert isinstance(out, str), out # error

        args = ['test_bdf', 'fem.bdf', '-xc']
        out = test_bdf_argparse(argv=args)
        #print(out)

        args = ['test_bdf', 'fem.bdf', '-xcdL']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '-xclL']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--encoding', 'latin1']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '-e', '100']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--dumplines']
        out = test_bdf_argparse(argv=args)
        args = ['test_bdf', 'fem.bdf', '--dictsort']
        out = test_bdf_argparse(argv=args)

        args = ['test_bdf', 'fem.bdf', '--profile']
        out = test_bdf_argparse(argv=args)
        str(out)

        #args = ['test_bdf', '-h']
        #out = test_bdf_argparse(argv=args)
        #print(out)
        os.remove('fem.bdf')


#def remove_args(dicti, *keys_to_remove):
    #"""removes keys from a dictionary to declutter the comparison"""
    #for key in keys_to_remove:
        #if key in dicti:
            #del dicti[key]

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
