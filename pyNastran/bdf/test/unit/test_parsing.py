"""defines various GUI unit tests"""
import os
import unittest

#from pyNastran.gui.arg_handling import get_inputs
from pyNastran.bdf.test.test_bdf import parse_args


class TestBDFParsing(unittest.TestCase):
    """tests parsing of the pyNastranGUI command line"""
    def test_parse_1(self):
        """tests parsing of the pyNastranGUI command line"""
        #keys_to_remove = ['noupdate', 'log', 'test', 'geomscript', 'postscript', 'qt',
                          #'plugin', 'is_groups', 'user_geom', 'user_points', 'debug']
        with open('fem.bdf', 'w') as unused_bdf_file:
            pass

        args = ['test_bdf', 'fem.bdf']
        out = parse_args(argv=args)
        #print(out)

        #args = ['test_bdf', 'fem.bdf', 'fem2.bdf']
        #out = parse_args(argv=args)
        #assert isinstance(out, str), out # error

        args = ['test_bdf', 'fem.bdf', '-x']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '--xref']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '--safe']
        out = parse_args(argv=args)

        #args = ['test_bdf', 'fem.bdf', '--safe', '--xref']
        #out = parse_args(argv=args)
        #assert isinstance(out, str), out # error

        args = ['test_bdf', 'fem.bdf', '-xc']
        out = parse_args(argv=args)
        #print(out)

        args = ['test_bdf', 'fem.bdf', '-xcd', '--skip_loads']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '-xcl', '--skip_loads']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '--encoding', 'latin1']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '-e', '100']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '--dumplines']
        out = parse_args(argv=args)
        args = ['test_bdf', 'fem.bdf', '--dictsort']
        out = parse_args(argv=args)

        args = ['test_bdf', 'fem.bdf', '--profile']
        out = parse_args(argv=args)
        str(out)

        #args = ['test_bdf', '-h']
        #out = parse_args(argv=args)
        #print(out)
        os.remove('fem.bdf')


#def remove_args(dicti, *keys_to_remove):
    #"""removes keys from a dictionary to declutter the comparison"""
    #for key in keys_to_remove:
        #if key in dicti:
            #del dicti[key]

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
