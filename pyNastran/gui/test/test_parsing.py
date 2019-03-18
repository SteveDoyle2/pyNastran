"""defines various GUI unit tests"""
from __future__ import print_function
#import os
import unittest
#import numpy as np

#import pyNastran
from pyNastran.gui.arg_handling import get_inputs


class GuiParsing(unittest.TestCase):
    """tests parsing of the pyNastranGUI command line"""
    def test_parse_1(self):
        """tests parsing of the pyNastranGUI command line"""
        keys_to_remove = ['noupdate', 'log', 'test', 'geomscript', 'postscript', 'qt',
                          'plugin', 'is_groups', 'user_geom', 'user_points', 'debug']
        with open('fem.bdf', 'w') as unused_bdf_file:
            pass
        with open('fem.op2', 'w') as unused_op2_file:
            pass
        with open('fem.tri', 'w') as unused_tri_file:
            pass

        args = ['pyNastranGUI']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': None, 'output': None, 'input': None}, out
        #print(out, '\n')

        args = ['pyNastranGUI', 'fem.bdf']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['nastran'], 'output': [], 'input': ['fem.bdf']}, out
        #print(out, '\n')

        args = ['pyNastranGUI', 'fem.tri']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['cart3d'], 'output': [], 'input': ['fem.tri']}, out
        #print(out, '\n')

        args = ['pyNastranGUI', '-f', 'nastran', 'fem.bdf']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['nastran'], 'output': [], 'input': ['fem.bdf']}, out
        #print(out, '\n')

        args = ['pyNastranGUI', '-f', 'nastran', '-i', 'fem.bdf']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['nastran'], 'output': [], 'input': ['fem.bdf']}, out
        #print(out, '\n')

        args = ['pyNastranGUI', '-f', 'nastran', '-i', 'fem.bdf', '-o', 'fem.op2']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['nastran'], 'output': ['fem.op2'], 'input': ['fem.bdf']}, out
        #print(out, '\n')

        args = ['pyNastranGUI', 'fem.bdf', 'fem.op2']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        assert out == {'format': ['nastran'], 'output': ['fem.op2'], 'input': ['fem.bdf']}, out

def remove_args(dicti, *keys_to_remove):
    """removes keys from a dictionary to declutter the comparison"""
    for key in keys_to_remove:
        if key in dicti:
            del dicti[key]

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
