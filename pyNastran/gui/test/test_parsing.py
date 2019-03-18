"""
defines various GUI unit tests
"""
from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.gui.arg_handling import get_inputs

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class GuiParsing(unittest.TestCase):
    def test_parse_1(self):
        """tests ``check_for_newer_version``"""
        keys_to_remove = ['noupdate', 'log', 'test', 'geomscript', 'postscript', 'qt', 'plugin', 'is_groups',
                          'user_geom', 'user_points', 'debug']
        with open('fem.bdf', 'w') as bdf_file:
            pass
        with open('fem.tri', 'w') as tri_file:
            pass

        #args = ['pyNastranGUI']
        #args = []
        #out = get_inputs(print_inputs=False, argv=args)
        #remove_args(out, *keys_to_remove)
        #print(out, '\n')

        args = ['C:\\Anaconda2\\Scripts\\pyNastranGUI', 'fem.bdf']
        #args = ['pyNastranGUI', 'fem.bdf']
        #args = ['fem.bdf']
        out = get_inputs(print_inputs=False, argv=args)
        remove_args(out, *keys_to_remove)
        print(out, '\n')

        #args = ['pyNastranGUI', 'fem.tri']
        #args = ['fem.tri']
        #out = get_inputs(print_inputs=False, argv=args)
        #remove_args(out, *keys_to_remove)
        #print(out, '\n')

        #args = ['pyNastranGUI', '-f', 'nastran', 'fem.bdf']
        #args = ['-f', 'nastran', 'fem.bdf']
        #out = get_inputs(print_inputs=False, argv=args)
        #remove_args(out, *keys_to_remove)
        #print(out, '\n')

        #args = ['pyNastranGUI', '-f', 'nastran', '-i', 'fem.bdf']
        #args = ['-f', 'nastran', '-i', 'fem.bdf']
        #out = get_inputs(print_inputs=False, argv=args)
        #remove_args(out, *keys_to_remove)
        #print(out, '\n')

        #args = ['pyNastranGUI', '-f', 'nastran', '-i', 'fem.bdf', '-o', 'fem.op2']
        #args = ['-f', 'nastran', '-i', 'fem.bdf', '-o', 'fem.op2']
        #out = get_inputs(print_inputs=False, argv=args)
        #remove_args(out, *keys_to_remove)
        #print(out, '\n')

def remove_args(dicti, *keys_to_remove):
    for key in keys_to_remove:
        if key in dicti:
            del dicti[key]

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
