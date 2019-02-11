import unittest
import os

import pyNastran
from pyNastran.utils.log import get_logger

#from pyNastran.bdf.bdf import BDF
#from pyNastran.op2.op2 import OP2, FatalError, read_op2
#from pyNastran.op2.op2_interface.op2_common import get_scode_word
from pyNastran.op2.op2_geom import read_op2_geom#, OP2Geom,
from pyNastran.op2.test.test_op2 import run_op2, read_op2
from pyNastran.op2.dev.op2_writer import OP2Writer

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.abspath(os.path.join(PKG_PATH, '..', 'models'))


class TestOP2Writer(unittest.TestCase):
    def test_write_1(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2_filename_debug = os.path.join(folder, 'solid_bending.debug.out')
        op2_filename_out = os.path.join(folder, 'solid_bending_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'solid_bending_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, include_results='displacements')
        #op2 = read_op2(op2_filename, debug_file=op2_filename_debug, include_results='displacements')

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False,
                       delete_objects=True)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        #op2b = read_op2(op2_filename_out, debug_file=op2_filename_debug_out)

    def test_write_2(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'solid_bending')
        op2_filename = os.path.join(folder, 'solid_bending.op2')
        op2_filename_debug = os.path.join(folder, 'solid_bending.debug.out')
        op2_filename_out = os.path.join(folder, 'solid_bending_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'solid_bending_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, include_results='stress')
        #op2 = read_op2(op2_filename, debug_file=op2_filename_debug, include_results=['displacements', 'stress'])

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False,
                       delete_objects=True)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        #op2b = read_op2(op2_filename_out, debug_file=op2_filename_debug_out)


if __name__ == '__main__':
    unittest.main()
