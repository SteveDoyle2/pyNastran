import os
import unittest
from numpy import array_equal, allclose
from cpylog import get_logger

import pyNastran
from pyNastran.converters.panair.panair_grid import PanairGrid

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'converters', 'panair')


class TestPanair(unittest.TestCase):

    def test_panair_io_01(self):
        """test the M100 model"""
        log = get_logger(level='warning')
        in_filename = os.path.join(TEST_PATH, 'M100', 'M100.inp')
        out_filename = os.path.join(TEST_PATH, 'M100', 'M100_out.inp')
        #with open(infile_name, 'w') as f:
            #f.write(lines)

        model = PanairGrid(log=log, debug=False)
        model.read_panair(in_filename)
        model.write_panair(out_filename)
        (points, elements, regions, kt, cp_nrom) = model.get_points_elements_regions()

        model.write_panair('junk_m100.inp')
        os.remove('junk_m100.inp')
        model.print_options()
        model.print_abutments()
        model.print_grid_summary()
        model.print_out_header()
        os.remove(out_filename)

    def test_panair_io_02(self):
        model = PanairGrid(log=None, debug=False)
        msg = (
            '$circular sections - nacelle with composite panels\n'
            # =kn\n'
            '1.\n'
            # =kt\n'
            '1.\n'
            # =nopt                                                                 netname\n'
            '0.                                                                    cowlu\n'
            # =nm\n'
            '20.\n'
            # =xs(1)    ri(1)     xs(2)     ri(2)     xs(*)     ri(*)\n'
            '    2.0000    2.3000    1.5756    2.3000    1.1486    2.3000\n'
            '    0.7460    2.3030    0.4069    2.3286    0.1624    2.3790\n'
            '    0.0214    2.4542   -0.0200    2.5485    0.0388    2.6522\n'
            '    0.2056    2.7554    0.4869    2.8522    0.8883    2.9413\n'
            '    1.4250    3.0178    2.1188    3.0656    2.9586    3.0658\n'
            '    3.8551    3.0175    4.6715    2.9439    5.3492    2.8700\n'
            '    6.0000    2.7842    6.4687    2.7442\n'
            # =nn\n'
            '5.\n'
            # =th(1)    th(2)     th(3)     th(4)     th(5)\n'
            '-90.      -45.      0.        45.       90.\n'
        )
        section = msg.split('\n')
        model._read_circular_section(section)
        model.write_panair('junk_circ.inp')
        os.remove('junk_circ.inp')

    def test_panair_io_03(self):
        """tests the SWB model"""
        log = get_logger(level='warning')
        in_filename = os.path.join(TEST_PATH, 'SWB', 'SWB.inp')
        #out_filename = os.path.join(TEST_PATH, 'M100', 'M100_out.inp')

        model = PanairGrid(log=log, debug=False)
        model.read_panair(in_filename)
        (points, elements, regions, kt, cp_nrom) = model.get_points_elements_regions()
        model.write_panair('junk_swb.inp')
        os.remove('junk_swb.inp')

if __name__ == '__main__':  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))
