import unittest
import os
from cpylog import get_logger

import pyNastran
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

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            include_results='displacements')

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        assert op2 == op2b

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

        op2 = read_op2(op2_filename, debug_file=op2_filename_debug,
                       #include_results=['displacements', 'stress']
                       )

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        assert op2 == op2b

    def _test_write_3(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'static_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'static_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'static_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            #include_results='displacements'
                            )

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)

    def test_write_4(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'static_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'static_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'static_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'static_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug)

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    def test_write_5(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'mode_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'mode_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'mode_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            '*_strain_energy',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results,
                            #include_results='eigenvectors',
                            #include_results=['crod_stress', 'cbar_stress'],
                            #include_results=['crod_force', 'cbar_force'],
                            #include_results='element_forces',
                            #include_results='stress',
                            )

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    def test_write_6(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'transient_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'transient_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'transient_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug)

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    def test_write_7(self):
        """tests basic op2 writing"""
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'freq_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'freq_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'freq_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug)

        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    def test_thermal_1(self):
        """tests basic op2 thermal writing"""
        folder = os.path.join(MODEL_PATH, 'elements')
        op2_filename = os.path.join(folder, 'time_thermal_elements.op2')
        op2_filename_debug = os.path.join(folder, 'time_thermal_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'time_thermal_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'time_thermal_elements.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            'chbdye_thermal_load',
            'chexa_thermal_load',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results,
                            #include_results='eigenvectors',
                            #include_results=['crod_stress', 'cbar_stress'],
                            #include_results=['crod_force', 'cbar_force'],
                            #include_results='element_forces',
                            #include_results='stress',
                            )
        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    def test_thermal_2(self):
        """tests basic op2 thermal writing"""
        folder = os.path.join(MODEL_PATH, 'other')
        op2_filename = os.path.join(folder, 'hd15306.op2')
        op2_filename_debug = os.path.join(folder, 'hd15306.debug.out')
        op2_filename_out = os.path.join(folder, 'hd15306_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'hd15306.debug.out')
        #debug_file = 'solid_bending.debug.out'
        model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            'chbdyg_thermal_load',
            'crod_thermal_load',
            'cquad4_thermal_load',
            #'chbdye_thermal_load',
            #'chexa_thermal_load',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, debug=True,
                            #include_results='eigenvectors',
                            #include_results=['crod_stress', 'cbar_stress'],
                            #include_results=['crod_force', 'cbar_force'],
                            #include_results='element_forces',
                            #include_results='stress',
                            )
        print(op2.get_op2_stats(short=True))
        op2w = OP2Writer(op2)
        op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)

    #def test_thermal_3(self):
        #"""tests basic op2 thermal writing"""
        #folder = os.path.join(MODEL_PATH, 'other')
        #op2_filename = os.path.join(folder, 'ofprand1.op2')
        #op2_filename_debug = os.path.join(folder, 'ofprand1.debug.out')
        #op2_filename_out = os.path.join(folder, 'ofprand1_out.op2')
        #op2_filename_debug_out = os.path.join(folder, 'ofprand1.debug.out')
        ##debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        ##debug_file = model + '.debug.out'

        #exclude_results = [
            ##'chbdyg_thermal_load',
            ##'crod_thermal_load',
            ##'cquad4_thermal_load',
            ##'chbdye_thermal_load',
            ##'chexa_thermal_load',
        #]
        #op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            #exclude_results=exclude_results,
                            ##include_results='eigenvectors',
                            ##include_results=['crod_stress', 'cbar_stress'],
                            ##include_results=['crod_force', 'cbar_force'],
                            ##include_results='element_forces',
                            ##include_results='stress',
                            #)
        #print(op2.get_op2_stats(short=True))
        #op2w = OP2Writer(op2)
        #op2w.write_op2(op2_filename_out, obj=op2, is_mag_phase=False)
        #op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        ##op2b = read_op2(op2_filename_out, debug_file=op2_filename_debug_out)
        #op2.assert_op2_equal(op2b,
                             #skip_results=['params', ],
                             #stop_on_failure=True, debug=False)

if __name__ == '__main__':
    unittest.main()
