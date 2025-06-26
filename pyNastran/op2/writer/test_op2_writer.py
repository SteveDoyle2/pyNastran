import unittest
import os
from pathlib import Path
from cpylog import SimpleLogger

import pyNastran
#from pyNastran.bdf.bdf import BDF
#from pyNastran.op2.op2 import FatalError
#from pyNastran.op2.op2_interface.op2_common import get_scode_word
from pyNastran.op2.op2_geom import read_op2_geom
from pyNastran.op2.op2 import OP2
#from pyNastran.op2.test.test_op2 import run_op2
#from pyNastran.op2.writer.op2_writer import OP2Writer

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = Path(os.path.abspath(os.path.join(PKG_PATH, '..', 'models')))


class TestOP2Writer(unittest.TestCase):
    def test_write_1(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = MODEL_PATH / 'solid_bending'
        op2_filename = folder / 'solid_bending.op2'
        op2_filename_debug = folder / 'solid_bending.debug.out'
        op2_filename_out = folder / 'solid_bending_out.op2'
        op2_filename_debug_out = folder / 'solid_bending_out.debug.out'
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            include_results='displacements', log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        assert op2 == op2b

    def test_write_2(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = MODEL_PATH / 'solid_bending'
        op2_filename = folder / 'solid_bending.op2'
        op2_filename_debug = folder / 'solid_bending.debug.out'
        op2_filename_out = folder / 'solid_bending_out.op2'
        op2_filename_debug_out = folder / 'solid_bending_out.debug.out'
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = OP2(debug=True, log=log, debug_file=op2_filename_debug, mode=None)
        op2.read_op2(op2_filename)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        assert op2 == op2b
        os.remove(op2_filename_debug_out)

    def test_write_3(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = MODEL_PATH / 'sol_101_elements'
        op2_filename = folder / 'static_solid_shell_bar.op2'
        op2_filename_debug = folder / 'static_solid_shell_bar.debug.out'
        op2_filename_out = folder / 'static_solid_shell_bar_out.op2'
        op2_filename_debug_out = folder / 'static_solid_shell_bar_out.debug.out'
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        unused_op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        os.remove(op2_filename_debug_out)

    def test_write_4(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = MODEL_PATH / 'sol_101_elements'
        op2_filename = folder / 'static_solid_shell_bar.op2'
        op2_filename_debug = folder / 'static_solid_shell_bar.debug.out'
        op2_filename_out = folder / 'static_solid_shell_bar_out.op2'
        op2_filename_debug_out = folder / 'static_solid_shell_bar_out.debug.out'
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_5(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'mode_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'mode_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'mode_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'mode_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            #'*_strain_energy',
            'grid_point_forces',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results,
                            log=log, )

        op2.write_op2(op2_filename_out) #is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_6(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'transient_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'transient_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'transient_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'transient_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = ['grid_point_forces']
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, log=log,
                            exclude_results=exclude_results)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_7(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'sol_101_elements')
        op2_filename = os.path.join(folder, 'freq_solid_shell_bar.op2')
        op2_filename_debug = os.path.join(folder, 'freq_solid_shell_bar.debug.out')
        op2_filename_out = os.path.join(folder, 'freq_solid_shell_bar_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'freq_solid_shell_bar_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_elements_1(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='warning', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'elements')
        op2_filename = os.path.join(folder, 'freq_elements.op2')
        op2_filename_debug = os.path.join(folder, 'freq_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'freq_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'freq_elements_out.debug.out')

        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_elements_2(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'elements')
        op2_filename = os.path.join(folder, 'freq_elements2.op2')
        op2_filename_debug = os.path.join(folder, 'freq_elements2.debug.out')
        op2_filename_out = os.path.join(folder, 'freq_elements_out2.op2')
        op2_filename_debug_out = os.path.join(folder, 'freq_elements_out2.debug.out')

        exclude_results = [
            'force.ctria6_force', 'force.ctriar_force', 'force.cshear_force',
            'force.cvisc_force', 'modal_contribution.cshear_stress',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, log=log)
        #print(op2.get_op2_stats())
        #from pyNastran.utils import object_stats
        #op2.op2_results.modal_contribution.celas1_stress = {}
        #op2.op2_results.modal_contribution.celas2_stress = {}
        #op2.op2_results.modal_contribution.celas3_stress = {}

        #op2.op2_results.modal_contribution.ctube_stress = {}
        #op2.op2_results.modal_contribution.crod_stress = {}
        #op2.op2_results.modal_contribution.conrod_stress = {}

        #op2.op2_results.modal_contribution.ctria3_stress = {}
        #op2.op2_results.modal_contribution.cquad4_stress = {}
        #op2.op2_results.modal_contribution.ctria6_stress = {}
        #op2.op2_results.modal_contribution.cquad8_stress = {}
        op2.op2_results.modal_contribution.ctriar_composite_stress = {}
        op2.op2_results.modal_contribution.cquadr_composite_stress = {}
        op2.op2_results.modal_contribution.cquad4_composite_stress = {}
        op2.op2_results.modal_contribution.ctria3_composite_stress = {}

        # ------------------------------------------------
        #op2.op2_results.modal_contribution.celas1_strain = {}
        #op2.op2_results.modal_contribution.celas2_strain = {}
        #op2.op2_results.modal_contribution.celas3_strain = {}

        #op2.op2_results.modal_contribution.ctube_strain = {}
        #op2.op2_results.modal_contribution.crod_strain = {}
        #op2.op2_results.modal_contribution.conrod_strain = {}

        #op2.op2_results.modal_contribution.ctria3_strain = {}
        #op2.op2_results.modal_contribution.cquad4_strain = {}
        #op2.op2_results.modal_contribution.ctria6_strain = {}
        #op2.op2_results.modal_contribution.cquad8_strain = {}
        op2.op2_results.modal_contribution.ctriar_composite_strain = {}
        op2.op2_results.modal_contribution.cquadr_composite_strain = {}
        op2.op2_results.modal_contribution.cquad4_composite_strain = {}
        op2.op2_results.modal_contribution.ctria3_composite_strain = {}
        #print(object_stats(op2.op2_results.modal_contribution))
        #aa

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        unused_op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        #op2.assert_op2_equal(op2b,
                             #skip_results=['params', ],
                             #stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_elements_3(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = MODEL_PATH / 'elements'
        op2_filename = os.path.join(folder, 'freq_random_elements.op2')
        op2_filename_debug = os.path.join(folder, 'freq_random_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'freq_random_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'freq_random_elements_out.debug.out')

        exclude_results = [
            'force.ctria6_force', 'force.ctriar_force', 'force.cshear_force',
            'force.cvisc_force', 'stress.cshear_stress', '*strain_energy',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results,
                            xref=False, log=log)
        op2.safe_cross_reference()

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out,
                             xref=False, log=log)
        op2b.safe_cross_reference()

        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_elements_4(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'elements')
        op2_filename = os.path.join(folder, 'modes_complex_elements.op2')
        op2_filename_debug = os.path.join(folder, 'modes_complex_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'modes_complex_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'modes_complex_elements_out.debug.out')

        exclude_results = [
            'force.ctria6_force', 'force.ctriar_force', 'force.cshear_force',
            'force.cvisc_force',
            'stress.cshear_stress',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_write_elements_5(self):
        """tests basic op2 writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = MODEL_PATH / 'elements'
        op2_filename = os.path.join(folder, 'time_elements.op2')
        op2_filename_debug = os.path.join(folder, 'time_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'time_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'time_elements_out.debug.out')
        #model = os.path.splitext(op2_filename)[0]

        exclude_results = [
            'force.cshear_force', 'force.cvisc_force',
            'grid_point_forces', '*strain_energy',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, log=log)

        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_thermal_1(self):
        """tests basic op2 thermal writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'elements')
        op2_filename = os.path.join(folder, 'time_thermal_elements.op2')
        op2_filename_debug = os.path.join(folder, 'time_thermal_elements.debug.out')
        op2_filename_out = os.path.join(folder, 'time_thermal_elements_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'time_thermal_elements.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            'thermal_load.chbdye_thermal_load',
            'thermal_load.chexa_thermal_load',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, log=log)
        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    def test_thermal_2(self):
        """tests basic op2 thermal writing"""
        log = SimpleLogger(level='info', encoding='utf-8')
        folder = os.path.join(MODEL_PATH, 'other')
        op2_filename = os.path.join(folder, 'hd15306.op2')
        op2_filename_debug = os.path.join(folder, 'hd15306.debug.out')
        op2_filename_out = os.path.join(folder, 'hd15306_out.op2')
        op2_filename_debug_out = os.path.join(folder, 'hd15306_out.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        exclude_results = [
            'thermal_load.chbdyg_thermal_load',
            'thermal_load.crod_thermal_load',
            'thermal_load.cquad4_thermal_load',
        ]
        op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            exclude_results=exclude_results, debug=True, log=log)
        str(op2.get_op2_stats(short=True))
        op2.write_op2(op2_filename_out) #, is_mag_phase=False)
        op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out, log=log)
        op2.assert_op2_equal(op2b,
                             skip_results=['params', ],
                             stop_on_failure=True, debug=False)
        os.remove(op2_filename_debug_out)

    #def test_thermal_3(self):
        #"""tests basic op2 thermal writing"""
        #folder = os.path.join(MODEL_PATH, 'other')
        #op2_filename = os.path.join(folder, 'ofprand1.op2')
        #op2_filename_debug = os.path.join(folder, 'ofprand1.debug.out')
        #op2_filename_out = os.path.join(folder, 'ofprand1_out.op2')
        #op2_filename_debug_out = os.path.join(folder, 'ofprand1.debug.out')
        #debug_file = 'solid_bending.debug.out'
        #model = os.path.splitext(op2_filename)[0]
        #debug_file = model + '.debug.out'

        #exclude_results = [
            #'thermal_load.chbdyg_thermal_load',
            #'thermal_load.crod_thermal_load',
            #'thermal_load.cquad4_thermal_load',
            #'thermal_load.chbdye_thermal_load',
            #'thermal_load.chexa_thermal_load',
        #]
        #op2 = read_op2_geom(op2_filename, debug_file=op2_filename_debug,
                            #exclude_results=exclude_results,
                            ##include_results='eigenvectors',
                            ##include_results=['crod_stress', 'cbar_stress'],
                            ##include_results=['crod_force', 'cbar_force'],
                            ##include_results='force',
                            ##include_results='stress',
                            #)
        #print(op2.get_op2_stats(short=True))
        #op2.write_op2(op2_filename_out, is_mag_phase=False)
        #op2b = read_op2_geom(op2_filename_out, debug_file=op2_filename_debug_out)
        ##op2b = read_op2(op2_filename_out, debug_file=op2_filename_debug_out)
        #op2.assert_op2_equal(op2b,
                             #skip_results=['params', ],
                             #stop_on_failure=True, debug=False)

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
