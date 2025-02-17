import os
import unittest
import pathlib
import numpy as np

import pyNastran
from pyNastran.op2.op2_geom import read_op2_geom
from pyNastran.converters.nastran.gui.result_objects.displacement_results import DisplacementResults2
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2

#from pyNastran.dev.bdf_vectorized3.test.test_bdf import main as test_bdf
#from pyNastran.dev.op2_vectorized3.test.test_op2 import main as test_op2

PKG_PATH = pathlib.Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'bdf' / 'test'
MODEL_PATH = PKG_PATH / '..' / 'models'

import vtkmodules

#from cpylog import SimpleLogger
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.dev.bdf_vectorized3.nastran_io3 import Nastran3

class NastranGUI(Nastran3, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        Nastran3.__init__(self, self)
        self.build_fmts(['nastran3'], stop_on_failure=True)

def run_nastran_gui(filename: str,
                    load_results: bool=True):
    assert os.path.exists(filename), filename
    filename = str(filename)
    test = NastranGUI()
    test.load_nastran3_geometry(filename)
    if filename.lower().endswith(('.op2', '.h5')) and load_results:
        test.load_nastran3_results(filename)
    test.cycle_results()
    test.on_rcycle_results()

def run_nastran_gui_results(input_filename: str,
                            output_filename: str):
    assert os.path.exists(input_filename), input_filename
    assert os.path.exists(output_filename), output_filename
    input_filename = str(input_filename)
    output_filename = str(output_filename)
    test = NastranGUI()
    test.load_nastran3_geometry(input_filename)
    test.load_nastran3_results(output_filename)
    test.cycle_results()
    test.on_rcycle_results()


class TestGuiModels(unittest.TestCase):
    def test_h5_freq(self):
        h5_filename = MODEL_PATH / 'elements' / 'freq_elements.h5'
        run_nastran_gui(h5_filename)

    def test_h5_transient(self):
        h5_filename = MODEL_PATH / 'elements' / 'time_elements.h5'
        run_nastran_gui(h5_filename)

    def test_h5_transient_thermal(self):
        h5_filename = MODEL_PATH / 'elements' / 'time_thermal_elements.h5'
        run_nastran_gui(h5_filename)

    def test_h5_static(self):
        h5_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.h5'
        run_nastran_gui(h5_filename)

    def test_h5_modes(self):
        h5_filename = MODEL_PATH / 'elements' / 'modes_elements.h5'
        run_nastran_gui(h5_filename)

    def test_h5_modes_complex(self):
        h5_filename = MODEL_PATH / 'elements' / 'modes_complex_elements.h5'
        run_nastran_gui(h5_filename)

    def _test_h5_buckling(self):
        h5_filename = MODEL_PATH / 'sol_101_elements' / 'buckling_solid_shell_bar.h5'
        run_nastran_gui(h5_filename)

    def _test_h5_buckling2(self):
        h5_filename = MODEL_PATH / 'sol_101_elements' / 'buckling2_solid_shell_bar.h5'
        run_nastran_gui(h5_filename)

    #def test_h5_bwb(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\bwb_saero_saved.h5'
        #run_nastran_gui(h5_filename)

    #def test_h5_mode_echo(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\msc\mode_echo.h5'
        #run_nastran_gui(h5_filename)

    #def _test_h5_mode_cfast(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\msc\test_model_cfast.h5'
        #run_nastran_gui(h5_filename)

    def test_beam_modes1(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'beam_modes.dat'
        run_nastran_gui(bdf_filename)
    def test_beam_modes2(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'cbarao_cbeam_static.bdf'
        run_nastran_gui(bdf_filename)
    def test_beam_modes3(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'cbarao_cbeam_modes.bdf'
        run_nastran_gui(bdf_filename)

    #def test_petite(self):
        #bdf_filename = MODEL_PATH / 'modele_petite_zone' / 'modele_petite_zone.dat'
        #args = ['test_bdf', str(bdf_filename), '--quiet']

    def test_random1(self):
        bdf_filename = MODEL_PATH / 'random' / 'random_test_bar_plus_tri.bdf'
        run_nastran_gui(bdf_filename)
    def test_random2(self):
        bdf_filename = MODEL_PATH / 'random' / 'rms_tri_oesrmx1.bdf'
        run_nastran_gui(bdf_filename)

    def test_freq_sine(self):
        bdf_filename = MODEL_PATH / 'freq_sine' / 'good_sine.dat'
        run_nastran_gui(bdf_filename)
    def test_support_structure(self):
        bdf_filename = MODEL_PATH / 'support_structure' / 'W1000BOstat.dat'
        run_nastran_gui(bdf_filename)
    def test_fsi(self):
        bdf_filename = MODEL_PATH / 'fsi' / 'fsi.bdf'
        run_nastran_gui(bdf_filename)
    def test_bwb(self):
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        run_nastran_gui(bdf_filename)
    def test_solid_bending(self):
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        op2_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.op2'
        #run_nastran_gui(bdf_filename)
        run_nastran_gui_results(bdf_filename, op2_filename)

    def test_solid_bending_disp(self):
        op2_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.op2'
        model = read_op2_geom(op2_filename, debug=None)
        subcase_id = 1

        nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array()
        node_id = nid_cp_cd[:, 0]
        assert xyz_cid0.shape == (72, 3), xyz_cid0.shape
        disp = DisplacementResults2.add_from_displacements(
            model.displacements, subcase_id, node_id, xyz_cid0,
        )
        #disp.set_coords(model.coords)
        #disp.set_coordiate_system(1)
        disp.set_translations([0, 1, 2], 'Magnitude')
        itime = 0
        res_name = 'junk'
        xyz, deflected_xyz = disp.get_vector_result(
            itime, res_name, scale=1.0, return_dense=True)
        assert xyz.shape == (72, 3), xyz.shape
        assert deflected_xyz.shape == (72, 3), deflected_xyz.shape

        #title = 'title'
        tetra_stress = model.op2_results.stress.ctetra_stress[subcase_id]
        element_id = np.array(list(model.elements))
        element_id.sort()
        #cases = [tetra_stress]

        #result_dict = {}
        # stress = SolidStrainStressResults2.add_stress(
        #     subcase_id, model, node_id, element_id, cases,
        #     result_dict, title)
        assert len(model.op2_results.stress.ctetra_stress) > 0
        stress = SolidStrainStressResults2.add_stress(
            subcase_id, model, model, element_id)
        #print(f'layer_indices={stress.layer_indices}')
        #print(f'min_max_method={stress.min_max_method!r}')
        #print(f'nodal_combine={stress.nodal_combine!r}')

        for result_name in ['xx', 'yy', 'zz', 'xy', 'yz', 'xz',
                            'von_mises']:
            case_tuple = stress.get_case_tuple(itime, result_name)

            #iresult = 0
            #header = ''
            #case_tuple = (itime, iresult, header)
            stress.set_to_centroid()
            fringe = stress.get_fringe_result(itime, case_tuple)
            assert len(fringe) == len(element_id), fringe.shape
            out = stress.get_annotation(itime, case_tuple)
        assert out == 'Solid von Mises (Centroid, ): von Mises', f'out={out!r}'

        stress.set_to_node(nodal_combine='Absolute Max')
        fringe2 = stress.get_fringe_result(itime, case_tuple)
        assert len(fringe2) == len(node_id), fringe2.shape
        out = stress.get_annotation(itime, case_tuple)
        assert out == 'Solid von Mises (Corner, Absolute Max, ): von Mises', f'out={out!r}'

        stress.set_to_node(nodal_combine='Mean')
        out = stress.get_annotation(itime, case_tuple)
        assert out == 'Solid von Mises (Corner, Mean, ): von Mises', f'out={out!r}'
        fringe3 = stress.get_fringe_result(itime, case_tuple)
        assert not(np.array_equal(fringe2, fringe3))

    def test_isat1(self):
        bdf_filename = MODEL_PATH / 'iSat' / 'ISat_Dploy_Sm.dat'
        run_nastran_gui(bdf_filename)
    def test_isat2(self):
        bdf_filename = MODEL_PATH / 'iSat' / 'iSat_launch_100Hz.dat'
        run_nastran_gui(bdf_filename)

    def test_contact_1(self):
        bdf_filename = MODEL_PATH / 'contact' / 'contact.bdf'
        run_nastran_gui(bdf_filename)
    def test_contact_2(self):
        bdf_filename = MODEL_PATH / 'contact' / '2bars_shell_s-contact.dat'
        run_nastran_gui(bdf_filename)
    def test_contact_3(self):
        bdf_filename = MODEL_PATH / 'contact' / '2bars_shell_s-contact.dat'
        run_nastran_gui(bdf_filename)

    def test_elements1(self):
        bdf_filename = MODEL_PATH / 'elements' / 'static_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'static_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements2(self):
        bdf_filename = MODEL_PATH / 'elements' / 'modes_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'modes_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements3(self):
        bdf_filename = MODEL_PATH / 'elements' / 'time_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'time_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements4(self):
        bdf_filename = MODEL_PATH / 'elements' / 'freq_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'freq_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements5(self):
        bdf_filename = MODEL_PATH / 'elements' / 'freq_elements2.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'freq_elements2.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements6(self):
        bdf_filename = MODEL_PATH / 'elements' / 'loadstep_elements.bdf'
        #h5_filename = MODEL_PATH / 'elements' / 'loadstep_elements.h5'
        run_nastran_gui(bdf_filename)
        #run_nastran_gui(h5_filename)
    def test_elements7(self):
        bdf_filename = MODEL_PATH / 'elements' / 'modes_complex_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'modes_complex_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)
    def test_elements8(self):
        bdf_filename = MODEL_PATH / 'elements' / 'time_thermal_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'time_thermal_elements.h5'
        run_nastran_gui(bdf_filename)
        run_nastran_gui(h5_filename)

    def test_thermal_1(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'hd15901.bdf'
        run_nastran_gui(bdf_filename)
    def test_thermal_2(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'htflw47.bdf'
        run_nastran_gui(bdf_filename)
    def test_thermal_3(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_elements.bdf'
        run_nastran_gui(bdf_filename)
    def test_thermal_4(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_elements2.bdf'
        run_nastran_gui(bdf_filename)
    def test_thermal_5(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_test_153.bdf'
        run_nastran_gui(bdf_filename)

    #def test_sol200_1(self):
        #bdf_filename = MODEL_PATH / 'sol200' / 'd200obus.bdf'
    def test_sol200_2(self):
        bdf_filename = MODEL_PATH / 'sol200' / 'model_200.bdf'
        run_nastran_gui(bdf_filename)
    def test_transfer_function(self):
        bdf_filename = MODEL_PATH / 'transfer_function' / 'actuator_tf_modeling.bdf'
        run_nastran_gui(bdf_filename)
    def test_aero_1(self):
        #h5_filename = MODEL_PATH / 'aero' / 'freedlm' / 'freedlm_msc.h5'
        #run_nastran_gui(h5_filename)
        bdf_filename = MODEL_PATH / 'aero' / 'freedlm' / 'freedlm.bdf'
        run_nastran_gui(bdf_filename)
    def test_aero_2(self):
        bdf_filename = MODEL_PATH / 'aero' / 'bah_plane' / 'bah_plane.bdf'
        run_nastran_gui(bdf_filename)

    def test_other_1(self):
        bdf_filename = MODEL_PATH / 'other' / 'ac10707a.bdf'
        run_nastran_gui(bdf_filename)
    def _test_other_2(self):
        bdf_filename = MODEL_PATH / 'other' / 'dbxdra2.bdf'
        run_nastran_gui(bdf_filename)

    def test_other_3(self):
        # missing GRID card
        #bdf_filename = MODEL_PATH / 'other' / 'ofprand1.bdf'

        # missing GRID cards
        #bdf_filename = MODEL_PATH / 'other' / 'ac10707a.bdf'

        # TSTEP
        #bdf_filename = MODEL_PATH / 'other' / 'ac10901a.bdf'

        # replication error
        #bdf_filename = MODEL_PATH / 'other' / 'v10111.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'ar29sadl.bdf'
        run_nastran_gui(bdf_filename)

        # spoint centroid bug
        #bdf_filename = MODEL_PATH / 'other' / 'randvar2.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'v12902.bdf'
        run_nastran_gui(bdf_filename)

        # replication
        #bdf_filename = MODEL_PATH / 'other' / 'mne7a.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'sdbush10.bdf'
        run_nastran_gui(bdf_filename)

        bdf_filename = MODEL_PATH / 'other' / 'v10112.bdf'
        run_nastran_gui(bdf_filename)

        # bad deck sections
        #bdf_filename = MODEL_PATH / 'other' / 'cbus129.bdf'

        # incorrect cbeam area
        bdf_filename = MODEL_PATH / 'other' / 'api3.bdf'
        #run_nastran_gui(bdf_filename)

        # incorrect cbeam area
        ##bdf_filename = MODEL_PATH / 'other' / 'v10601s.bdf'

        ## PLOT MODAL, 0, SET 1, ORIGIN 1, SET 2, SYMBOL 3
        ##bdf_filename = MODEL_PATH / 'other' / 'sbuckl2a.bdf'

        # bad parsing of GRID
        bdf_filename = MODEL_PATH / 'other' / 'dbxdra2.bdf'

        # ADAPT - bad pshell mass/area
        #bdf_filename = MODEL_PATH / 'other' / 'phsflux4.bdf'

        # bad cpenta volume
        #bdf_filename = MODEL_PATH / 'other' / 'cc508a.bdf'

        # missing GRID for superelement
        #bdf_filename = MODEL_PATH / 'other' / 'see101nd.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'see101ta.bdf'
        run_nastran_gui(bdf_filename)

        # CGEN model - no nodes
        #bdf_filename = MODEL_PATH / 'other' / 'gpst17.bdf'

        # = parsing - replication
        #bdf_filename = MODEL_PATH / 'other' / 'cqra00366.bdf'

        # = parsing - replication
        #bdf_filename = MODEL_PATH / 'other' / 'dbxdra7.bdf'

        # axisymmetric
        #bdf_filename = MODEL_PATH / 'other' / 'ehbus69.bdf'

        # missing nodes crash
        #bdf_filename = MODEL_PATH / 'other' / 'tst1d3.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'trncomp12.bdf'
        run_nastran_gui(bdf_filename)

        # EGRID/SPCG
        #bdf_filename = MODEL_PATH / 'other' / 'tr1091x.bdf'

        # overflow
        bdf_filename = MODEL_PATH / 'other' / 'sdr11se_s2dclg.bdf'
        run_nastran_gui(bdf_filename)

        # EIGC parsing
        #bdf_filename = MODEL_PATH / 'other' / 'rot063akd2s_107.bdf'

        bdf_filename = MODEL_PATH / 'other' / 'htrussx.bdf'
        run_nastran_gui(bdf_filename)

        bdf_filename = MODEL_PATH / 'nx' / 'contact_model.bdf'
        run_nastran_gui(bdf_filename)

        # missing properties - PCOMPS
        bdf_filename = MODEL_PATH / 'nx' / 'composite_solids' / 'test.bdf'
        run_nastran_gui(bdf_filename)
