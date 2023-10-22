#import os
import unittest
import pathlib
#import numpy as np

import pyNastran
from pyNastran.dev.bdf_vectorized3.test.test_bdf import (
    read_bdfv, main as test_bdf, IS_PYTABLES)
from pyNastran.dev.op2_vectorized3.test.test_op2 import main as test_op2

PKG_PATH = pathlib.Path(pyNastran.__path__[0])
TEST_PATH = PKG_PATH / 'bdf' / 'test'
MODEL_PATH = PKG_PATH / '..' / 'models'

#from cpylog import SimpleLogger


class TestModels(unittest.TestCase):
    def test_h5_freq(self):
        h5_filename = MODEL_PATH / 'elements' / 'freq_elements.h5'
        args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)

        #args2 = ['test_op2', str(h5_filename), '-ctfo', '--quiet']
        #test_op2(args2, show_args=False)

    def test_h5_transient(self):
        #bdf_filename = MODEL_PATH / 'elements' / 'time_elements.bdf'
        h5_filename = MODEL_PATH / 'elements' / 'time_elements.h5'
        #args = ['test_bdf', str(bdf_filename), '--skip_nominal', '--quiet']
        #test_bdf(args, show_args=False)

        args2 = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args2, show_args=False)

        args3 = ['test_op2', str(h5_filename), '-ctfg', '--quiet'] # o
        test_op2(args3, show_args=False)

    def test_h5_transient_thermal(self):
        h5_filename = MODEL_PATH / 'elements' / 'time_thermal_elements.h5'

        args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        #test_bdf(args, show_args=False)

        args2 = ['test_op2', str(h5_filename), '-ctfo', '--quiet']
        #test_op2(args2, show_args=False)

    def test_h5_static(self):
        #bdf_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.bdf'
        #args = ['test_bdf', str(bdf_filename), '--quiet']
        #test_bdf(args, show_args=False)

        if IS_PYTABLES:
            h5_filename = MODEL_PATH / 'sol_101_elements' / 'static_solid_shell_bar.h5'
            #args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
            #test_bdf(args, show_args=False)
            args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet'] # o
            test_op2(args2, show_args=False)

    def test_h5_modes(self):
        h5_filename = MODEL_PATH / 'elements' / 'modes_elements.h5'
        args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)
        args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet'] # o
        test_op2(args2, show_args=False)

    def test_h5_modes_complex(self):
        h5_filename = MODEL_PATH / 'elements' / 'modes_complex_elements.h5'
        args1 = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args1, show_args=False)
        args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet'] # o
        test_op2(args2, show_args=False)

    def _test_h5_buckling(self):
        h5_filename = MODEL_PATH / 'sol_101_elements' / 'buckling_solid_shell_bar.h5'
        args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)
        args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet'] # o
        test_op2(args2, show_args=False)

    def _test_h5_buckling2(self):
        h5_filename = MODEL_PATH / 'sol_101_elements' / 'buckling2_solid_shell_bar.h5'
        args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)

        args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet'] # o
        test_op2(args2, show_args=False)

    #def test_h5_bwb(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\bwb_saero_saved.h5'
        #args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        #test_bdf(args, show_args=False)
        #args2 = ['test_op2', str(h5_filename), '-gctf', '--quiet']
        #test_op2(args2, show_args=False)

    #def test_h5_mode_echo(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\msc\mode_echo.h5'
        #args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        #test_bdf(args, show_args=False)
        #args2 = ['test_op2', str(h5_filename), '-gctfo', '--quiet']
        #test_op2(args2, show_args=False)

    #def _test_h5_mode_cfast(self):
        #h5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\msc\test_model_cfast.h5'
        #args = ['test_bdf', str(h5_filename), '--skip_nominal', '--quiet']
        #test_bdf(args, show_args=False)

    def test_beam_modes1(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'beam_modes.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_beam_modes2(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'cbarao_cbeam_static.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_beam_modes3(self):
        bdf_filename = MODEL_PATH / 'beam_modes' / 'cbarao_cbeam_modes.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def _test_petite(self):
        bdf_filename = MODEL_PATH / 'modele_petite_zone' / 'modele_petite_zone.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_random1(self):
        bdf_filename = MODEL_PATH / 'random' / 'random_test_bar_plus_tri.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_random2(self):
        bdf_filename = MODEL_PATH / 'random' / 'rms_tri_oesrmx1.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_freq_sine(self):
        bdf_filename = MODEL_PATH / 'freq_sine' / 'good_sine.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_support_structure(self):
        bdf_filename = MODEL_PATH / 'support_structure' / 'W1000BOstat.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_fsi(self):
        bdf_filename = MODEL_PATH / 'fsi' / 'fsi.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_bwb(self):
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_solid_bending(self):
        bdf_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

        op2_filename = MODEL_PATH / 'solid_bending' / 'solid_bending.op2'
        args = ['test_op2', str(op2_filename), '-gnt']
        test_op2(args, show_args=False)

    def test_isat1(self):
        bdf_filename = MODEL_PATH / 'iSat' / 'ISat_Dploy_Sm.dat'
        #from pyNastran.bdf.bdf import read_bdf as read_bdf_old
        #model_old = read_bdf_old(bdf_filename)
        #elem_old = model_old.elements[3910]
        #elem_old.get_axes(model_old)

        #modelv = read_bdfv(bdf_filename)
        #cbar = modelv.cbar.slice_card_by_id(3910)
        #print(cbar.write())

        #m = np.array([
        #    [ 0.70710678 -0.70710678  0.        ],
        #    [ 0.          0.         -1.        ],
        #    [ 0.70710678  0.70710678  0.        ], ])
        #p = np.array([-2.76770038,  1.72858252, -1.59188988])

        #p3 = p @ matrix
        #p3 = np.array([-3.08269584,  0.83142358, -1.72858252])
        #xyz1, xyz2 = cbar.get_xyz()
        #v, cd = cbar.get_bar_vector(xyz1)
        #cbar.get_axes(xyz1, xyz2)

        #builtins.RuntimeError: CBAR eid=3910
        #v=[-3.08269584 0.83142358 -1.72858252] v_new=[-2.76770038 1.72858252 -1.59188988]

        #args = ['test_bdf', str(bdf_filename), '--quiet', '--skip_nominal']

        #xform2 = np.array([[
        #    [ 0.70710678, -0.70710678,  0.        ],
        #    [ 0.        ,  0.        , -1.        ],
        #    [ 0.70710678,  0.70710678,  0.        ]]])

        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_isat2(self):
        bdf_filename = MODEL_PATH / 'iSat' / 'iSat_launch_100Hz.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_contact_1(self):
        bdf_filename = MODEL_PATH / 'contact' / 'contact.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet', '--skip_nominal']
        test_bdf(args, show_args=False)
    def test_contact_2(self):
        bdf_filename = MODEL_PATH / 'contact' / '2bars_shell_s-contact.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_contact_3(self):
        bdf_filename = MODEL_PATH / 'contact' / '2bars_shell_s-contact.dat'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_elements_static(self):
        bdf_filename = MODEL_PATH / 'elements' / 'static_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

        #if IS_PYTABLES:
            #h5_filename = MODEL_PATH / 'elements' / 'static_elements.h5'
            #args = ['test_bdf', str(h5_filename), '--quiet']
            #test_bdf(args, show_args=False)
    def test_elements_modes(self):
        bdf_filename = MODEL_PATH / 'elements' / 'modes_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
        #if IS_PYTABLES:
            #h5_filename = MODEL_PATH / 'elements' / 'modes_elements.h5'
            #args = ['test_bdf', str(h5_filename), '--quiet']
            #test_bdf(args, show_args=False)

    def test_elements_transient(self):
        bdf_filename = MODEL_PATH / 'elements' / 'time_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_elements_freq(self):
        bdf_filename = MODEL_PATH / 'elements' / 'freq_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_elements_freq2(self):
        bdf_filename = MODEL_PATH / 'elements' / 'freq_elements2.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
        if IS_PYTABLES:
            h5_filename = MODEL_PATH / 'elements' / 'freq_elements2.h5'
            args = ['test_bdf', str(h5_filename), '--quiet']
            test_bdf(args, show_args=False)

    def test_elements_loadstep(self):
        bdf_filename = MODEL_PATH / 'elements' / 'loadstep_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
        #if IS_PYTABLES:
            #h5_filename = MODEL_PATH / 'elements' / 'loadstep_elements.h5'
            #args = ['test_bdf', str(h5_filename), '--quiet']
            #test_bdf(args, show_args=False)
    def test_elements_modes_complex(self):
        bdf_filename = MODEL_PATH / 'elements' / 'modes_complex_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_elements_transient_thermal(self):
        bdf_filename = MODEL_PATH / 'elements' / 'time_thermal_elements.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

    def test_thermal_1(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'hd15901.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)
    def test_thermal_2(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'htflw47.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)
    def test_thermal_3(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_elements.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)
    def test_thermal_4(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_elements2.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal', '--quiet']
        test_bdf(args, show_args=False)
    def test_thermal_5(self):
        bdf_filename = MODEL_PATH / 'thermal' / 'thermal_test_153.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def _test_sol200_1(self):
        bdf_filename = MODEL_PATH / 'sol200' / 'd200obus.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)
    def _test_sol200_2(self):
        bdf_filename = MODEL_PATH / 'sol200' / 'model_200.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def test_transfer_function(self):
        bdf_filename = MODEL_PATH / 'transfer_function' / 'actuator_tf_modeling.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
    def _test_aero_1(self):
        bdf_filename = MODEL_PATH / 'aero' / 'freedlm' / 'freedlm.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)
        if IS_PYTABLES:
            h5_filename = MODEL_PATH / 'aero' / 'freedlm' / 'freedlm_msc.h5'
            args = ['test_bdf', str(h5_filename), '--quiet']
            test_bdf(args, show_args=False)
    def _test_aero_2(self):
        bdf_filename = MODEL_PATH / 'aero' / 'bah_plane' / 'bah_plane.bdf'
        args = ['test_bdf', str(bdf_filename), '--quiet']
        test_bdf(args, show_args=False)

    def test_other_1(self):
        bdf_filename = MODEL_PATH / 'other' / 'ac10707a.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal']
        test_bdf(args, show_args=False)
    def test_other_2(self):
        bdf_filename = MODEL_PATH / 'other' / 'dbxdra2.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal']
        test_bdf(args, show_args=False)

    def _test_overflow(self):
        # overflow
        bdf_filename = MODEL_PATH / 'other' / 'sdr11se_s2dclg.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal']
        test_bdf(args, show_args=False)

    def test_other_3(self):
        # missing GRID card
        #bdf_filename = MODEL_PATH / 'other' / 'ofprand1.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # missing GRID cards
        #bdf_filename = MODEL_PATH / 'other' / 'ac10707a.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # TSTEP
        #bdf_filename = MODEL_PATH / 'other' / 'ac10901a.bdf'
        #args = ['test_bdf', str(bdf_filename), '--skip_nominal']
        #test_bdf(args, show_args=False)

        # replication error
        #bdf_filename = MODEL_PATH / 'other' / 'v10111.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'ar29sadl.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        # spoint centroid bug
        #bdf_filename = MODEL_PATH / 'other' / 'randvar2.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'v12902.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        # replication
        #bdf_filename = MODEL_PATH / 'other' / 'mne7a.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'sdbush10.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'v10112.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        # bad deck sections
        #bdf_filename = MODEL_PATH / 'other' / 'cbus129.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # incorrect cbeam area
        bdf_filename = MODEL_PATH / 'other' / 'api3.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_mass', '--skip_nominal']
        test_bdf(args, show_args=False)

        # incorrect cbeam area
        ##bdf_filename = MODEL_PATH / 'other' / 'v10601s.bdf'
        ##args = ['test_bdf', str(bdf_filename)]
        ##test_bdf(args, show_args=False)

        ## PLOT MODAL, 0, SET 1, ORIGIN 1, SET 2, SYMBOL 3
        ##bdf_filename = MODEL_PATH / 'other' / 'sbuckl2a.bdf'
        ##args = ['test_bdf', str(bdf_filename)]
        ##test_bdf(args, show_args=False)

        # bad parsing of GRID
        bdf_filename = MODEL_PATH / 'other' / 'dbxdra2.bdf'
        args = ['test_bdf', str(bdf_filename), '--skip_nominal']
        test_bdf(args, show_args=False)

        # ADAPT - bad pshell mass/area
        #bdf_filename = MODEL_PATH / 'other' / 'phsflux4.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # bad cpenta volume
        #bdf_filename = MODEL_PATH / 'other' / 'cc508a.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # missing GRID for superelement
        #bdf_filename = MODEL_PATH / 'other' / 'see101nd.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'see101ta.bdf'
        args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # CGEN model - no nodes
        #bdf_filename = MODEL_PATH / 'other' / 'gpst17.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # = parsing - replication
        #bdf_filename = MODEL_PATH / 'other' / 'cqra00366.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # = parsing - replication
        #bdf_filename = MODEL_PATH / 'other' / 'dbxdra7.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # axisymmetric
        #bdf_filename = MODEL_PATH / 'other' / 'ehbus69.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # missing nodes crash
        #bdf_filename = MODEL_PATH / 'other' / 'tst1d3.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'trncomp12.bdf'
        args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # EGRID/SPCG
        #bdf_filename = MODEL_PATH / 'other' / 'tr1091x.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # EGRID/SPCG
        #bdf_filename = MODEL_PATH / 'other' / 'tr1091x.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        # EIGC parsing
        #bdf_filename = MODEL_PATH / 'other' / 'rot063akd2s_107.bdf'
        #args = ['test_bdf', str(bdf_filename)]
        #test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'other' / 'htrussx.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        bdf_filename = MODEL_PATH / 'nx' / 'contact_model.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)

        # missing properties - PCOMPS
        bdf_filename = MODEL_PATH / 'nx' / 'composite_solids' / 'test.bdf'
        args = ['test_bdf', str(bdf_filename)]
        test_bdf(args, show_args=False)


if __name__ == '__main__':
    import warnings
    warnings.simplefilter('always', DeprecationWarning)
    unittest.main()
