from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.gui.gui_utils.utils import load_csv, load_deflection_csv, load_user_geom
from pyNastran.gui.gui_utils.write_gif import setup_animation
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class GuiUtils(unittest.TestCase):
    def test_gui_csv_01(self):
        """tests solid_bending.txt"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending.txt')
        with self.assertRaises(ValueError):
            load_deflection_csv(csv_filename) # bad shape
        load_csv(csv_filename)

    def test_gui_csv_02(self):
        """tests solid_bending_multi_node.txt"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node.txt')
        with self.assertRaises(ValueError):
            load_deflection_csv(csv_filename) # bad shape
        load_csv(csv_filename)


    def test_gui_csv_03a(self):
        """tests solid_bending_multi_node.csv with deflection loader"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node.csv')
        with self.assertRaises(ValueError):
            load_deflection_csv(csv_filename) # bad shape

    def test_gui_csv_03b(self):
        """tests solid_bending_multi_node.csv"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node.csv')
        load_csv(csv_filename)

    def test_gui_deflection_csv_01a(self):
        """tests solid_bending_multi_node.csv with deflection loader"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_node.csv')
        with self.assertRaises(ValueError):
            load_deflection_csv(csv_filename) # bad shape

    def test_gui_deflection_csv_01b(self):
        """tests solid_bending_multi_deflection_node.txt with deflection loader"""
        csv_filename = os.path.join(MODEL_PATH, 'solid_bending', 'solid_bending_multi_deflection_node.txt')
        load_deflection_csv(csv_filename)

    def test_gui_custom_geom_01(self):
        """tests custom_geom.csv"""
        csv_filename = os.path.join(MODEL_PATH, 'custom_geom.csv')
        load_user_geom(csv_filename)

    def test_animation_scale_1(self):
        """0 to scale"""
        scale = 1.0
        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, animation_profile='0 to scale to 0',
            fps=30)
        assert len(np.unique(icases)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales

    def test_animation_scale_2(self):
        """-scale to scale"""
        scale = 1.0
        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, animation_profile='-scale to scale',
            fps=30)
        assert len(np.unique(icases)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales

    def test_animation_scale_3(self):
        """-scale to scale to -scale"""
        scale = 1.0
        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, animation_profile='-scale to scale to -scale',
            fps=30)
        assert len(np.unique(icases)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales

    def test_animation_bad_1(self):
        """animate_scale/phase/time must be specified"""
        scale = 1.0
        with self.assertRaises(NotImplementedError):
            phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
                scale, istep=None,
                animate_scale=False, animate_phase=False, animate_time=False,
                icase=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=2.0,  animation_profile='-scale to scale',
                fps=30)

    def test_animation_phase_1(self):
        """phase plot"""
        scale = 1.0
        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=True, animate_time=False,
            icase=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        assert len(np.unique(icases)) == 1
        assert len(np.unique(scales)) == 1
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases[0], 0.), phases
        assert np.allclose(phases[-1], 354.), phases

    def test_animation_time_1(self):
        """time plot"""
        scale = 1.0
        phases, icases, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=False, animate_time=True,
            icase=42,
            icase_start=1, icase_end=10, icase_delta=2,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        assert np.array_equal(icases, [1, 3, 5, 7, 9]), icases
        assert len(np.unique(scales)) == 1, scales
        assert len(np.unique(isteps)) > 1, isteps
        assert len(np.unique(phases)) == 1, phases
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases.max(), 0.), phases


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
