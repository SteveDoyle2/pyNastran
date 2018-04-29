from __future__ import print_function
import os
import unittest
import numpy as np

import pyNastran
from pyNastran.gui.utils.load_results import (
    load_csv, load_deflection_csv, load_user_geom, create_res_obj)
from pyNastran.gui.utils.write_gif import setup_animation, make_two_sided, make_symmetric
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
        A, fmt_dict, headers = load_deflection_csv(csv_filename)
        result_type = 'node'

        header0 = headers[0]
        result0 = A[header0]
        #nrows = result0.shape[0]
        #assert nrows == self.nnodes, 'nrows=%s nnodes=%s' % (nrows, self.nnodes)
        header = header0
        islot = 0
        create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                       dim_max=1.0, xyz_cid0=None)

    def test_gui_custom_geom_01(self):
        """tests custom_geom.csv"""
        csv_filename = os.path.join(MODEL_PATH, 'custom_geom.csv')
        load_user_geom(csv_filename)

    def test_animation_scale_0(self):
        """0 to scale"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale',
            fps=5)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        np.testing.assert_almost_equal(scales, [0.0, 0.25*scale, 0.5*scale, 0.75*scale, scale])

    def test_animation_scale_1a(self):
        """0 to scale to 0"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to 0',
            fps=11)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        np.testing.assert_almost_equal(scales, [0., 0.5*scale, scale])

    def test_animation_scale_1b(self):
        """0 to scale to 0"""
        scale = 1.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to 0',
            fps=10)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        np.testing.assert_almost_equal(scales, [0., 0.5, 1.])

    def test_animation_scale_2(self):
        """-scale to scale"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale',
            fps=5)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        np.testing.assert_almost_equal(scales, [-scale, -0.5*scale, 0., 0.5*scale, scale])

    def test_animation_scale_3a(self):
        """-scale to scale to -scale"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=11)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        np.testing.assert_almost_equal(scales, [-scale, 0., scale])

    def test_animation_scale_3b(self):
        """-scale to scale to -scale"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=10)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        assert len(scales) == 3, scales
        np.testing.assert_almost_equal(scales, [-scale, 0., scale])

    def test_animation_scale_3c(self):
        """-scale to scale to -scale"""
        scale = 2.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=1)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        np.testing.assert_almost_equal(scales, [-scale, scale])
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        assert len(scales) == 2, scales

    def test_animation_scale_4a(self):
        """0 to scale to -scale to 0"""
        scale_atols = [(1., 0.00000001), (10.0, 0.00000001)]
        for scale, atol in scale_atols:
            phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
                scale, istep=None,
                animate_scale=True, animate_phase=False, animate_time=False,
                icase_disp=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=1.0, animation_profile='0 to scale to -scale to 0',
                fps=5, )
            assert len(np.unique(icases_disp)) == 1
            assert len(np.unique(phases)) == 1
            assert np.allclose(analysis_time, 1.0), analysis_time
            np.testing.assert_almost_equal(scales, [0., scale, 0., -scale])

            msg = ''
            if not np.allclose(scales.min(), -scale, atol=atol):
                msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
            if not np.allclose(scales.max(), scale, atol=atol):
                msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
            if msg:  # pragma: no cover
                raise ValueError(msg)

    def test_animation_scale_4b(self):
        """0 to scale to -scale to 0"""
        scale, atol = (1., 0.00000001)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to -scale to 0',
            fps=7, animate_in_gui=True)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        np.testing.assert_almost_equal(scales, [0., 0.5*scale, scale, 0.5*scale,
                                                0., -0.5*scale, -scale, -0.5*scale, 0.])

        scales, phases, icases_fringe, icases_disp, icases_vector, isteps = make_two_sided(
            scales, phases, icases_fringe, icases_disp, icases_vector, isteps, onesided)
        np.testing.assert_almost_equal(scales, [0., 0.5*scale, scale, 0.5*scale, 0., -0.5*scale, -scale, -0.5*scale,
                                                0., 0.5*scale, scale, 0.5*scale, 0., -0.5*scale, -scale, -0.5*scale, ] )

        msg = ''
        if not np.allclose(scales.min(), -scale, atol=atol):
            msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
        if not np.allclose(scales.max(), scale, atol=atol):
            msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
        if msg:  # pragma: no cover
            raise ValueError(msg)

    def test_animation_scale_5(self):
        """sinusoidal: 0 to scale to -scale to 0"""
        scale_atols = [(1., 0.00001), (10.0, 0.0001)]
        for scale, atol in scale_atols:
            phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
                scale, istep=None,
                animate_scale=True, animate_phase=False, animate_time=False,
                icase_disp=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=2.0, animation_profile='sinusoidal: 0 to scale to -scale to 0',
                fps=5)
            assert len(np.unique(icases_disp)) == 1
            assert len(np.unique(phases)) == 1
            assert np.allclose(analysis_time, 2.0), analysis_time

            msg = ''
            if not np.allclose(scales.min(), -scale, atol=atol):
                msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
            if not np.allclose(scales.max(), scale, atol=atol):
                msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
            if msg:  # pragma: no cover
                raise ValueError(msg)

    def test_animation_bad_1(self):
        """animate_scale/phase/time must be specified"""
        scale = 1.0
        with self.assertRaises(NotImplementedError):
            setup_animation(
                scale, istep=None,
                animate_scale=False, animate_phase=False, animate_time=False,
                icase_disp=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=2.0,  animation_profile='-scale to scale',
                fps=30)

    def test_animation_phase_1(self):
        """phase plot"""
        scale = 1.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=True, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(scales)) == 1
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases[0], 0.), phases
        assert np.allclose(phases[-1], 354.), phases

    def test_animation_time_1(self):
        """time plot"""
        scale = 1.0
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=False, animate_time=True,
            icase_disp=42,
            icase_start=1, icase_end=10, icase_delta=2,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        assert np.array_equal(icases_disp, [1, 3, 5, 7, 9]), icases_disp
        assert len(np.unique(scales)) == 1, scales
        assert len(np.unique(isteps)) > 1, isteps
        assert len(np.unique(phases)) == 1, phases
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases.max(), 0.), phases


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
