"""defines various GUI unit tests"""
import os
import shutil
import unittest
import numpy as np

try:
    import matplotlib.pyplot as plt
    IS_MATPLOTLIB = True
except ImportError:
    IS_MATPLOTLIB = False

import pyNastran
from pyNastran.gui.utils.load_results import (
    load_csv, load_deflection_csv, load_user_geom, create_res_obj)
from pyNastran.gui.utils.version import check_for_newer_version
from pyNastran.gui.utils.utils import find_next_value_in_sorted_list

from pyNastran.gui.menus.legend.write_gif import (
    setup_animation, make_two_sided, make_symmetric, write_gif, IS_IMAGEIO)
from pyNastran.gui.menus.results_sidebar_utils import get_cases_from_tree, build_pruned_tree

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')
from pyNastran.gui.gui_objects.gui_result import GuiResult


class GuiUtils(unittest.TestCase):
    def test_gui_result(self):
        """tests GuiResult"""
        subcase_id = 1
        header = 'cat'
        title = 'dog'
        location = 'node'
        scalar = np.ones(10)
        x = GuiResult(subcase_id, header, title, location, scalar,
                      mask_value=None, nlabels=None, labelsize=None, ncolors=None,
                      colormap='jet', data_map=None, data_format=None, uname='GuiResult')
        y = x + x
        #print(y.scalar)

        y = x + 3
        #print(y.scalar)

        y = x - 3
        #print(y.scalar)

        y = x * 3
        y = x / 3

        y = 2 * x
        y = 2 - x
        y = 2 + x
        #print(y.scalar)

        y = (-x - 1) ** 3
        #print(y.scalar)

        y = (+x + 1) ** 3
        #print(y.scalar)

        scalar = np.ones(8) * 2.
        x2 = GuiResult(subcase_id, header, title, location, scalar,
                      mask_value=None, nlabels=None, labelsize=None, ncolors=None,
                      colormap='jet', data_map=None, data_format=None, uname='GuiResult')
        y2 = 3. / x2
        #print(x2.scalar)
        #print(y2.scalar)

        x2 + y2
        x2 / y2
        x2 * y2
        x2 ** y2
        x2 % 3
        x2 % y2



    def test_check_version_fake(self):
        """
        Tests ``check_for_newer_version``

        we're faking the version for the purpose of the test
        """
        # no dev versions
        version_current_test = '1.0.0'
        version_latest_test = '1.1.0'
        version_latest, version_current, is_newer = check_for_newer_version(
            version_current=version_current_test,
            version_latest=version_latest_test,
            quiet=True)
        #print(version_latest, version_current, is_newer)
        assert version_current == version_current_test
        assert version_latest == version_latest_test
        assert is_newer is True, (version_latest, version_current, is_newer)

        # ------------------------
        # a bigger number version takes priority
        version_current_test = '1.4.0+dev.5378fd363'
        version_latest_test = '1.0.0'
        (version_latest, version_current, is_newer) = check_for_newer_version(
            version_current=version_current_test,
            version_latest=version_latest_test,
            quiet=True)
        assert version_current == version_current_test
        assert version_latest == version_latest_test
        assert is_newer is False, (version_latest, version_current, is_newer)
        # ------------------------
        # a dev version is newer than a non-dev version
        version_current_test = '1.4.0+dev.5378fd363'
        version_latest_test = '1.4.0'
        (version_latest, version_current, is_newer) = check_for_newer_version(
            version_current=version_current_test,
            version_latest=version_latest_test,
            quiet=True)
        assert version_current == version_current_test
        assert version_latest == version_latest_test
        assert is_newer is True, (version_latest, version_current, is_newer)

    def test_check_version_actual(self):
        """tests ``check_for_newer_version`` with actual data"""
        version_latest, version_current, is_newer = check_for_newer_version(
            '1.0.0', quiet=True)
        assert is_newer is True, (version_latest, version_current, is_newer)

        # current/dev release version -> False
        (version_latest, version_current, is_newer) = check_for_newer_version(
            version_current=None,
            version_latest=None,
            quiet=True)
        assert is_newer is False, (version_latest, version_current, is_newer)

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
        A, nids_index, fmt_dict, headers = load_deflection_csv(csv_filename)
        result_type = 'node'

        header0 = headers[0]
        unused_result0 = A[header0]
        #nrows = result0.shape[0]
        #assert nrows == self.nnodes, 'nrows=%s nnodes=%s' % (nrows, self.nnodes)
        header = header0
        islot = 0
        with self.assertRaises(RuntimeError):
            create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                           dim_max=1.0, xyz_cid0=None,
                           is_deflection=False, is_force=False)

        create_res_obj(islot, headers, header, A, fmt_dict, result_type,
                       dim_max=1.0, xyz_cid0=None,
                       is_deflection=True, is_force=False)

    def test_gui_custom_geom_01(self):
        """tests custom_geom.csv"""
        csv_filename = os.path.join(MODEL_PATH, 'custom_geom.csv')
        load_user_geom(csv_filename)

    def test_animation_scale_0(self):
        """0 to scale"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale',
            fps=5)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        expected_scales = [0.0, 0.25*scale, 0.5*scale, 0.75*scale, scale]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_1a(self):
        """0 to scale to 0"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to 0',
            fps=11)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        expected_scales = [0., 0.5*scale, scale]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_1b(self):
        """0 to scale to 0"""
        scale = 1.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to 0',
            fps=10)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), 0.), scales
        assert np.allclose(scales.max(), scale), scales
        expected_scales = [0., 0.5, 1.]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_2(self):
        """-scale to scale"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale',
            fps=5)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        expected_scales = [-scale, -0.5*scale, 0., 0.5*scale, scale]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_3a(self):
        """-scale to scale to -scale"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=11)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        expected_scales = [-scale, 0., scale]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_3b(self):
        """-scale to scale to -scale"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=10)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        assert len(scales) == 3, scales
        expected_scales = [-scale, 0., scale]
        assert_array(scales, expected_scales, 'scales')

    def test_animation_scale_3c(self):
        """-scale to scale to -scale"""
        scale = 2.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='-scale to scale to -scale',
            fps=1)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        expected_scales = [-scale, scale]
        assert_array(scales, expected_scales, 'scales')

        assert np.allclose(analysis_time, 0.5), analysis_time
        assert np.allclose(scales.min(), -scale), scales
        assert np.allclose(scales.max(), scale), scales
        assert len(scales) == 2, scales

    def test_animation_scale_4a(self):
        """0 to scale to -scale to 0"""
        scale_atols = [(1., 0.00000001), (10.0, 0.00000001)]
        for scale, atol in scale_atols:
            out = setup_animation(
                scale, istep=None,
                animate_scale=True, animate_phase=False, animate_time=False,
                icase_disp=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=1.0, animation_profile='0 to scale to -scale to 0',
                fps=5, )
            phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
            assert len(np.unique(icases_disp)) == 1
            assert len(np.unique(phases)) == 1
            assert np.allclose(analysis_time, 1.0), analysis_time
            expected_scales = [0., scale, 0., -scale]
            assert_array(scales, expected_scales, 'scales')

            msg = ''
            if not np.allclose(scales.min(), -scale, atol=atol):  # pragma: no cover
                msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
            if not np.allclose(scales.max(), scale, atol=atol):  # pragma: no cover
                msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
            if msg:  # pragma: no cover
                raise ValueError(msg)

    def test_animation_scale_4b(self):
        """0 to scale to -scale to 0"""
        scale, atol = (1., 0.00000001)
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=1.0, animation_profile='0 to scale to -scale to 0',
            fps=7, animate_in_gui=True)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 1.0), analysis_time

        expected_scales = [0., 0.5*scale, scale, 0.5*scale,
                           0., -0.5*scale, -scale, -0.5*scale]
        assert_array(scales, expected_scales, 'scales')

        scales, phases, icases_fringe, icases_disp, icases_vector, isteps = make_two_sided(
            scales, phases, icases_fringe, icases_disp, icases_vector, isteps, onesided)
        expected_scales = [0., 0.5*scale, scale, 0.5*scale, 0., -0.5*scale, -scale, -0.5*scale,
                           0., 0.5*scale, scale, 0.5*scale, 0., -0.5*scale, -scale, -0.5*scale, ]
        assert_array(scales, expected_scales, 'scales')

        msg = ''
        if not np.allclose(scales.min(), -scale, atol=atol):  # pragma: no cover
            msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
        if not np.allclose(scales.max(), scale, atol=atol):  # pragma: no cover
            msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
        if msg:  # pragma: no cover
            raise ValueError(msg)

    def test_animation_scale_4c(self):
        """0 to scale to -scale to 0"""
        scale, atol = (1., 0.00000001)
        out = setup_animation(
            scale, istep=None,
            animate_scale=True, animate_phase=False, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, animation_profile='0 to scale to -scale to 0',
            fps=7, animate_in_gui=True)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(phases)) == 1
        assert np.allclose(analysis_time, 2.0), analysis_time
        expected_scales = [0., 0.25, 0.5, 0.75, 1., 0.75, 0.5, 0.25, 0.,
                           -0.25, -0.5, -0.75, -1., -0.75, -0.5, -0.25]
        assert_array(scales, expected_scales, 'scales')

    #msg = '%s = %s\n' % (actual_name, actual_array)
    #msg += '%s = %s\n' % (expected_name, expected_array)

        msg = ''
        if not np.allclose(scales.min(), -scale, atol=atol):  # pragma: no cover
            msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
        if not np.allclose(scales.max(), scale, atol=atol):  # pragma: no cover
            msg += 'scales=%s max=%s expected=%s' % (scales, scales.max(), scale)
        if msg:  # pragma: no cover
            raise ValueError(msg)

    def test_animation_scale_5(self):
        """sinusoidal: 0 to scale to -scale to 0"""
        scale_atols = [(1., 0.00001), (10.0, 0.0001)]
        for scale, atol in scale_atols:
            out = setup_animation(
                scale, istep=None,
                animate_scale=True, animate_phase=False, animate_time=False,
                icase_disp=42,
                icase_start=None, icase_end=None, icase_delta=None,
                time=2.0, animation_profile='sinusoidal: 0 to scale to -scale to 0',
                fps=5)
            phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
            assert len(np.unique(icases_disp)) == 1
            assert len(np.unique(phases)) == 1
            assert np.allclose(analysis_time, 2.0), analysis_time

            msg = ''
            if not np.allclose(scales.min(), -scale, atol=atol):  # pragma: no cover
                msg += 'scales=%s min=%s expected=%s\n' % (scales, scales.min(), -scale)
            if not np.allclose(scales.max(), scale, atol=atol):  # pragma: no cover
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
                time=2.0, animation_profile='-scale to scale',
                fps=30)

    def test_animation_phase_1(self):
        """phase plot"""
        scale = 1.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=True, animate_time=False,
            icase_disp=42,
            icase_start=None, icase_end=None, icase_delta=None,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert len(np.unique(icases_disp)) == 1
        assert len(np.unique(scales)) == 1
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases[0], 0.), phases
        assert np.allclose(phases[-1], 354.), phases

    def test_animation_time_1(self):
        """time plot"""
        scale = 1.0
        out = setup_animation(
            scale, istep=None,
            animate_scale=False, animate_phase=False, animate_time=True,
            icase_disp=42,
            icase_start=1, icase_end=10, icase_delta=2,
            time=2.0, #animation_profile='0 to scale',
            fps=30)
        phases, icases_fringe, icases_disp, icases_vector, isteps, scales, analysis_time, onesided, endpoint = out
        assert np.array_equal(icases_disp, [1, 3, 5, 7, 9]), icases_disp
        assert len(np.unique(scales)) == 1, scales
        assert len(np.unique(isteps)) > 1, isteps
        assert len(np.unique(phases)) == 1, phases
        assert np.allclose(analysis_time, 2.0), analysis_time
        assert np.allclose(phases.max(), 0.), phases

    def test_cases_from_tree(self):
        """tests ``get_cases_from_tree``"""
        form = [
            [u'Geometry', None, [
                (u'NodeID', 0, []),
                (u'ElementID', 1, []),
                (u'PropertyID', 2, []),
                (u'MaterialID', 3, []),
                (u'E', 4, []),
                (u'Element Checks', None, [
                    (u'ElementDim', 5, []),
                    (u'Min Edge Length', 6, []),
                    (u'Min Interior Angle', 7, []),
                    (u'Max Interior Angle', 8, [])],
                ),],
            ],
        ]
        cases = get_cases_from_tree(form)
        assert np.array_equal(cases, [0, 1, 2, 3, 4, 5, 6, 7, 8]), cases

    def test_cases_from_tree_spaces(self):
        """tests ``get_cases_from_tree``"""
        form = [
            [u'Geometry', '', [
                (u'NodeID', 0, []),
                (u'ElementID', 1, []),
                (u'PropertyID', 2, []),
                (u'MaterialID', 3, []),
                (u'E', 4, []),
                (u'Element Checks', '', [
                    (u'ElementDim', 5, []),
                    (u'Min Edge Length', 6, []),
                    (u'Min Interior Angle', 7, []),
                    (u'Max Interior Angle', 8, [])],
                ),],
            ],
        ]
        cases = get_cases_from_tree(form)
        assert np.array_equal(cases, [0, 1, 2, 3, 4, 5, 6, 7, 8]), cases

    def test_build_pruned_tree(self):
        """tests ``get_cases_from_tree``"""
        tree = [
            [u'Geometry', None, [
                (u'NodeID', 0, []),
                (u'ElementID', 1, []),
                (u'PropertyID', 2, []),
                (u'MaterialID', 3, []),
                (u'E', 4, []),
                (u'Element Checks', None, [
                    (u'ElementDim', 5, []),
                    (u'Min Edge Length', 6, []),
                    (u'Min Interior Angle', 7, []),
                    (u'Max Interior Angle', 8, [])],
                ),],
            ],
        ]
        cases = [0, 1, 5]
        tree2 = build_pruned_tree(tree, cases)
        assert tree2 == [(u'Geometry', None, [(u'NodeID', 0, []), (u'ElementID', 1, []), (u'Element Checks', None, [(u'ElementDim', 5, [])])])], tree2

        cases = [0, 1, 7]
        tree2 = build_pruned_tree(tree, cases)
        assert tree2 == [(u'Geometry', None, [(u'NodeID', 0, []), (u'ElementID', 1, []), (u'Element Checks', None, [(u'Min Interior Angle', 7, [])])])], tree2

    def test_write_gif(self):
        """tests the write_gif method"""
        png_filenames = []
        gif_filename = os.path.join('mydir', 'line.gif')
        png_filenames = ['line.png', 'line.png']
        if os.path.exists('line.png'):
            os.remove('line.png')
        success = write_gif(gif_filename, png_filenames, time=2.0,
                            onesided=True, nrepeat=0,
                            delete_images=False, make_gif=False)
        if IS_IMAGEIO:
            self.assertTrue(success)
        else:
            self.assertFalse(success)

        if IS_MATPLOTLIB and IS_IMAGEIO:
            fig = plt.figure(1)
            ax = fig.gca()
            ax.clear()
            xrange = np.linspace(0., 1)
            ax.plot(xrange)
            fig.savefig('line.png')
            png_filenames = ['line.png', 'line.png']

            gif_filename = os.path.join('mydir', 'line.gif')
            success = write_gif(gif_filename, png_filenames, time=2.0,
                                onesided=True, nrepeat=0,
                                delete_images=True, make_gif=True)
            self.assertTrue(success)
            shutil.rmtree('mydir')


def assert_array(actual_array, expected_array, name):
    """checks the arrays for equivalence and the correct length"""
    expected_array = np.asarray(expected_array)
    actual_name = 'actual_%s  ' % (name)
    expected_name = 'expected_%s' % (name)
    if not len(actual_array) == len(expected_array):
        raise AssertionError('Invalid Length\n'
                             '%s=%s; len=%s\n'
                             '%s=%s; len=%s' % (
                                 actual_name, actual_array.tolist(), len(actual_array),
                                 expected_name, expected_array.tolist(), len(expected_array)))

    if not np.allclose(actual_array, expected_array):
        raise AssertionError('%s=scales=%s; len=%s\nexpected=%s; len=%s' % (
            actual_name, actual_array.tolist(), len(actual_array),
            expected_name, expected_array.tolist(), len(expected_array)))
    #np.testing.assert_almost_equal(scales, expected_scales)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
