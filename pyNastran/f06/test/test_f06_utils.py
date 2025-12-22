"""
tests:
 - plot sol_145
"""
import os
import warnings
from pathlib import Path
import unittest

import numpy as np
from scipy.spatial import KDTree
from scipy.interpolate import interp1d, splrep, splev

from cpylog import SimpleLogger, get_logger
try:
    import matplotlib  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    #matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt

    #try:  # pragma: no cover
        #plt.figure()
        #plt.close()
    #except Exception:  # pragma: no cover
    plt.switch_backend('Agg')

import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh

from pyNastran.f06.utils import (
    split_float_colons, split_int_colon,
    cmd_line_plot_flutter, cmd_line as cmd_line_f06)
from pyNastran.f06.parse_flutter import (
    plot_flutter_f06, make_flutter_plots, make_flutter_response,
    FlutterResponse,
)
from pyNastran.f06.flutter_response import _reshape_eigenvectors
from pyNastran.f06.parse_trim import read_f06_trim
from pyNastran.f06.f06_matrix_parser import _parse_real_row_lines, _parse_complex_row_lines
from pyNastran.f06.f06_to_pressure_loads import f06_to_pressure_loads
from pyNastran.op2.op2 import OP2

IS_DEV = pyNastran.DEV
if IS_DEV: # pragma: no cover
    from pyNastran.f06.dev.read_sol_200 import plot_sol_200  # read_sol_200

DIRNAME = os.path.dirname(__file__)
PKG_PATH = Path(pyNastran.__path__[0])
MODEL_PATH = PKG_PATH / '..' / 'models'
AERO_PATH = MODEL_PATH / 'aero'

#pyNastran\bdf\cards\aero\examples\flutter\case6
AERO_EXAMPLES = PKG_PATH / 'bdf' / 'cards' / 'aero' / 'examples'
assert AERO_EXAMPLES.exists(), print_bad_path(AERO_EXAMPLES)


class TestF06Flutter(unittest.TestCase):

    def test_reshape_eigenvectors(self):
        """
        tests reshape_eigenvectors, but hacks the
        standard input in order to confirm the dimensions
        and organization are correct. This is done using
        an extra row in the input matrix.
        """
        nmode = 2
        nvel = 1

        # 2x2 -> 2x2x1
        #(imode, jmode, ivel)
        eigenvectors = np.array([
            [1.+0.j, -0.0223484-0.00115277j],
            [-0.00381052-0.00061724j, 1.+0.j],
            [0.1, 0.2],   # fake row
        ])
        assert eigenvectors.shape == (nmode+1, nmode*nvel), eigenvectors.shape

        # (ivel, imode, mpfs)
        ivel = 0
        eigenvectors_expected = np.zeros((nvel, nmode, nmode+1), dtype=eigenvectors.dtype)
        eigenvectors_expected[ivel, 0, :] = eigenvectors[:, 0]
        eigenvectors_expected[ivel, 1, :] = eigenvectors[:, 1]

        #(nmode*nvel*nvel) -> (nmode, nvel)?
        eigr_eigi_vel = np.array([
            [-9.88553e-02, 1.71977e+01, 1.52383e+02],
            [-1.71903e-01, 6.60547e+01, 1.52383e+02]])
        eigr_eigi_vel_expected = np.zeros((nvel, nmode, 3))
        assert eigr_eigi_vel.shape == (nmode*nvel, 3), eigr_eigi_vel.shape

        eigenvectors2, eigr_eigi_vel2 = _reshape_eigenvectors(
            eigenvectors, eigr_eigi_vel, incorrect_shape=True)
        assert eigenvectors2.dtype == eigenvectors_expected.dtype, (eigenvectors2.dtype, eigenvectors_expected.dtype)
        assert eigenvectors2.shape == eigenvectors_expected.shape, (eigenvectors2.shape, eigenvectors_expected.shape)
        assert np.allclose(eigenvectors2, eigenvectors_expected)

        assert eigr_eigi_vel2.shape == eigr_eigi_vel_expected.shape, (eigr_eigi_vel2.dtype, eigr_eigi_vel_expected.dtype)
        assert eigr_eigi_vel2.shape == eigr_eigi_vel_expected.shape, (eigr_eigi_vel2.shape, eigr_eigi_vel_expected.shape)

    def test_reshape_eigenvectors2(self):
        """
        we use a faked number of modes (modes+1) to
        make sure we're swapping the correct axes
        """
        # 2 modes, 4 vel
        ivel = 0
        nmode = 2
        nvel = 4
        eigenvectors = np.array([
            [1, 2, 3, 4, 5, 6, 7, 8],
            [10, 20, 30, 40, 50, 60, 70, 80],
        ])
        eigr_eigi_vel = np.array([
            [1, 10, 100],
            [2, 20, 200],
            [3, 30, 300],
            [4, 40, 400],
            [5, 50, 500],
            [6, 60, 600],
            [7, 70, 700],
            [8, 80, 800],
        ])
        eigenvectors_expected = np.zeros((nvel, nmode, nmode), dtype=eigenvectors.dtype)
        eigenvectors_expected[ivel+0, 0, :] = eigenvectors[:, 0]
        eigenvectors_expected[ivel+0, 1, :] = eigenvectors[:, 1]
        eigenvectors_expected[ivel+1, 0, :] = eigenvectors[:, 2]
        eigenvectors_expected[ivel+1, 1, :] = eigenvectors[:, 3]
        eigenvectors_expected[ivel+2, 0, :] = eigenvectors[:, 4]
        eigenvectors_expected[ivel+2, 1, :] = eigenvectors[:, 5]
        eigenvectors_expected[ivel+3, 0, :] = eigenvectors[:, 6]
        eigenvectors_expected[ivel+3, 1, :] = eigenvectors[:, 7]

        eigr_eigi_vel_expected = np.zeros((nvel, nmode, 3))
        eigr_eigi_vel_expected[ivel+0, 0, :] = eigr_eigi_vel[0, :]
        eigr_eigi_vel_expected[ivel+0, 1, :] = eigr_eigi_vel[1, :]
        eigr_eigi_vel_expected[ivel+1, 0, :] = eigr_eigi_vel[2, :]
        eigr_eigi_vel_expected[ivel+1, 1, :] = eigr_eigi_vel[3, :]
        eigr_eigi_vel_expected[ivel+2, 0, :] = eigr_eigi_vel[4, :]
        eigr_eigi_vel_expected[ivel+2, 1, :] = eigr_eigi_vel[5, :]
        eigr_eigi_vel_expected[ivel+3, 0, :] = eigr_eigi_vel[6, :]
        eigr_eigi_vel_expected[ivel+3, 1, :] = eigr_eigi_vel[7, :]
        eigenvectors2, eigr_eigi_vel2 = _reshape_eigenvectors(
            eigenvectors, eigr_eigi_vel)
        assert np.allclose(eigr_eigi_vel2.shape, eigr_eigi_vel_expected.shape)
        assert np.allclose(eigenvectors2, eigenvectors_expected)
        assert np.allclose(eigr_eigi_vel2, eigr_eigi_vel_expected)

    def test_make_grid_point_singularity_table(self):
        model = OP2()
        failed = [
            (1, 1), (1, 2), (1, 7),
            (2, 1), (2, 2), (2, 7),
        ]
        table = model.make_grid_point_singularity_table(failed)
        msg = '\n'.join(table.split('\n')[:-2]) + '\n'
        expected = (
            '0                                         G R I D   P O I N T   S I N G U L A R I T Y   T A B L E\n'
            '0                             POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET\n'
            '                               ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION\n'
            '                                1        G      1         0.00E+00          B        F         SB       SB   *\n'
            '                                1        G      2         0.00E+00          B        F         SB       SB   *\n'
            '                                1        G      7         0.00E+00          B        F         SB       SB   *\n'
            '                                2        G      1         0.00E+00          B        F         SB       SB   *\n'
            '                                2        G      2         0.00E+00          B        F         SB       SB   *\n'
            '                                2        G      7         0.00E+00          B        F         SB       SB   *\n'
            #'1                                                                           JANUARY   1, 2000  pyNastran v1.5.0+dev.140b9348d  PAGE     1\n'
        )
        #print(msg)
        assert msg == expected

        failed = []
        msg = model.make_grid_point_singularity_table(failed)
        assert msg == '', msg
        #print(msg)

    def test_f06_pt145(self):
        """tests read_f06_trim"""
        log = get_logger(log=None, level=None, encoding='utf-8')
        f06_filename = AERO_PATH / 'pt145.f06'

        plot_methods = [
            '--eas', '--tas', '--rho', '--alt',
            '--mach', '--q', '--index',
        ]
        for plot_method in plot_methods:
            argv = ['f06', 'plot_145', str(f06_filename),
                    plot_method]
            cmd_line_plot_flutter(argv=argv, plot=False,
                                  show=False, log=log)

        argv = [
            'f06', 'plot_145', str(f06_filename),
            '--tas',
            '--in_units', 'si',
            '--out_units', 'english_in',
            '--freq_tol', '0.02',
            '--freq_tol_remove', '0.01',
            '--vd_limit', '100.',
            '--ylimfreq', '0:',
            '--damping_limit', '0.001',
        ]
        cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB,
                              show=False, log=log)

    def test_plot_flutter_bah(self):
        """tests plot_flutter_f06"""
        f06_filename = AERO_PATH / 'bah_plane' / 'bah_plane.f06'
        log = get_logger(log=None, level=None, encoding='utf-8')
        flutters, data = plot_flutter_f06(
            f06_filename, show=False, close=True,
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True, plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB,
            log=log)
        matrices = data['matrices']
        matrices['freq']
        matrices['KHH']
        matrices['BHH']
        matrices['MHH']
        assert len(flutters) == 2, list(flutters.keys())
        assert flutters[1].results.shape == (10, 22, 7), flutters[1].results.shape
        assert flutters[2].results.shape == (10, 30, 7), flutters[2].results.shape
        fix_modes_2024(flutters[1])
        fix_modes_2024(flutters[2])
        # plot_flutter_f06(f06_filename, show=True, close=False,
        #                 plot_vg=False, plot_vg_vf=True, plot_root_locus=True, plot_kfreq_damping=False,
        #                 subcases=[1, 3],
        #                 log=log)

    def test_plot_flutter_0012(self):
        """
        tests plot_flutter_f06

        has issues with writing the subcase...
        """
        dirname = AERO_PATH / '2_mode_flutter'
        f06_filename = dirname / '0012_flutter.f06'

        #log = get_logger(log=None, level=None, encoding='utf-8')
        # log = get_logger(log=None, level=False, encoding='utf-8')
        log = SimpleLogger(level='warning')
        plot_flutter_f06(
            f06_filename,
            modes=[2],
            plot_type='alt',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, log=log,
            close=True,
        )

        flutters, data = plot_flutter_f06(
            f06_filename,
            f06_units='si', out_units=None,
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            export_csv_filename=dirname/'nastran.csv',
            export_f06_filename=dirname/'nastran.f06',
            export_veas_filename=dirname/'nastran.veas',
            export_zaero_filename=dirname/'zaero.f06',
            vg_filename='vg_subcase_%i.png',
            vg_vf_filename='vg_vf_subcase_%i.png',
            kfreq_damping_filename='kfreq_damping_subcase_%i.png',
            root_locus_filename='root_locus_subcase_%i.png',
            modal_participation_filename='modal_participation_subcase_%i.png',
            plot=IS_MATPLOTLIB, show=False, log=log,
            close=True,
        )
        matrices = data['matrices']
        matrices['freq']
        matrices['MHH']
        matrices['BHH']
        matrices['KHH']
        flutter = flutters[1]
        assert flutter.results.shape == (2, 93, 12), flutter.results.shape
        fix_modes_2024(flutter)

        if IS_MATPLOTLIB:
            os.remove('vg_subcase_1.png')
            os.remove('vg_vf_subcase_1.png')
            os.remove('kfreq_damping_subcase_1.png')
            os.remove('root_locus_subcase_1.png')

            os.remove(dirname/'nastran.csv')
            os.remove(dirname/'nastran.f06')
            os.remove(dirname/'nastran.veas')
            os.remove(dirname/'zaero.f06')

        with self.assertRaises(NotImplementedError):
            plot_flutter_f06(
                f06_filename,
                plot_type='tas',
                f06_units='cat', out_units=None,
                plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
                plot_kfreq_damping=True,
                plot=IS_MATPLOTLIB, show=False, close=True, log=log)
        with self.assertRaises(NotImplementedError):
            plot_flutter_f06(
                f06_filename,
                plot_type='density',
                f06_units='si', out_units='english_ft',
                plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
                plot_kfreq_damping=True,
                show=True, close=True, log=log)

        flutters, mass = plot_flutter_f06(
            f06_filename,
            plot_type='rho',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, close=True, log=log)
        #flutters[1].plot_zimmerman([1, 2], show=True)
        flutters[1].plot_zimmerman([1, 2], show=False)

        plot_flutter_f06(
            f06_filename,
            plot_type='freq',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, clear=True, close=True, log=log)

        plot_flutter_f06(
            f06_filename,
            plot_type='kfreq',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, clear=True, close=True, log=log)

        plot_flutter_f06(
            f06_filename,
            plot_type='ikfreq',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, clear=True, close=True, log=log)

        plot_flutter_f06(
            f06_filename,
            modes=[2],
            plot_type='damp',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, clear=True, close=True, log=log)

        plot_type = 'eas'
        modes = None
        # xlim = [0., 500000.] # in/s
        xlim = None
        ylim_damping = [-0.2, 0.2]
        ylim_freq = None
        ylim_kfreq = None
        if IS_MATPLOTLIB:
            make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                               plot_type,
                               plot_vg=True, plot_vg_vf=True,
                               plot_root_locus=True, plot_kfreq_damping=True,
                               nopoints=True, noline=False,
                               # export_zaero_filename=export_zaero_filename,
                               # export_veas_filename=export_veas_filename,
                               # export_f06_filename=export_f06_filename,
                               # vg_filename=vg_filename,
                               # vg_vf_filename=vg_vf_filename,
                               # root_locus_filename=root_locus_filename,
                               # kfreq_damping_filename=kfreq_damping_filename,
                               show=False, clear=True, close=True)

            export_csv_filename = os.path.join(DIRNAME, 'csv_%d.csv')
            export_zaero_filename = os.path.join(DIRNAME, 'zaero_%d.f06')
            export_veas_filename = os.path.join(DIRNAME, 'flutter_%d.veas')
            export_f06_filename = os.path.join(DIRNAME, 'flutter_%d.f06')
            make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                               plot_type,
                               plot_vg=True, plot_vg_vf=True,
                               plot_root_locus=True, plot_kfreq_damping=True,
                               nopoints=True, noline=False,
                               export_csv_filename=export_csv_filename,
                               export_zaero_filename=export_zaero_filename,
                               export_veas_filename=export_veas_filename,
                               export_f06_filename=export_f06_filename,
                               # vg_filename=vg_filename,
                               # vg_vf_filename=vg_vf_filename,
                               # root_locus_filename=root_locus_filename,
                               # kfreq_damping_filename=kfreq_damping_filename,
                               show=False, clear=True, close=True)

    def test_cmd_line_plot_flutter_0012(self):
        log = get_logger(log=None, level=None, encoding='utf-8')
        f06_filename = AERO_PATH / '2_mode_flutter' / '0012_flutter.f06'
        ivel = '0'
        mode = '1:2'
        argv = [
            'f06', 'plot_145', str(f06_filename), '--eas',
            '--in_units', 'si', '--out_units', 'english_in',
            '--modal', ivel, mode,
            '--mag_tol', '0.1',
            '--modes', '1:', '--ylimdamp', '-.3:', '--export_csv',
            '--ncol', '2',
        ]
        cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)

    def test_cmd_line_plot_flutter_no_input_0012(self):
        """no input???"""
        log = get_logger(log=None, level=None, encoding='utf-8')
        f06_filename = AERO_PATH / '2_mode_flutter' / '0012_flutter.f06'
        argv = ['f06', 'plot_145', str(f06_filename), '--eas',
                '--out_units', 'english_in']
        flutters = cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)

    @unittest.skipIf(not IS_MATPLOTLIB, 'no matplotlib')
    def test_plot_func_0012(self):
        log = SimpleLogger(level='warning')
        dirname = AERO_PATH / '2_mode_flutter'
        f06_filename = dirname / '0012_flutter.f06'
        zaero_filename = dirname / 'junk.zaero'
        flutters, data = make_flutter_response(f06_filename, f06_units='si', log=log)
        flutter = flutters[1]
        flutter.set_plot_settings(
            #figsize=(4,4),
            xtick_major_locator_multiple=[50.0, 50.0],
            ytick_major_locator_multiple=[0.05, None],
        )
        flutter.set_plot_options(noline=False)
        flutter.set_symbol_settings(
            nopoints=False,
            show_mode_number=False,
            point_spacing=3,
            markersize=None,
            #markersize=0,
        )
        flutter.plot_vg_vf(plot_type='eas')
        flutter.plot_vg(plot_type='eas')
        #---------------------------------
        flutter.set_plot_settings(figsize=(4, 4),)
        flutter.set_symbol_settings(
            nopoints=False,
            show_mode_number=True,
            point_spacing=3,
            markersize=None,
        )
        flutter.plot_vg_vf(plot_type='eas')
        flutter.plot_vg_vf(plot_type='eas', show_detailed_mode_info=True)
        flutter.plot_vg_vf(plot_type='eas', mode_switch_method='freq')
        flutter.plot_vg_vf(plot_type='eas', mode_switch_method='damping')
        with self.assertRaises(RuntimeError):
            flutter.plot_vg_vf(plot_type='eas', mode_switch_method='cat')
        flutter.plot_vg(plot_type='eas')
        flutter.export_to_zaero(zaero_filename)
        str(flutter.object_methods())

        flutter.set_symbol_settings(
            nopoints=False,
            show_mode_number=False,
            point_spacing=3,
            markersize=8,
        )
        flutter.plot_vg_vf(plot_type='eas')
        flutter.plot_vg(plot_type='eas')
        flutter.plot_kfreq_damping()
        flutter.plot_kfreq_damping2()

    def test_fix_modes_0012(self):
        log = SimpleLogger(level='warning')
        f06_filename = AERO_PATH / '2_mode_flutter' / '0012_flutter.f06'
        flutters, mass = make_flutter_response(f06_filename, f06_units='si', log=log)
        flutter = flutters[1]
        fix_modes_2024(flutter)

    def test_fix_modes_constant(self):
        """constant to check dxyz"""
        # nmodes, nvel
        eigr = np.array([
            [0., 1., 2.],
            [0., 1., 2.],
        ])
        eigi = np.array([
            [1., 1., 1.],
            [10., 10., 10.],
        ])
        nmodes, nvel = eigr.shape
        out = _fix_modes(eigr, eigi, nmodes, nvel, kmodes=0)

    def test_fix_modes_linear(self):
        """constant to check dxyz"""
        # nmodes, nvel
        eigr = np.array([
            [0., 1., 2.],
            [0., 1., 2.],
        ])
        eigi = np.array([
            [1., 1., 1.],
            [0., .5, 1.1],
        ])
        nmodes, nvel = eigr.shape
        out = _fix_modes(eigr, eigi, nmodes, nvel, kmodes=0)

        def fix_modes_dumb(self):
            """
            steps
            1. sort by frequency
            """
            coeffs = [0, -40, -2, 0.5]
            x = np.linspace(-10, 13., num=100)
            y1 = coeffs[0] + coeffs[1] * x + coeffs[2] * x**2 + coeffs[3] * x**3

            coeffs2 = [-20, 5, 1, 0.1]
            y2 = coeffs[0] + coeffs[1] * x + coeffs[2] * x ** 2 + coeffs[3] * x ** 3

            fig = plt.Figure(1)
            ax = fig.gca()
            ax.plot(x, y1, label='y1')
            ax.plot(x, y2, label='y2')
            plt.grid(True)
            #plt.show()
            #flutter.sort_modes_by_freq(freq)


def fix_modes_2024(flutter: FlutterResponse,
                   kmodes: int=0,
                   debug: bool=False) -> None:
    ieigr = flutter.ieigr
    ieigi = flutter.ieigi
    eigr = flutter.results[:,]

    # results[imode, ivelocity, iresponse]
    eigr = flutter.results[:, :, ieigr]
    eigi = flutter.results[:, :, ieigi]
    nmodes, nvel, nresponses = flutter.results.shape

    out = _fix_modes(eigr, eigi, nmodes, nvel, kmodes=kmodes, debug=debug)
    outi_expected = np.arange(nmodes)
    out_expected = np.repeat(outi_expected, nvel).reshape(nmodes, nvel)
    if np.array_equal(out, out_expected):
        if debug:
            print('no mode switching')
    else:
        raise RuntimeError('mode_switching')
    return


def _fix_modes(eigr: np.ndarray,
               eigi: np.ndarray,
               nmodes: int, nvel: int,
               kmodes: int=0,
               debug: bool=False) -> np.ndarray:
    if kmodes == 0:
        kmodes = max(2, nmodes // 2)
    kmodes = min(kmodes, nmodes)
    # scipy.interpolate.interp1d(eigr, eigi)
    # scipy.interpolate.interp2d

    # all_data = [nmode, nvel, 2]
    all_data = np.dstack([eigr, eigi])
    data0 = np.column_stack([eigr[:, 0], eigi[:, 0]])
    assert np.array_equal(data0, all_data[:, 0, :])
    assert all_data.shape == (nmodes, nvel, 2), all_data.shape
    imode_expected = np.arange(nmodes, dtype='int32')

    # floats because ????
    ivelocity = np.arange(nvel, dtype='float64')

    data0 = all_data[:, 0, :]
    out = np.full((nmodes, nvel), -1, dtype='int32')
    out[:, 0] = np.arange(nmodes)
    tree0 = KDTree(data0)
    for ivel in range(1, nvel):
        # datai = [imode, 2=real/imag]
        datai = np.column_stack([eigr[:, ivel], eigi[:, ivel]])
        dxyz, ixyz = tree0.query(datai, k=kmodes)
        if debug:  # pragma: no cover
            print(f'\ndata0[{ivel}]:\n{data0}')
            print(f'datai[{ivel}]:\n{datai}')
            print(f'dxyz[{ivel}]:\n{dxyz}')
            print(f'ixyz[{ivel}]:\n{ixyz}')

        imode_next = ixyz[:, 0]
        iunique, idx, inv, counts = np.unique(
            imode_next,
            return_index=True,
            return_inverse=True,
            return_counts=True)

        assert counts.min() > 0, counts  # no points...
        isort = np.full(nmodes, -1)
        i1 = np.where(counts == 1)[0]
        i2 = np.where(counts > 1)[0]

        if len(i1):
            # fill the points with no mode switching
            #idx0_i1 = idx[i1]
            idx_i1 = inv[i1]
            try:
                iunique_idx = iunique[idx_i1]
                isort[idx_i1] = iunique[idx_i1]
            except IndexError:
                bad_index

        if len(i2):
            if debug:
                print(f'*ivel={ivel}: switching; i2={i2}')

            # multiple modes found the same closest point
            if ivel == 1:
                # first point
                raise NotImplementedError('not handling 2nd point')
            else:
                nx = ivel
                x = ivelocity[:ivel]
                assert len(x) == nx
                xi = ivel
                s = min(ivel-1, 5)
                eigri_predicted = np.empty((nmodes, 2), dtype='float64')

                # pack up the next predicted_eigris for all the modes
                # for imodei, imode in enumerate(i2):
                for imode in range(nmodes):
                    # all_data = [nmode, nvel, 2]
                    eigr_mode = all_data[imode, :ivel, 0]
                    eigi_mode = all_data[imode, :ivel, 1]
                    assert len(eigr_mode) == nx
                    if nx == 2:
                        # y2 = (y1-y0)/(x1-x0) * (x2-x0) + y0
                        #    = (y1-y0) * (x2-x0)/(x1-x0) + y0
                        #    = (y1-y0) * m + y0
                        #    = y1*m + y0 - y0*m
                        #    = y1*m + y0 * (1-m)
                        #
                        # let x[i] = [0, 1, 2, 3, ...]
                        # so m = (2-0)/(1-0) = (3-1)/(2-1) = 2
                        m = 2.0
                        yr0 = eigr_mode[-2]
                        yr1 = eigr_mode[-1]
                        eigrii = eigr_mode[-1]*m + eigr_mode[-2]*(1-m)
                        yi0 = eigi_mode[-2]
                        yi1 = eigi_mode[-1]
                        eigiii = eigi_mode[-1]*m + eigi_mode[-2]*(1-m)
                        #func_eigr = interp1d(x, eigr_mode, fill_value="extrapolate")
                        #func_eigi = interp1d(x, eigi_mode, fill_value="extrapolate")
                        #eigrii = func_eigr(xi)
                        #eigiii = func_eigi(xi)
                        if debug:
                            print(f'  eig_predicted[{imode+1}] = {eigrii} + {eigiii}j')
                        x = 1
                    else:
                        spline_eigr = splrep(x, eigr_mode, s=s)
                        spline_eigi = splrep(x, eigi_mode, s=s)
                        eigrii = splev(xi, spline_eigr)
                        eigiii = splev(xi, spline_eigi)
                    # eigri_predicted[imodei, :] = [eigrii, eigiii]
                    eigri_predicted[imode, :] = [eigrii, eigiii]

                # using the predicted_eigris, find the closest point
                etree = KDTree(eigri_predicted)
                adxyz, aixyz = tree0.query(datai, k=kmodes)
                aimode_next = aixyz[:, 0]

                isort_predicted = np.full(nmodes, -1)
                for imodei, imode in enumerate(aimode_next):
                    # if the predicted mode = the actual mode,
                    # default to not changing it (e.g., constant)
                    if imodei == imode:
                        isort_predicted[imodei] = imode

                imissing = np.where(isort_predicted == -1)[0]
                nmissing = len(imissing)
                if nmissing:
                    tree1 = KDTree(eigri_predicted[imissing, :])
                    # datai = [imode, 2=real/imag]
                    kmodesb = min(kmodes, nmissing)
                    bdxyz, bixyz = tree1.query(datai[imissing, :], k=kmodesb)
                    if nmissing == 1:
                        isort_predicted = imode_expected
                    else:
                        nmissing2
                        for imodei, imiss in enumerate(imissing):
                            asdf
                        isort_predicted
                else:
                    if 0:
                        aiunique, aidx, ainv, acounts = np.unique(
                            aimode_next,
                            return_index=True,
                            return_inverse=True,
                            return_counts=True)
                        assert acounts.min() > 0, acounts  # no points...
                        a1 = np.where(acounts == 1)[0]
                        a2 = np.where(acounts != 1)[0]
                        assert len(a1) == len(i2)
                        assert len(a2) == 0, a2
                        bbb
            isort = isort_predicted
        elif debug:
            print(f'ivel={ivel}: no switching')

        out[:, ivel] = isort
        imode_next
        data0 = all_data[:, ivel, :]
        tree0 = KDTree(data0)
    return out


class TestZaeroFlutter(unittest.TestCase):
    def test_zaero_gafa(self):
        from pyNastran.f06.dev.flutter.read_zaero_out import read_zaero_out
        #bdf_filename = AERO_EXAMPLES / 'flutter' / 'case6' / 'agardztran.bdf'
        f06_filename = AERO_EXAMPLES / 'flutter' / 'case6' / 'agardztran.out'
        png_filename = AERO_EXAMPLES / 'flutter' / 'case6' / 'agardztran.png'
        responses, mass = read_zaero_out(f06_filename, debug=False)
        # plot_sol_200(f06_filename, png_filename=png_filename,
        #              show=True)
        if len(responses) == 1:
            warnings.warn('nresponses=1 and should be 2')
        else:
            assert len(responses) == 2, list(responses.keys())


class TestF06Utils(unittest.TestCase):
    @unittest.skipIf(not IS_DEV, 'skipping plot_sol_200')
    def test_opt_aerobeam(self):
        """tests optimization"""
        log = SimpleLogger('warning')
        f06_filename = AERO_PATH / 'aerobeam.f06'
        png_filename = AERO_PATH / 'aerobeam.png'
        plot_sol_200(f06_filename, png_filename=png_filename,
                     show=True, log=log)
        argv = ['f06', 'plot_200', str(f06_filename)]
        cmd_line_f06(argv, plot=False, show=False, log=log)

    @unittest.skipIf(not IS_DEV, 'skipping plot_sol_200')
    def test_opt_mdb200(self):
        """tests optimization"""
        f06_filename = MODEL_PATH / 'other' / 'mdb200.f06'
        png_filename = MODEL_PATH / 'other' / 'mdb200.png'
        log = SimpleLogger(level='warning')
        plot_sol_200(
            f06_filename, png_filename=png_filename,
            log=log, show=True)

    def test_f06_trim_bwb(self):
        """tests f06_to_pressure_loads"""
        bdf_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.bdf'
        aerobox_caero_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.caero.bdf'
        f06_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.f06'
        loads_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.blk'
        nid_csv_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.nid'
        eid_csv_filename = MODEL_PATH / 'bwb' / 'bwb_saero_trim.eid'
        log = SimpleLogger(level='warning')
        model = read_bdf(bdf_filename, log=log)
        # export_caero_mesh(
        #     model,
        #     caero_bdf_filename=aerobox_caero_filename,
        #     is_aerobox_model=True,
        #     pid_method='caero',
        #     rotate_panel_angle_deg=45.,
        #     write_panel_xyz=False)
        export_caero_mesh(
            model,
            caero_bdf_filename=aerobox_caero_filename,
            is_aerobox_model=True,
            pid_method='caero',
            write_panel_xyz=True)

        # auto-generated files
        trim_results = f06_to_pressure_loads(
            f06_filename, aerobox_caero_filename,
            loads_filename,
            log=log, nlines_max=1_000_000)

        # set the file names?
        trim_results = f06_to_pressure_loads(
            f06_filename, aerobox_caero_filename,
            loads_filename,
            nid_csv_filename=nid_csv_filename,
            eid_csv_filename=eid_csv_filename,
            log=log, nlines_max=1_000_000)

        argv = [
            'f06', 'plot_144', str(f06_filename),
            '--aerobox', str(aerobox_caero_filename)]
        #log = get_logger(log=None, level=None, encoding='utf-8')
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB,
                     show=False, log=log)

        argv = ['f06', 'plot_144', str(f06_filename)]
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB,
                     show=False, log=log)

    def test_f06_trim_freedlm(self):
        """tests read_f06_trim"""
        bdf_filename = AERO_PATH / 'freedlm' / 'freedlm.bdf'
        caero_filename = AERO_PATH / 'freedlm' / 'freedlm_caero.bdf'
        model = read_bdf(bdf_filename, debug=None)
        export_caero_mesh(
            model,
            caero_bdf_filename=caero_filename,
            is_aerobox_model=True,
            pid_method='caero',
            write_panel_xyz=True)
        model2 = read_bdf(bdf_filename, debug=None)
        #print(f'nnodes = {len(model2.nodes)}')
        #print(f'nelements = {len(model2.elements)}')

        f06_filename = AERO_PATH / 'freedlm' / 'freedlm.f06'
        trim_results = read_f06_trim(
            f06_filename,
            log=None, nlines_max=1_000_000, debug=None)['trim_results']
        assert len(trim_results.aero_force.keys()) == 2
        assert len(trim_results.aero_pressure.keys()) == 2
        assert len(trim_results.controller_state.keys()) == 2
        assert len(trim_results.trim_variables.keys()) == 2
        assert len(trim_results.structural_monitor_loads.keys()) == 2
        keys = list(trim_results.aero_pressure)
        key0 = keys[0]
        #print(f'keys = {list(trim_results.aero_pressure.keys())}')
        apress = trim_results.aero_pressure[key0]
        print(trim_results.aero_force[key0])
        nids = apress.nodes
        press = apress.pressure

        aforce = trim_results.aero_force[key0]
        nids = aforce.nodes
        force = aforce.force
        #print(eids)
        #print(press)
        #print(f'npressure = {len(press)}')
        #print(f'nforce = {len(force)}')

    def test_f06_trim_aerobeam(self):
        """tests read_f06_trim"""
        f06_filename = AERO_PATH / 'aerobeam.f06'
        trim_results = read_f06_trim(f06_filename,
                                     log=None, nlines_max=1_000_000, debug=None)['trim_results']
        assert len(trim_results.aero_force.keys()) == 0
        assert len(trim_results.aero_pressure.keys()) == 0
        assert len(trim_results.controller_state.keys()) == 2
        assert len(trim_results.trim_variables.keys()) == 2
        assert len(trim_results.structural_monitor_loads.keys()) == 2

    #def test_f06_trim_cpmopt(self):
        #"""tests read_f06_trim"""
        #f06_filename = MODEL_PATH / 'aero' / 'cpmopt.f06'
        #trim_results = read_f06_trim(f06_filename,
        #                             log=None, nlines_max=1_000_000, debug=None)['trim_results']
        #assert len(trim_results.aero_force.keys()) == 0
        #assert len(trim_results.aero_pressure.keys()) == 0
        #assert len(trim_results.controller_state.keys()) == 4
        #assert len(trim_results.trim_variables.keys()) == 4
        #assert len(trim_results.structural_monitor_loads.keys()) == 4

    def test_split_float_colon(self):
        """tests split_float_colon"""
        a = split_float_colons('1:')
        b = split_float_colons('1:5')
        c = split_float_colons(':4')
        d = split_float_colons(None)

        assert a == [1.0, None], a
        assert b == [1.0, 5.0], b
        assert c == [None, 4.0], c
        assert d is None, d
        with self.assertRaises(AssertionError):
            split_float_colons('1:5:2')

    def test_split_int_colon(self):
        """tests split_int_colon"""
        a = split_int_colon('1:5')
        assert a == [1, 2, 3, 4, 5], a

        b = split_int_colon('1:5:2')
        assert b == [1, 3, 5], b

        ## TODO: wrong?
        b2 = split_int_colon('1:6:2')
        assert b2 == [1, 3, 5], b2

        c = split_int_colon(':4')
        assert c == [0, 1, 2, 3, 4], c

        d = split_int_colon('1:5,10:15')
        assert d == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], d

        e = split_int_colon('10:15,1:5')
        assert e == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], e

        f = split_int_colon('1,3,5')
        assert f == [1, 3, 5], f

    def test_parse_real_row_lines(self):
        lines = [
            '        1)    1.1010E+01  1.3762E+00 -4.2021E+00 -5.2526E-01',
        ]
        row1, data1 = _parse_real_row_lines(lines)
        assert np.array_equal(row1, [1, 2, 3, 4])
        assert np.array_equal(data1, [1.1010E+01, 1.3762E+00, -4.2021E+00, -5.2526E-01])

    def test_parse_complex_row_lines(self):
        lines = [
            '        1) -3.5846E+01,-1.3275E+02  -1.5510E+01, 2.3578E-01  -3.2339E+01,-4.9373E+00   6.8078E+01, 1.3428E+01   3.0262E+01, 2.4554E+01',
            '        6)  1.5360E-04,-1.1042E-04  -4.7606E-04, 2.3069E-04   1.0359E-03,-1.5668E-04  -1.3075E-03, 7.8472E-04   2.3471E-04,-4.8359E-04',
        ]
        row1, data1 = _parse_complex_row_lines(lines)
        assert np.array_equal(row1, [1, 2, 3, 4, 5,
                                        6, 7, 8, 9, 10]), row1
        reals = [
            -3.5846E+01, -1.5510E+01, -3.2339E+01,  6.8078E+01, 3.0262E+01,
            1.5360E-04,  -4.7606E-04,  1.0359E-03, -1.3075E-03, 2.3471E-04]
        imags = [
            -1.3275E+02, 2.3578E-01, -4.9373E+00, 1.3428E+01,  2.4554E+01,
            -1.1042E-04, 2.3069E-04, -1.5668E-04, 7.8472E-04, -4.8359E-04]
        assert np.array_equal(data1.real, reals)
        assert np.array_equal(data1.imag, imags)

        lines = [
            '        1) -3.5846E+01,-1.3275E+02',
            '        6)  1.5360E-04,-1.1042E-04',
        ]
        row2, data2 = _parse_complex_row_lines(lines)
        assert np.array_equal(row2, [1, 6])
        assert np.array_equal(data2.real, [-3.5846E+01, 1.5360E-04])
        assert np.array_equal(data2.imag, [-1.3275E+02, -1.1042E-04])


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
