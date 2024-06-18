"""
tests:
 - plot sol_145
"""
import os
import unittest

import numpy as np
from scipy.spatial import KDTree
from scipy.interpolate import interp1d, splrep, splev

from cpylog import SimpleLogger, get_logger2
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
from pyNastran.f06.utils import (
    split_float_colons, split_int_colon,
    cmd_line_plot_flutter, cmd_line as cmd_line_f06)
from pyNastran.f06.parse_flutter import (
    plot_flutter_f06, make_flutter_plots, make_flutter_response,
    FlutterResponse,
)
from pyNastran.f06.parse_trim import read_f06_trim
from pyNastran.f06.dev.read_sol_200 import plot_sol_200  # read_sol_200

DIRNAME = os.path.dirname(__file__)
PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


class TestFlutter(unittest.TestCase):

    def test_f06_pt145(self):
        """tests read_f06_trim"""
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        f06_filename = os.path.join(MODEL_PATH, 'aero', 'pt145.f06')
        #trim_results = read_f06_trim(f06_filename,
        #                             log=None, nlines_max=1_000_000, debug=None)
        #assert len(trim_results.aero_force.keys()) == 0
        #assert len(trim_results.aero_pressure.keys()) == 0
        #assert len(trim_results.controller_state.keys()) == 0
        #assert len(trim_results.trim_variables.keys()) == 0
        #assert len(trim_results.structural_monitor_loads.keys()) == 4
        argv = ['f06', 'plot_145', f06_filename, '--tas',
                '--out_units', 'english_in']
        cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB,
                              show=False, log=log)

    def test_plot_flutter_bah(self):
        """tests plot_flutter_f06"""
        f06_filename = os.path.join(MODEL_PATH, 'aero', 'bah_plane', 'bah_plane.f06')
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        flutters = plot_flutter_f06(
            f06_filename, show=False, close=True,
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True, plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB,
            log=log)
        assert len(flutters) == 2, list(flutters.keys())
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
        f06_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.f06')
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        plot_flutter_f06(
            f06_filename, make_alt=True,
            modes=[2],
            plot_type='alt',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, log=log,
            close=True,
        )

        flutters = plot_flutter_f06(
            f06_filename,
            f06_units=None, out_units=None,
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            export_csv_filename='nastran.csv',
            export_f06_filename='nastran.f06',
            export_veas_filename='nastran.veas',
            export_zona_filename='zona.f06',
            vg_filename='vg_subcase_%i.png',
            vg_vf_filename='vg_vf_subcase_%i.png',
            kfreq_damping_filename='kfreq_damping_subcase_%i.png',
            root_locus_filename='root_locus_subcase_%i.png',
            plot=IS_MATPLOTLIB, show=False, log=log,
            close=True,
        )
        flutter = flutters[1]
        fix_modes_2024(flutter)

        if IS_MATPLOTLIB:
            os.remove('vg_subcase_1.png')
            os.remove('vg_vf_subcase_1.png')
            os.remove('kfreq_damping_subcase_1.png')
            os.remove('root_locus_subcase_1.png')

            os.remove('nastran.csv')
            os.remove('nastran.f06')
            os.remove('nastran.veas')
            os.remove('zona.f06')

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

        plot_flutter_f06(
            f06_filename,
            plot_type='rho',
            f06_units='si', out_units='english_ft',
            plot_vg=True, plot_vg_vf=True, plot_root_locus=True,
            plot_kfreq_damping=True,
            plot=IS_MATPLOTLIB, show=False, close=True, log=log)

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
                               # export_zona_filename=export_zona_filename,
                               # export_veas_filename=export_veas_filename,
                               # export_f06_filename=export_f06_filename,
                               # vg_filename=vg_filename,
                               # vg_vf_filename=vg_vf_filename,
                               # root_locus_filename=root_locus_filename,
                               # kfreq_damping_filename=kfreq_damping_filename,
                               show=False, clear=True, close=True)

            export_csv_filename = os.path.join(DIRNAME, 'csv_%d.csv')
            export_zona_filename = os.path.join(DIRNAME, 'zona_%d.f06')
            export_veas_filename = os.path.join(DIRNAME, 'flutter_%d.veas')
            export_f06_filename = os.path.join(DIRNAME, 'flutter_%d.f06')
            make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                               plot_type,
                               plot_vg=True, plot_vg_vf=True,
                               plot_root_locus=True, plot_kfreq_damping=True,
                               nopoints=True, noline=False,
                               export_csv_filename=export_csv_filename,
                               export_zona_filename=export_zona_filename,
                               export_veas_filename=export_veas_filename,
                               export_f06_filename=export_f06_filename,
                               # vg_filename=vg_filename,
                               # vg_vf_filename=vg_vf_filename,
                               # root_locus_filename=root_locus_filename,
                               # kfreq_damping_filename=kfreq_damping_filename,
                               show=False, clear=True, close=True)

    def test_cmd_line_plot_flutter_0012(self):
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        f06_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.f06')
        argv = ['f06', 'plot_145', f06_filename, '--eas',
                '--in_units', 'si', '--out_units', 'english_in',
                '--modes', '1:', '--ylimdamp', '-.3:', '--export_csv', '--ncol', '2']
        cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)

    def test_cmd_line_plot_flutter_no_input_0012(self):
        """no input???"""
        log = get_logger2(log=None, debug=None, encoding='utf-8')
        f06_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.f06')
        argv = ['f06', 'plot_145', f06_filename, '--eas',
                '--out_units', 'english_in']
        flutters = cmd_line_plot_flutter(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)
        cmd_line_f06(argv=argv, plot=IS_MATPLOTLIB, show=False, log=log)

    def test_fix_modes_0012(self):
        log = SimpleLogger(level='warning')
        f06_filename = os.path.join(MODEL_PATH, 'aero', '2_mode_flutter', '0012_flutter.f06')
        flutter = make_flutter_response(f06_filename, log=log)[1]
        fix_modes_2024(flutter)

    def test_fix_modes_constant(self):
        """constant to check dxyz"""
        # nmodes, nvel
        eigr = np.array([
            [0.,1., 2.],
            [0.,1., 2.],
        ])
        eigi = np.array([
            [1.,1.,1.],
            [10.,10,10],
        ])
        nmodes, nvel = eigr.shape
        out = _fix_modes(eigr, eigi, nmodes, nvel, kmodes=0)

    def test_fix_modes_linear(self):
        """constant to check dxyz"""
        # nmodes, nvel
        eigr = np.array([
            [0.,1., 2.],
            [0.,1., 2.],
        ])
        eigi = np.array([
            [1.,1.,1.],
            [0.,.5,1.1],
        ])
        nmodes, nvel = eigr.shape
        out = _fix_modes(eigr, eigi, nmodes, nvel, kmodes=0)


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
        mode_switching
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
    data0 = np.column_stack([eigr[:,0],eigi[:, 0]])
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

        imode_next = ixyz[:,0]
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

        out[:,ivel] = isort
        imode_next
        data0 = all_data[:, ivel, :]
        tree0 = KDTree(data0)
    return out


class TestF06Utils(unittest.TestCase):
    def test_opt_aerobeam(self):
        """tests optimization"""
        f06_filename = os.path.join(MODEL_PATH, 'aero', 'aerobeam.f06')
        png_filename = os.path.join(MODEL_PATH, 'aero', 'aerobeam.png')
        plot_sol_200(f06_filename, png_filename=png_filename,
                     show=True)
        #read_sol_200(f06_filename)

    def test_opt_mdb200(self):
        """tests optimization"""
        f06_filename = os.path.join(MODEL_PATH, 'other', 'mdb200.f06')
        png_filename = os.path.join(MODEL_PATH, 'other', 'mdb200.png')
        plot_sol_200(f06_filename, png_filename=png_filename,
                     show=True)

    def test_f06_trim_freedlm(self):
        """tests read_f06_trim"""
        f06_filename = os.path.join(MODEL_PATH, 'aero', 'freedlm', 'freedlm.f06')
        trim_results = read_f06_trim(f06_filename,
                                     log=None, nlines_max=1_000_000, debug=None)
        assert len(trim_results.aero_force.keys()) == 2
        assert len(trim_results.aero_pressure.keys()) == 2
        assert len(trim_results.controller_state.keys()) == 2
        assert len(trim_results.trim_variables.keys()) == 2
        assert len(trim_results.structural_monitor_loads.keys()) == 2

    def test_f06_trim_aerobeam(self):
        """tests read_f06_trim"""
        f06_filename = os.path.join(MODEL_PATH, 'aero', 'aerobeam.f06')
        trim_results = read_f06_trim(f06_filename,
                                     log=None, nlines_max=1_000_000, debug=None)
        assert len(trim_results.aero_force.keys()) == 0
        assert len(trim_results.aero_pressure.keys()) == 0
        assert len(trim_results.controller_state.keys()) == 2
        assert len(trim_results.trim_variables.keys()) == 2
        assert len(trim_results.structural_monitor_loads.keys()) == 2

    #def test_f06_trim_cpmopt(self):
        #"""tests read_f06_trim"""
        #f06_filename = os.path.join(MODEL_PATH, 'aero', 'cpmopt.f06')
        #trim_results = read_f06_trim(f06_filename,
        #                             log=None, nlines_max=1_000_000, debug=None)
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

        assert a == [1.0, None], a
        assert b == [1.0, 5.0], b
        assert c == [None, 4.0], c
        with self.assertRaises(AssertionError):
            split_float_colons('1:5:2')

    def test_split_int_colon(self):
        """tests split_int_colon"""
        a = split_int_colon('1:5')
        assert a == [1, 2, 3, 4, 5], a

        b = split_int_colon('1:5:2')
        assert b == [1, 3, 5], b

        c = split_int_colon(':4')
        assert c == [0, 1, 2, 3, 4], c

        d = split_int_colon('1:5,10:15')
        assert d == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], d

        d = split_int_colon('10:15,1:5')
        assert d == [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15], d


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
