# encoding: utf-8
"""
SOL 145 plotter

kfreq = ωc/(2V)
"""
from typing import Optional, cast
import numpy as np
#import PySide
try:
    import matplotlib.pyplot as plt  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

# if you're on linux and you don't have a backend, add this...
# we'd add it here, but it breaks PySide/PySide2's QApplication...
#try:  # pragma: no cover
    #plt.figure()
    #plt.close()
#except Exception:  # pragma: no cover
    #plt.switch_backend('Agg')


from cpylog import get_logger2, SimpleLogger
from pyNastran.f06.flutter_response import FlutterResponse, get_flutter_units
from pyNastran.utils import PathLike
from pyNastran.utils.numpy_utils import float_types


def make_flutter_response(f06_filename: PathLike,
                          f06_units=None, out_units=None,
                          log: Optional[SimpleLogger]=None) -> dict[int, FlutterResponse]:
    """
    Creates the FlutterResponse object

    Parameters
    ----------
    f06_filename : str
        the filename to plot
    f06_units : dict[name]=unit; default=None
        f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
    out_units : dict[name]=unit; default=None
        out_units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf'}
    Returns
    -------
    flutters : dict
        key : int
           subcase_id
        value : FlutterResponse()

    """
    f06_units = get_flutter_units(f06_units)
    out_units = get_flutter_units(out_units)

    if log is None:
        log = get_logger2(log=None, debug=True, encoding='utf-8')
    flutters = {}
    iline = 0

    # 1 is the default subcase number
    subcase = 1
    results = []
    modes = []

    configuration = None
    xysym = None
    xzsym = None
    mach = None
    density_ratio = None
    method = None
    found_flutter_summary = False

    log.info('f06_filename = %r' % f06_filename)
    with open(f06_filename, 'r') as f06_file:
        while 1:
            nblank = 0
            line = f06_file.readline()
            iline += 1
            #log.debug('line%ia = %r' % (iline, line))
            while 'SUBCASE ' not in line and 'FLUTTER  SUMMARY' not in line:
                line = f06_file.readline()
                iline += 1
                if not line:
                    nblank += 1
                if nblank == 100:
                    #print(line.strip())
                    break
            if nblank == 100:
                break
            #if 'FLUTTER  SUMMARY' in line:
                #found_flutter_summary = True

            #log.debug('line%ib = %r' % (iline, line))
            if 'SUBCASE' in line[109:]:
                sline = line.strip().split()
                isubcase = sline.index('SUBCASE')
                new_subcase = int(sline[isubcase + 1])
                #print('subcasei=%r iline=%s' % (new_subcase, iline))
                if new_subcase > subcase:
                    log.debug('subcase=%s -> new_subcase=%s' % (subcase, new_subcase))
                    log.debug('modes1 = %s' % modes)
                    flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                              mach, density_ratio, method,
                                              modes, results,
                                              in_units=f06_units)
                    flutter.set_out_units(out_units)
                    #_remove_neutrinos(flutter, log)
                    flutters[subcase] = flutter
                    modes = []
                    results = []

                    subcase = new_subcase
                    #break
                continue

            #log.debug('line%i_FSa = %r' % (iline, line))
            last_line = None
            while 'FLUTTER  SUMMARY' not in line:
                last_line = line
                line = f06_file.readline()
                #log.debug('i=%s %s' % (iline, line.strip().replace('   ', ' ')))

                if '* * * END OF JOB * * *' in line:
                    last_line = None
                    break
                iline += 1
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    log.warning('breaking on nblank=100')
                    break

            if '* * * END OF JOB * * *' in line:
                break
            if nblank == 100:
                log.warning('breaking on nblank=100')
                break
            if 'FLUTTER  SUMMARY' in line:
                found_flutter_summary = True

            # pulls the subcase id for the first subcase
            if last_line is not None:
                #log.debug('line%i_FSb = %r' % (iline, line))
                #log.debug('line%i_FSb = %r' % (iline-1, last_line.replace('     ', ' ')))
                sline = last_line.strip().split()
                try:
                    isubcase = sline.index('SUBCASE')
                    subcase = int(sline[isubcase + 1])
                except ValueError:
                    log.error(f"expected: SUBCASE line -> ['SUBCASE', 1]")
                    log.error(f'found: i={iline} sline={sline}')
                    log.error(f'assuming subcase=1')
                    subcase = 1
                    #raise
                log.debug('subcase = %s' % subcase)

            configuration_sline = f06_file.readline().split()
            #log.error(f'configuration_sline={configuration_sline}')
            iline += 1
            configuration = configuration_sline[2]
            xysym = configuration_sline[5]
            xzsym = configuration_sline[8]
            #print(configuration, xysym, xzsym)

            # ['POINT', '=', '30', 'METHOD', '=', 'PKNL']
            point_sline = f06_file.readline().split()
            iline += 1
            mode = int(point_sline[2])
            method = point_sline[-1]  # 13 for PN, 5 for PK

            #log.debug(point_sline)
            if method == 'PK':
                mach = float(point_sline[6])
                density_ratio = float(point_sline[10])
                #method = point_sline[13]
                #if mode == 1:
                    #print('# iline mode mach density_ratio method')
                #print(iline, mode, mach, density_ratio, method)

            elif method == 'PKNL':
                mach = None
                density_ratio = None
                #if mode == 1:
                    #print('# iline mode method')
                #print(iline, mode, method)
                f06_file.readline()
                iline += 1
            elif method == 'KE':
                #KFREQ       1./KFREQ       VELOCITY          DAMPING       FREQUENCY          COMPLEX   EIGENVALUE
                    #0.5500  1.8181818E+00   2.2589194E+01  -2.4541089E-03   3.3374373E+00  -5.6140173E-01     4.5752050E+02
                mach = None
                density_ratio = None
                #if mode == 1:
                    #print('# iline mode method')
                #print(iline, mode, method)
                f06_file.readline()
                iline += 1
            else:  # pragma: no cover
                raise NotImplementedError(f'method={method!r} point_sline={point_sline}')

            found_existing_mode = mode in modes
            #if found_existing_mode:
                #log.warning('found existing mode %s...' % mode)
                #print('nresults = %s' % len(results))
                #continue
            #else:
                #modes.append(mode)

            # blanks
            f06_file.readline()
            f06_file.readline()
            iline += 2

            lines = []

            # KFREQ  1./KFREQ                      VELOCITY  DAMPING  FREQUENCY   COMPLEX EIGENVALUE - PK
            # KFREQ  1./KFREQ  DENSITY   MACH NO.  VELOCITY  DAMPING  FREQUENCY   COMPLEX EIGENVALUE - PKNL
            if method == 'PK':
                nvalues = 7
            elif method == 'PKNL':
                nvalues = 9
            elif method == 'KE':
                #KFREQ       1./KFREQ       VELOCITY          DAMPING       FREQUENCY          COMPLEX   EIGENVALUE
                    #0.5500  1.8181818E+00   2.2589194E+01  -2.4541089E-03   3.3374373E+00  -5.6140173E-01     4.5752050E+02
                nvalues = 7
            else:  # pragma: no cover
                raise NotImplementedError(method)

            sline = [None] * nvalues
            while len(sline) == nvalues:
                sline = f06_file.readline().split()
                iline += 1

                is_line = (
                    sline and
                    'PAGE' not in sline and
                    'INFORMATION' not in sline and
                    'EIGENVALUE' not in sline and
                    'USER' not in sline
                )
                if is_line:
                    #print('sline = %s' % sline)
                    lines.append(sline)

            if found_existing_mode:
                results[mode-1].extend(lines)
            else:
                results.append(lines)
                modes.append(mode)
            #print('')

        log.debug('modes = %s' % modes)
        #if not found_flutter_summary:
            #print(line)
            #raise RuntimeError("failed to find 'FLUTTER SUMMARY'")
        flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                  mach, density_ratio, method,
                                  modes, results,
                                  in_units=f06_units)
        flutter.set_out_units(out_units)
        flutters[subcase] = flutter
    return flutters

def plot_flutter_f06(f06_filename: PathLike,
                     f06_units: Optional[dict[str, str]]=None,
                     out_units: Optional[dict[str, str]]=None,
                     plot_type: str='tas',
                     modes: Optional[list[int]]=None,
                     plot_vg: bool=False,
                     plot_vg_vf: bool=False,
                     plot_root_locus: bool=False,
                     plot_kfreq_damping: bool=False,
                     xlim: Optional[list[float]]=None,
                     ylim_damping: Optional[list[float]]=None,
                     ylim_freq: Optional[list[float]]=None,
                     ylim_kfreq: Optional[list[float]]=None,
                     vd_limit: Optional[float]=None,
                     damping_limit: Optional[float]=None,
                     nopoints: bool=False,
                     noline: bool=False,
                     export_csv_filename: Optional[str]=None,
                     export_zona_filename: Optional[str]=None,
                     export_veas_filename: Optional[str]=None,
                     export_f06_filename: Optional[str]=None,
                     vg_filename: Optional[str]=None,
                     vg_vf_filename: Optional[str]=None,
                     root_locus_filename: Optional[str]=None,
                     kfreq_damping_filename: Optional[str]=None,
                     subcases: Optional[list[int]]=None,
                     ncol: int=0,
                     plot: bool=True, show: bool=True, clear: bool=False, close: bool=False,
                     log: Optional[SimpleLogger]=None) -> dict[int, FlutterResponse]:
    """
    Plots a flutter (SOL 145) deck

    Parameters
    ----------
    f06_filename : str
        the filename to plot
    f06_units : dict[name]=unit; default=None
        f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
    out_units : dict[name]=unit; default=None
        out_units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf'}
    modes : bool; default=None
        specifies the modes to plot
    plot_type : str; default='tas'
        'tas' : true airspeed
        'eas' : equivalent airspeed
        'alt' : altitude
        'dynamic_pressure' : dynamic pressure
        'mach' : Mach number
    plot_vg : bool; default=False
        make a V-damping plot
    plot_vg_vf : bool; default=False
        make a V-damping/V-freq plot
    plot_root_locus : bool; default=False
        make a Re/Imag plot
    plot_kfreq_damping : bool; default=False
        make a kfreq-damping plot
    show : bool; default=True
        call plt.show()
    xlim : bool; default=None
        xlimits for the V-g and V-g/V-f plots
    ylim_damping : bool; default=None
        ylimits for the V-g plots
    ylim_freq : bool; default=None
        ylimits for the V-f plots
    nopoints : bool; default=False
        suppress the points
    noline : bool; default=False
        suppress the lines
    subcases: list[int]; default=None
        the list of subcases that should be considered
    export_csv_filename : Optional[str]; default=None -> no csv
        the csv filename to dump

    Returns
    -------
    flutters : dict
        key : int
           subcase_id
        value : FlutterResponse()

    Supports:
     o single subcase
     o single subcase, no subcase marker
     o multiple subcases
     o PK
     o PKNL
       o calculation of:
         - equivalent airspeed
         - dynamic pressure
         - altitude

    Doesn't support:
     o SOL 200
     o fixing mode switching problem
     o fixing unconverged points

    """
    assert vd_limit is None or isinstance(vd_limit, float_types), vd_limit
    assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit
    flutters = make_flutter_response(
        f06_filename, f06_units=f06_units, out_units=out_units, log=log)

    if plot:
        make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                           plot_type,
                           plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                           nopoints, noline, ncol=ncol,
                           vd_limit=vd_limit, damping_limit=damping_limit,
                           export_csv_filename=export_csv_filename,
                           export_zona_filename=export_zona_filename,
                           export_veas_filename=export_veas_filename,
                           export_f06_filename=export_f06_filename,
                           vg_filename=vg_filename,
                           vg_vf_filename=vg_vf_filename,
                           root_locus_filename=root_locus_filename,
                           kfreq_damping_filename=kfreq_damping_filename,
                           subcases=subcases,
                           show=show, clear=clear, close=close)
    return flutters

def make_flutter_plots(modes: list[int], flutters: dict[int, FlutterResponse],
                       xlim: Optional[list[float]],
                       ylim_damping: Optional[list[float]],
                       ylim_freq: Optional[list[float]],
                       ylim_kfreq: Optional[list[float]],
                       plot_type: str,
                       plot_vg: bool,
                       plot_vg_vf: bool,
                       plot_root_locus: bool,
                       plot_kfreq_damping: bool,
                       nopoints: bool,
                       noline: bool,
                       ncol: int=0,
                       legend: bool=True,
                       vd_limit: Optional[float]=None,
                       damping_limit: Optional[float]=None,
                       export_csv_filename: Optional[str]=None,
                       export_zona_filename: Optional[str]=None,
                       export_veas_filename: Optional[str]=None,
                       export_f06_filename: Optional[str]=None,
                       vg_filename: Optional[str]=None,
                       vg_vf_filename: Optional[str]=None,
                       root_locus_filename: Optional[str]=None,
                       kfreq_damping_filename: Optional[str]=None,
                       subcases: Optional[list[int]]=None,
                       show: bool=True, clear: bool=False, close: bool=False,
                       log: Optional[SimpleLogger]=None) -> None:
    """actually makes the flutter plots"""
    assert vd_limit is None or isinstance(vd_limit, float_types), vd_limit
    assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit

    assert len(flutters) > 0, flutters
    subcases_flutter_set = set(list(flutters.keys()))
    if subcases is None:
        subcases_set = subcases_flutter_set
    else:
        subcases_set = set(subcases)
    missing_cases_set = subcases_set - subcases_flutter_set

    log = get_logger2(log=log, debug=True, encoding='utf-8')
    if missing_cases_set:
        missing_cases_list = list(missing_cases_set)
        missing_cases_list.sort()
        log.warning(f'missing subcases={missing_cases_list}')

    for subcase, flutter in sorted(flutters.items()):
        if subcase not in subcases_set:
            continue
        flutter = cast(FlutterResponse, flutter)
        _make_flutter_subcase_plot(
            modes, flutter, subcase, xlim, ylim_damping, ylim_freq, ylim_kfreq,
            plot_type, plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
            nopoints, noline,
            ncol=ncol, legend=legend,
            vd_limit=vd_limit, damping_limit=damping_limit,
            vg_filename=vg_filename,
            vg_vf_filename=vg_vf_filename,
            root_locus_filename=root_locus_filename,
            kfreq_damping_filename=kfreq_damping_filename,
            show=show, clear=clear, close=close, log=log)

        if export_csv_filename:
            flutter.export_to_csv(export_csv_filename, modes=modes)
        if export_zona_filename:
            flutter.export_to_zona(export_zona_filename, modes=modes,
                                   xlim=xlim, plot_type=plot_type)
        if export_veas_filename:
            flutter.export_to_veas(export_veas_filename, modes=modes)
        if export_f06_filename:
            flutter.export_to_f06(export_f06_filename, modes=modes)

    if show:
        plt.show()
    #if close:
        #plt.close()


def _make_flutter_subcase_plot(modes, flutter: FlutterResponse, subcase: int,
                               xlim, ylim_damping, ylim_freq, ylim_kfreq,
                               plot_type: str,
                               plot_vg: bool,
                               plot_vg_vf: bool,
                               plot_root_locus: bool,
                               plot_kfreq_damping: bool,
                               nopoints: bool,
                               noline: bool,
                               ncol: int=0,
                               legend: bool=True,
                               vd_limit: Optional[float]=None,
                               damping_limit: Optional[float]=None,
                               vg_filename: Optional[str]=None,
                               vg_vf_filename: Optional[str]=None,
                               root_locus_filename: Optional[str]=None,
                               kfreq_damping_filename: Optional[str]=None,
                               show: bool=True, clear: bool=False, close: bool=False,
                               log: Optional[SimpleLogger]=None):
        #_remove_neutrinos(flutter, log)
        flutter.nopoints = nopoints
        flutter.noline = noline
        if plot_vg:
            filenamei = None if vg_filename is None else (vg_filename % subcase)
            flutter.plot_vg(modes=modes,
                            plot_type=plot_type,
                            xlim=xlim, ylim_damping=ylim_damping,
                            ncol=ncol,
                            #vd_limit=vd_limit,
                            png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_vg_vf:
            filenamei = None if vg_vf_filename is None else (vg_vf_filename % subcase)
            flutter.plot_vg_vf(modes=modes,
                               plot_type=plot_type,
                               xlim=xlim,
                               ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                               vd_limit=vd_limit, damping_limit=damping_limit,
                               ncol=ncol,
                               legend=legend,
                               png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_root_locus:
            filenamei = None if root_locus_filename is None else (root_locus_filename % subcase)
            flutter.plot_root_locus(modes=modes,
                                    fig=None, axes=None,
                                    eigr_lim=None, eigi_lim=None,
                                    ncol=ncol,
                                    clear=clear, legend=True,
                                    png_filename=filenamei,
                                    show=False, close=close)

        if plot_kfreq_damping:
            filenamei = None if kfreq_damping_filename is None else (kfreq_damping_filename % subcase)
            flutter.plot_kfreq_damping(modes=modes,
                                       plot_type=plot_type,
                                       ylim_damping=ylim_damping,
                                       ylim_kfreq=ylim_kfreq,
                                       vd_limit=vd_limit, damping_limit=damping_limit,
                                       ncol=ncol,
                                       png_filename=filenamei, show=False, clear=clear, close=close)


def _remove_neutrinos(flutter: FlutterResponse, log: SimpleLogger):  # pragma: no cover
    str(flutter)
    isave, ifilter = (
        _find_modes_to_keep(flutter, log, tol=1e-8))

    # fix mode switching
    results = flutter.results
    real = results[isave, :, flutter.ieigr]
    imag = results[isave, :, flutter.ieigi]

    # find modes that cross (5 modes, 22 points)
    #imag.shape = (5, 22)
    velocity = results[:, :, flutter.ivelocity]
    damping = results[:, :, flutter.idamping]
    from scipy.interpolate import CubicSpline
    for imode in isave:
        V = velocity[imode]
        g = damping[imode]
        #isort = np.argsort(V)
        V, isort = np.unique(V, return_index=True)
        V = V[isort]
        g = g[isort]
        Vmin = V.min()
        Vmax = V.max()
        #CubicSpline(x,y,bc_type='natural')
        func = CubicSpline(V, g, axis=0, bc_type='not-a-knot', extrapolate=None)
        num = 200
        vx = np.linspace(Vmin, Vmax, num=num)
        gx = func(vx)
        y = 1
        #if gx.max() < -0.1:
            #results[imode, :, :] = np.nan

    radius = np.sqrt(real ** 2 + imag ** 2)

    # R sin(t) = i
    # R cos(t) = r
    # tan(t) = i/r
    #
    # i
    # ^     * (R, t) / (real, imag)
    # |   /
    # | /
    # +------> r
    theta = np.arctan2(imag, real)  # x / y
    theta_deg = np.degrees(theta)
    flutter.results[ifilter, :, :] = np.nan
    #flutter.results = flutter.results[isave, :, :]
    return


def _find_modes_to_keep(flutter: FlutterResponse,
                        log: SimpleLogger,
                        tol: float=1e-8) -> tuple[np.ndarray, np.ndarray]:  # pragma: no cover
    """
    FlutterResponse:
        subcase= 1
        xysym  = 'ASYMMETRIC'
        xzsym  = 'SYMMETRIC'
        f06_units  = {'velocity': 'in/s', 'density': 'slinch/in^3', 'altitude': 'ft', 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        out_units  = {'velocity': 'in/s', 'density': 'slinch/in^3', 'altitude': 'ft', 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        names  = ['kfreq', '1/kfreq', 'velocity', 'damping', 'freq', 'eigr', 'eigi']
        method  = PK
        modes  = [ 1  2  3  4  5  6  7  8  9 10]; n=10
        results.shape = (10, 22, 7); (nmodes, npoint, nresults)
    """
    #flutter.results.shape = (10, 22, 7)

    # find the delta for each mode for each result
    results = flutter.results
    mini = results.min(axis=1)
    maxi = results.max(axis=1)
    delta = maxi - mini

    #floatmodestr='fixed',
    np.set_printoptions(precision=4, suppress=True)

    # lets reduce this downn to to 10 delta reals and delta imaginary eigenvalues
    #delta.shape = (10, 7)
    abs_dreal = np.abs(delta[:, flutter.ieigr])
    abs_dimag = np.abs(delta[:, flutter.ieigi])

    ddamp = delta[:, flutter.idamping]
    damping = results[:, :, flutter.idamping]
    abs_damp = np.abs(damping).max(axis=1)
    abs_ddamp = np.abs(ddamp)

    dfreq = delta[:, flutter.ifreq]
    freq = results[:, :, flutter.ifreq]
    abs_freq = np.abs(freq).max(axis=1)
    abs_dfreq = np.abs(dfreq)

    # filter the neutrino cases (dreal=0 and dimag=0)
    # find the cases where the delta
    damp_ztol = 1e-4
    damp_atol = 0.2

    freq_ztol = 1e-4
    freq_atol = -0.1
    ifilter = np.where(
        ((abs_dreal < tol) & (abs_dimag < tol)) |  # remove root-locus dots
        ((abs_ddamp < damp_ztol) & (abs_damp > damp_atol)) | # remove flat damping lines far from the damping axis
        ((abs_dfreq < freq_ztol) & (abs_freq > freq_atol))   # remove flat freq lines far from the damping axis
    )
    nmodes = results.shape[0]
    iall = np.arange(nmodes)
    isave = np.setdiff1d(iall, ifilter)
    #isave[ineutrino] = 0
    #isave = np.where(
        #((abs_dreal > tol) | (abs_dimag > tol)) |
    #)


    #return isave
    #iimag = np.where(dimag != 0.)[0]
    #log
    return isave, ifilter


if __name__ == '__main__':  # pragma: no cover
    plot_flutter_f06('bah_plane.f06')
