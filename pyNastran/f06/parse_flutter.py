# encoding: utf-8
"""
SOL 145 plotter

kfreq = ωc/(2V)
"""
from typing import Optional, TextIO, cast, Any
import numpy as np

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
from pyNastran.utils.numpy_utils import float_types, integer_types
from pyNastran.f06.f06_matrix_parser import read_real_eigenvalues
Crossing = tuple[float, float, float]


def make_flutter_response(f06_filename: PathLike,
                          f06_units=None, out_units=None,
                          use_rhoref: bool=False,
                          log: Optional[SimpleLogger]=None) -> tuple[dict[int, FlutterResponse],
                                                               dict[str, Any]]:
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
    use_rhoref: bool; default=False
        False: assume the density in the table is absolute density
        True: assume the density should be defined by sea level density,
              so density is a density ratio


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
    methodi = ''
    found_flutter_summary = False
    ieigenvector = -1
    eigenvectors = []
    eigenvectors_array = None
    eigr_eigi_velocity_list = []
    eigr_eigi_velocity = None
    load_eigenvalues = True
    data = {}
    matrices = {}

    log.info('f06_filename = %r' % f06_filename)
    with open(f06_filename, 'r') as f06_file:
        while 1:
            nblank = 0
            line = f06_file.readline()
            iline += 1
            #log.debug(f'A: line[{iline:d}] = {line!r}')
            if 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R' in line:
                print(f'line = {line}')
                asdf
            elif 'R E A L   E I G E N V A L U E S' in line and load_eigenvalues:
                Mhh, Bhh, Khh = read_real_eigenvalues(f06_file, log, line, iline)
                asdf
                isort = np.argsort(Khh)
                # matrices['MHH'].append(np.diag(Mhh[isort]))
                # matrices['BHH'].append(np.diag(Bhh[isort]))
                # matrices['KHH'].append(np.diag(Khh[isort]))
                matrices['MHH'] = np.diag(Mhh[isort])
                matrices['BHH'] = np.diag(Bhh[isort])
                matrices['KHH'] = np.diag(Khh[isort])
                del Mhh, Bhh, Khh

            while ('SUBCASE ' not in line and
                   'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R' not in line and
                   #'R E A L   E I G E N V A L U E S' not in line and
                   'FLUTTER  SUMMARY' not in line and
                   'EIGENVECTOR FROM THE' not in line):
                line = f06_file.readline()
                iline += 1
                #print(f'i={iline} {line.rstrip()}')
                if not line:
                    nblank += 1
                else:
                    nblank = 0
                if nblank == 100:
                    #print(line.strip())
                    break
            if nblank == 100:
                break
            #if 'FLUTTER  SUMMARY' in line:
                #found_flutter_summary = True

            #log.debug(f'B: line[{iline}] = {line!r}')
            if 'SUBCASE' in line[109:]:
                #log.debug(f'B1: subcase')
                ieigenvector = -1
                eigenvectors = []
                sline = line.strip().split()
                isubcase = sline.index('SUBCASE')
                new_subcase = int(sline[isubcase + 1])
                #print('subcasei=%r iline=%s' % (new_subcase, iline))
                if new_subcase > subcase:
                    log.debug('subcase=%s -> new_subcase=%s' % (subcase, new_subcase))
                    log.debug('modes1 = %s' % modes)
                    response = FlutterResponse(
                        f06_filename, subcase, configuration, xysym, xzsym,
                        mach, density_ratio, method,
                        modes, results,
                        in_units=f06_units,
                        use_rhoref=use_rhoref,
                        eigenvector=eigenvectors_array,
                        eigr_eigi_velocity=eigr_eigi_velocity)
                    response.set_out_units(out_units)
                    #_remove_neutrinos(response, log)
                    flutters[subcase] = response
                    modes = []
                    results = []
                    subcase = new_subcase
                    #break
                continue
            elif 'R E A L   E I G E N V A L U E S' in line:
                Mhh, Bhh, Khh = read_real_eigenvalues(f06_file, log, line, iline)
                isort = np.argsort(Khh)
                # matrices['MHH'].append(np.diag(Mhh[isort]))
                # matrices['BHH'].append(np.diag(Bhh[isort]))
                # matrices['KHH'].append(np.diag(Khh[isort]))
                matrices['MHH'] = np.diag(Mhh[isort])
                matrices['BHH'] = np.diag(Bhh[isort])
                matrices['KHH'] = np.diag(Khh[isort])
            elif 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R' in line:
                iline, line, opgwg = _read_opgwg(f06_file, iline, line)
                data['opgwg'] = opgwg

            #log.debug(f'C: line[{iline:d}]_FSa = {line}')
            last_line = None
            while 'FLUTTER  SUMMARY' not in line:
                #log.debug('i=%s %s' % (iline, line.strip().replace('   ', ' ')))

                iline, line, ieigenvector, methodi = _check_for_eigenvector(
                    f06_file, iline, line, eigr_eigi_velocity_list,
                    eigenvectors, ieigenvector, log)

                #short_line = line.strip().replace('   ', ' ')
                #log.debug(f'i={iline} {short_line!r}')
                last_line = line
                line = f06_file.readline()

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
            #if methodi:  # true for MSC, but not NX
                #assert methodi == method, f'methodi={methodi!r}; method={method!r}'

            if len(eigenvectors):
                eigr_eigi_velocity = np.array(eigr_eigi_velocity_list, dtype='float64') # eigr, eigi, velo
                assert eigr_eigi_velocity.ndim == 2, eigr_eigi_velocity
                eigenvectors_array = np.column_stack(eigenvectors)
                #print(f'eigr_eigi_velocity.shape = {eigr_eigi_velocity.shape}')
                eigenvectors = []
                eigr_eigi_velocity_list = []
                #print(f'eigr_eigi_velocity:\n{eigr_eigi_velocity}')

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
                    'USER' not in sline and
                    '^^^' not in sline
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
        response = FlutterResponse(
            f06_filename, subcase, configuration, xysym, xzsym,
            mach, density_ratio, method,
            modes, results,
            in_units=f06_units,
            use_rhoref=use_rhoref,
            eigenvector=eigenvectors_array,
            eigr_eigi_velocity=eigr_eigi_velocity)
        response.set_out_units(out_units)
        flutters[subcase] = response

    if len(matrices):
        data['matrices'] = matrices
    return flutters, data

def _read_opgwg(f06_file: TextIO, iline: int,
                line: str) -> tuple[int, str, dict[str, np.ndarray]]:
    assert 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R' in line, line
    # REFERENCE POINT = 0
    sline = f06_file.readline().split('=')
    iline += 1
    ref_point = int(sline[1])

    # MO
    line = f06_file.readline()
    #print(f'MO = {line.strip()}')
    iline += 1

    mo_list = []
    for i in range(6):
        #  * 1 2 3 4 5 6 *
        line = f06_file.readline().strip()
        iline += 1
        #print(line)
        assert '*' in line[0], line
        sline = line[1:-1].split()
        assert len(sline) == 6, sline
        mo_list.append(sline)

    # S
    line = f06_file.readline()
    iline += 1
    #print(f'S = {line.strip()}')
    s_list = []
    for i in range(3):
        #  * 1 2 36 *
        line = f06_file.readline().strip()
        iline += 1
        assert '*' in line[0], line
        sline = line[1:-1].split()
        assert len(sline) == 3, sline
        s_list.append(sline)

    #---------------------------------------------
    # DIRECTION
    line = f06_file.readline()
    #print(f'dir = {line.strip()}')
    iline += 1

    # MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.
    line = f06_file.readline()
    iline += 1

    mass_cg_list = []
    for i in range(3):
        #  X 1 2 3
        #  Y 1 2 3
        #  Z 1 2 3
        line = f06_file.readline().strip()
        iline += 1
        sline = line[1:].split()
        assert len(sline) == 4, sline
        mass_cg_list.append(sline)

    # I(S)
    line = f06_file.readline()
    iline += 1
    is_list = []
    #print('IS = ', line.strip())
    for i in range(3):
        #  * 1 2 3 *
        #  * 1 2 3 *
        #  * 1 2 3 *
        line = f06_file.readline().strip()
        iline += 1
        assert '*' in line[0], line
        sline = line[1:-1].split()
        assert len(sline) == 3, sline
        #print(f'is[{i}] = {sline}')
        is_list.append(sline)

    # I(Q)
    line = f06_file.readline()
    iline += 1
    iq_list = []
    #print('IQ = ', line.strip())
    for i in range(3):
        #  * 1     *
        #  *   2   *
        #  *     3 *
        line = f06_file.readline().strip()
        iline += 1
        assert '*' in line[0], line
        sline = line[1:-1].split()
        assert len(sline) == 1, sline
        iq_list.append(sline)

    # Q
    line = f06_file.readline()
    iline += 1
    q_list = []
    for i in range(3):
        #  * 1 2 3 *
        #  * 1 2 3 *
        #  * 1 2 3 *
        line = f06_file.readline().strip()
        iline += 1
        assert '*' in line[0], line
        sline = line[1:-1].split()
        assert len(sline) == 3, sline
        q_list.append(sline)

    MO = np.array(mo_list, dtype='float64')
    mass_cg = np.array(mass_cg_list, dtype='float64')
    IS = np.array(is_list, dtype='float64')
    IQ = np.array(iq_list, dtype='float64')
    S = np.array(s_list, dtype='float64')
    Q = np.array(q_list, dtype='float64')

    assert MO.shape == (6, 6), MO.shape
    assert mass_cg.shape == (3, 4), mass_cg.shape
    assert IS.shape == (3, 3), IS.shape
    assert IQ.shape == (3, 1), IQ.shape
    assert S.shape == (3, 3), S.shape
    assert Q.shape == (3, 3), Q.shape
    #print(MO.shape, mass_cg.shape, IS.shape, IQ.shape, S.shape, Q.shape)
    line = f06_file.readline().strip()
    iline += 1

    mass_full = mass_cg[:, 0]
    cg_full = mass_cg[:, 1:]
    mass = mass_full.mean()
    mass_error = mass_full - mass
    xcg = (cg_full[1, 0] + cg_full[2, 0]) / 2.
    ycg = (cg_full[0, 1] + cg_full[2, 1]) / 2.
    zcg = (cg_full[0, 2] + cg_full[1, 2]) / 2.
    xyz_cg = np.array([xcg, ycg, zcg])

    cg_error = cg_full.copy()
    cg_error[1, 0] = cg_error[1, 0] - xcg
    cg_error[2, 0] = cg_error[2, 0] - xcg

    # ycg
    cg_error[0, 1] = cg_error[0, 1] - ycg
    cg_error[2, 1] = cg_error[2, 1] - ycg

    # zcg
    cg_error[0, 2] = cg_error[0, 2] - zcg
    cg_error[1, 2] = cg_error[1, 2] - zcg

    opgwg = {
        'ref_point': ref_point,
        'MO': MO,
        # 'mass_full': mass_full,
        # 'cg_full': cg_full,
        'cg': xyz_cg,
        'mass': mass,
        'mass_error': mass_error.sum(),
        'cg_error': cg_error.sum(),
        'I(S)': IS,
    }
    #print(opgwg)
    return iline, line, opgwg

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
                     freq_tol: float=-1.0,
                     mag_tol: float=-1.0,
                     ivelocity: Optional[int]=None,
                     mode: Optional[int]=None,
                     use_rhoref: bool=False,
                     nopoints: bool=False,
                     noline: bool=False,
                     export_csv_filename: Optional[str]=None,
                     export_zona_filename: Optional[str]=None,
                     export_veas_filename: Optional[str]=None,
                     export_f06_filename: Optional[str]=None,
                     vg_filename: Optional[str]=None,
                     vg_vf_filename: Optional[str]=None,
                     root_locus_filename: Optional[str]=None,
                     modal_participation_filename: Optional[str]=None,
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
    use_rhoref: bool; default=False
        assume rho is density ratio (rho/rhoSL)
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
    assert ivelocity is None or isinstance(ivelocity, integer_types), ivelocity
    #assert mode is None or isinstance(mode, integer_types), mode
    flutters, mass = make_flutter_response(
        f06_filename, f06_units=f06_units, out_units=out_units,
        use_rhoref=use_rhoref, log=log)

    if plot:
        make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                           plot_type,
                           plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                           nopoints, noline, ncol=ncol,
                           ivelocity=ivelocity, mode=mode,
                           vd_limit=vd_limit, damping_limit=damping_limit,
                           freq_tol=freq_tol, mag_tol=mag_tol,
                           export_csv_filename=export_csv_filename,
                           export_zona_filename=export_zona_filename,
                           export_veas_filename=export_veas_filename,
                           export_f06_filename=export_f06_filename,
                           vg_filename=vg_filename,
                           vg_vf_filename=vg_vf_filename,
                           root_locus_filename=root_locus_filename,
                           modal_participation_filename=modal_participation_filename,
                           kfreq_damping_filename=kfreq_damping_filename,
                           subcases=subcases,
                           show=show, clear=clear, close=close)
    return flutters

def make_flutter_plots(modes: list[int],
                       flutters: dict[int, FlutterResponse],
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
                       freq_tol: float=-1.0,
                       mag_tol: float=-1.0,
                       ivelocity: Optional[int]=None,
                       mode: Optional[int]=None,
                       export_csv_filename: Optional[str]=None,
                       export_zona_filename: Optional[str]=None,
                       export_veas_filename: Optional[str]=None,
                       export_f06_filename: Optional[str]=None,
                       vg_filename: Optional[str]=None,
                       vg_vf_filename: Optional[str]=None,
                       root_locus_filename: Optional[str]=None,
                       modal_participation_filename: Optional[str]=None,
                       kfreq_damping_filename: Optional[str]=None,
                       subcases: Optional[list[int]]=None,
                       show: bool=True, clear: bool=False, close: bool=False,
                       log: Optional[SimpleLogger]=None) -> None:
    """actually makes the flutter plots"""
    assert vd_limit is None or isinstance(vd_limit, float_types), vd_limit
    assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit
    assert ivelocity is None or isinstance(ivelocity, integer_types), ivelocity
    #assert mode is None or isinstance(mode, integer_types), mode

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

    for subcase, response in sorted(flutters.items()):
        if subcase not in subcases_set:
            continue
        response = cast(FlutterResponse, response)
        _make_flutter_subcase_plot(
            modes, response, subcase, xlim, ylim_damping, ylim_freq, ylim_kfreq,
            ivelocity, mode,
            plot_type, plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
            nopoints, noline,
            ncol=ncol, legend=legend,
            vd_limit=vd_limit,
            damping_limit=damping_limit,
            freq_tol=freq_tol, mag_tol=mag_tol,
            vg_filename=vg_filename,
            vg_vf_filename=vg_vf_filename,
            root_locus_filename=root_locus_filename,
            modal_participation_filename=modal_participation_filename,
            kfreq_damping_filename=kfreq_damping_filename,
            show=show, clear=clear, close=close, log=log)

        if export_csv_filename:
            response.export_to_csv(export_csv_filename, modes=modes)
        if export_zona_filename:
            response.export_to_zona(export_zona_filename, modes=modes,
                                    xlim=xlim, plot_type=plot_type)
        if export_veas_filename:
            response.export_to_veas(export_veas_filename, modes=modes)
        if export_f06_filename:
            response.export_to_f06(export_f06_filename, modes=modes)

    if show:
        plt.show()
    #if close:
        #plt.close()


def _make_flutter_subcase_plot(modes, response: FlutterResponse,
                               subcase: int,
                               xlim, ylim_damping, ylim_freq, ylim_kfreq,
                               ivelocity: Optional[int],
                               mode: Optional[int],
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
                               freq_tol: float=-1.0,
                               mag_tol: float=-1.0,
                               vg_filename: Optional[str]=None,
                               vg_vf_filename: Optional[str]=None,
                               root_locus_filename: Optional[str]=None,
                               kfreq_damping_filename: Optional[str]=None,
                               modal_participation_filename: Optional[str]=None,
                               show: bool=True, clear: bool=False, close: bool=False,
                               log: Optional[SimpleLogger]=None):
        #_remove_neutrinos(response, log)
        response.nopoints = nopoints
        response.noline = noline
        v_lines = _vd_limit_to_v_lines(vd_limit)
        if plot_vg:
            filenamei = None if vg_filename is None else (vg_filename % subcase)
            response.plot_vg(modes=modes,
                             plot_type=plot_type,
                             xlim=xlim, ylim_damping=ylim_damping,
                             ncol=ncol,
                             v_lines=v_lines,
                             png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_vg_vf:
            filenamei = None if vg_vf_filename is None else (vg_vf_filename % subcase)
            response.plot_vg_vf(modes=modes,
                                plot_type=plot_type,
                                xlim=xlim,
                                ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                                damping_limit=damping_limit,
                                v_lines=v_lines,
                                freq_tol=freq_tol,
                                ivelocity=ivelocity,
                                ncol=ncol, legend=legend,
                                png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_root_locus:
            filenamei = None if root_locus_filename is None else (root_locus_filename % subcase)
            response.plot_root_locus(modes=modes,
                                     fig=None, axes=None,
                                     eigr_lim=None, eigi_lim=None,
                                     freq_tol=freq_tol,
                                     ivelocity=ivelocity,
                                     ncol=ncol,
                                     clear=clear, legend=True,
                                     png_filename=filenamei,
                                     show=False, close=close)
        is_modal_participation = len(response.eigr_eigi_velocity) and ivelocity is not None
        if is_modal_participation:
            assert modal_participation_filename is not None
            filenamei = None if modal_participation_filename is None else (modal_participation_filename % subcase)
            mode_ = [mode] if isinstance(mode, integer_types) or mode is None else mode
            for modei in mode_:
                response.plot_modal_participation(
                    ivelocity, modei,
                    modes=modes,
                    fig=None, axes=None,
                    freq_tol=freq_tol,
                    mag_tol=mag_tol,
                    ncol=ncol,
                    clear=clear, legend=True,
                    png_filename=filenamei,
                    show=False, close=close)

        if plot_kfreq_damping:
            filenamei = None if kfreq_damping_filename is None else (kfreq_damping_filename % subcase)
            response.plot_kfreq_damping(
                modes=modes,
                plot_type=plot_type,
                ylim_damping=ylim_damping,
                ylim_kfreq=ylim_kfreq,
                v_lines=v_lines, damping_limit=damping_limit,
                ncol=ncol,
                png_filename=filenamei,
                show=False, clear=clear, close=close)


def _remove_neutrinos(response: FlutterResponse,
                      log: SimpleLogger) -> None:  # pragma: no cover
    str(response)
    isave, ifilter = (
        _find_modes_to_keep(response, log, tol=1e-8))

    # fix mode switching
    results = response.results
    real = results[isave, :, response.ieigr]
    imag = results[isave, :, response.ieigi]

    # find modes that cross (5 modes, 22 points)
    #imag.shape = (5, 22)
    velocity = results[:, :, response.ivelocity]
    damping = results[:, :, response.idamping]
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
    response.results[ifilter, :, :] = np.nan
    #response.results = response.results[isave, :, :]
    return


def _find_modes_to_keep(response: FlutterResponse,
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
    #response.results.shape = (10, 22, 7)

    # find the delta for each mode for each result
    results = response.results
    mini = results.min(axis=1)
    maxi = results.max(axis=1)
    delta = maxi - mini

    #floatmodestr='fixed',
    np.set_printoptions(precision=4, suppress=True)

    # lets reduce this downn to to 10 delta reals and delta imaginary eigenvalues
    #delta.shape = (10, 7)
    abs_dreal = np.abs(delta[:, response.ieigr])
    abs_dimag = np.abs(delta[:, response.ieigi])

    ddamp = delta[:, response.idamping]
    damping = results[:, :, response.idamping]
    abs_damp = np.abs(damping).max(axis=1)
    abs_ddamp = np.abs(ddamp)

    dfreq = delta[:, response.ifreq]
    freq = results[:, :, response.ifreq]
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


def _check_for_eigenvector(f06_file: TextIO, iline: int, line: str,
                           eigr_eigi_velocity_list: list[Crossing],
                           eigenvectors: list[np.ndarray],
                           ieigenvector: int,
                           log: SimpleLogger) -> tuple[int, str, int]:
    methodi = ''
    if 'EIGENVECTOR FROM THE' in line:
        # '                                               EIGENVECTOR FROM THE  PKNL METHOD
        # '   EIGENVALUE =    -9.88553E-02    1.71977E+01       VELOCITY =     1.52383E+02
        # '
        # '   EIGENVECTOR
        # '                    1.00000E+00    0.00000E+00
        # '                   -3.81052E-03   -6.17235E-04
        ieigenvector += 1
        # 'EIGENVECTOR FROM THE  PKNL METHOD'
        # print(line)
        #log.debug(f'eigenvector iline: {iline}')
        iline, line, eigenvector, eigr, eigi, velocity, methodi = _get_eigenvector(f06_file, iline, line)
        eigr_eigi_velocity_list.append((eigr, eigi, velocity))
        eigenvectors.append(eigenvector)
    return iline, line, ieigenvector, methodi


def _get_eigenvector(f06_file: TextIO, iline: int,
                     line: str) -> tuple[int, str, np.ndarray, float, float, float, str]:
    methodi = line.split('EIGENVECTOR FROM THE')[1].split('METHOD')[0].strip()
    # print(f'imode={ieigenvector+1}; methodi={methodi!r}')
    line = f06_file.readline()
    iline += 1

    if 'VELOCITY' in line:
        # '   EIGENVALUE =    -9.88553E-02    1.71977E+01       VELOCITY =     1.52383E+02
        eig_value, velocity_str = line.split('VELOCITY =')
        eigr_str, eigi_str = eig_value.strip().split('EIGENVALUE =')[1].split()
        velocity = float(velocity_str)
    else:
        eigr_str, eigi_str = line.strip().split('EIGENVALUE =')[1].split()
        velocity = np.nan

    eigr = float(eigr_str)
    eigi = float(eigi_str)
    # print(f'eigr,eigi,velo = {eigr}, {eigi}, {velocity}')

    while 'EIGENVECTOR' not in line:
        line = f06_file.readline()
        iline += 1
    line = f06_file.readline()
    iline += 1

    spline_point_complex_eigenvector_slines = []
    sline = line.split()
    while len(sline) == 2:
        # print('***', iline, line, sline)
        spline_point_complex_eigenvector_slines.append(sline)
        line = f06_file.readline()
        iline += 1
        sline = line.split()
    real_imag = np.array(spline_point_complex_eigenvector_slines, dtype='float64')
    # eigenvector = real_imag[:, 0].flatten() + 1j*real_imag[:, 1].flatten()
    nspline_points = len(spline_point_complex_eigenvector_slines)
    eigenvector = (real_imag[:, 0] + real_imag[:, 1] * 1j).reshape(nspline_points, 1)
    return iline, line, eigenvector, eigr, eigi, velocity, methodi


def _vd_limit_to_v_lines(vd_limit: Optional[float]=None) -> list:
    v_lines = []
    if vd_limit:
        # (name, value, color, linestyle)
        v_lines.append(('VD', vd_limit, 'k', '--'))
        v_lines.append(('1.15*VD', vd_limit, 'k', '-'))
    return v_lines


if __name__ == '__main__':  # pragma: no cover
    plot_flutter_f06('bah_plane.f06')
