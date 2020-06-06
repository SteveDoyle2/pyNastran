# encoding: utf-8
"""
SOL 145 plotter

kfreq = ?c/(2V)
"""
from typing import  Optional, Dict, Union
#import PySide
try:
    import matplotlib.pyplot as plt  # pylint: disable=unused-import
    IS_MATPLOTLIB = True
except ImportError:  # pragma: no cover
    IS_MATPLOTLIB = False

# if you're on linux and you don't have a backend, add this...
# we'd add it here, but it breaks PySide/PySide2's QApplication...
#try:  # pragma: no cover
    #plt.figure()
    #plt.close()
#except:  # pragma: no cover
    #plt.switch_backend('Agg')


from cpylog import get_logger2
from pyNastran.f06.flutter_response import FlutterResponse


def make_flutter_response(f06_filename, f06_units=None, out_units=None, make_alt=False, log=None):
    """
    Creates the FlutterResponse object

    Parameters
    ----------
    f06_filename : str
        the filename to plot
    f06_units : Dict[name]=unit; default=None
        f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
    out_units : Dict[name]=unit; default=None
        out_units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf'}
    Returns
    -------
    flutters : dict
        key : int
           subcase_id
        value : FlutterResponse()

    """
    f06_units = _get_units(f06_units)
    out_units = _get_units(out_units)

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
                                              f06_units=f06_units, out_units=out_units,
                                              make_alt=make_alt)
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
                    log.warning('breaking on nblank=100 a')
                    break

            if '* * * END OF JOB * * *' in line:
                break
            if nblank == 100:
                log.warning('breaking on nblank=100 b')
                break

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
            else:
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
            else:
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
        flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                  mach, density_ratio, method,
                                  modes, results,
                                  f06_units=f06_units, out_units=out_units,
                                  make_alt=make_alt)
        flutters[subcase] = flutter
    return flutters

def _get_units(units):
    # type: (Optional[Union[str, Dict[str, str]]]) -> Optional[Union[str, Dict[str, str]]]
    """gets the units"""
    if units is None:
        units = 'english_in'
        #units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                 #'altitude' : 'ft', 'dynamic_pressure' : 'psf', 'eas':'ft/s'}

    if isinstance(units, str):
        units = units.lower()
        if units == 'si':
            units = {'velocity' : 'm/s', 'density' : 'kg/m^3',
                     'altitude' : 'm', 'dynamic_pressure' : 'Pa', 'eas':'m/s'}
        elif units == 'english_in':
            units = {'velocity' : 'in/s', 'density' : 'slinch/in^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psi', 'eas':'in/s'}
        elif units == 'english_ft':
            units = {'velocity' : 'ft/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf', 'eas':'ft/s'}
        elif units == 'english_kt':
            units = {'velocity' : 'knots', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf', 'eas':'knots'}
        else:
            raise NotImplementedError('units=%r must be in [si, english_in, '
                                      'english_ft, english_kt]' % units)
    else:
        assert isinstance(units, dict), 'units=%r' % (units)
    return units


def plot_flutter_f06(f06_filename, f06_units=None, out_units=None, make_alt=False,
                     plot_type='tas', modes=None,
                     plot_vg=False, plot_vg_vf=False, plot_root_locus=False,
                     plot_kfreq_damping=False,
                     xlim=None, ylim_damping=None, ylim_freq=None, ylim_kfreq=None,
                     nopoints=False, noline=False,
                     export_zona_filename=None,
                     export_veas_filename=None,
                     export_f06_filename=None,
                     vg_filename=None,
                     vg_vf_filename=None,
                     root_locus_filename=None,
                     kfreq_damping_filename=None,
                     plot=True, show=True, clear=False, close=False,
                     log=None):
    """
    Plots a flutter (SOL 145) deck

    Parameters
    ----------
    f06_filename : str
        the filename to plot
    f06_units : Dict[name]=unit; default=None
        f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
    out_units : Dict[name]=unit; default=None
        out_units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf'}
    modes : bool; default=None
        specifies the modes to plot
    plot_type : str; default='tas'
        'tas' : true airspeed
        'eas' : equivalent airspeed
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
    flutters = make_flutter_response(
        f06_filename, f06_units=f06_units, out_units=out_units, make_alt=make_alt, log=log)

    if plot:
        make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                           plot_type,
                           plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                           nopoints, noline,
                           export_zona_filename=export_zona_filename,
                           export_veas_filename=export_veas_filename,
                           export_f06_filename=export_f06_filename,
                           vg_filename=vg_filename,
                           vg_vf_filename=vg_vf_filename,
                           root_locus_filename=root_locus_filename,
                           kfreq_damping_filename=kfreq_damping_filename,
                           show=show, clear=clear, close=close)
    return flutters

def make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq, ylim_kfreq,
                       plot_type,
                       plot_vg: bool,
                       plot_vg_vf: bool,
                       plot_root_locus: bool,
                       plot_kfreq_damping: bool,
                       nopoints: bool,
                       noline: bool,
                       legend: bool=True,
                       export_zona_filename: Optional[str]=None,
                       export_veas_filename: Optional[str]=None,
                       export_f06_filename: Optional[str]=None,
                       vg_filename: Optional[str]=None,
                       vg_vf_filename: Optional[str]=None,
                       root_locus_filename: Optional[str]=None,
                       kfreq_damping_filename: Optional[str]=None,
                       show: bool=True, clear: bool=False, close: bool=False):
    """actually makes the flutter plots"""
    for subcase, flutter in sorted(flutters.items()):
        if plot_vg:
            filenamei = None if vg_filename is None else vg_filename % subcase
            flutter.plot_vg(modes=modes,
                            plot_type=plot_type,
                            xlim=xlim, ylim_damping=ylim_damping,
                            png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_vg_vf:
            filenamei = None if vg_vf_filename is None else vg_vf_filename % subcase
            flutter.plot_vg_vf(modes=modes,
                               plot_type=plot_type,
                               xlim=xlim,
                               ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                               nopoints=nopoints, noline=noline,
                               legend=legend,
                               png_filename=filenamei, show=False, clear=clear, close=close)
        if plot_root_locus:
            filenamei = None if root_locus_filename is None else root_locus_filename % subcase
            flutter.plot_root_locus(modes=modes,
                                    fig=None, axes=None,
                                    xlim=None, ylim=None,
                                    clear=clear, legend=True,
                                    png_filename=filenamei,
                                    show=False, close=close)

        if plot_kfreq_damping:
            filenamei = None if kfreq_damping_filename is None else kfreq_damping_filename % subcase
            flutter.plot_kfreq_damping(modes=modes,
                                       plot_type=plot_type,
                                       ylim_damping=ylim_damping,
                                       ylim_kfreq=ylim_kfreq,
                                       nopoints=nopoints, noline=noline,
                                       png_filename=filenamei, show=False, clear=clear, close=close)
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

if __name__ == '__main__':  # pragma: no cover
    plot_flutter_f06('bah_plane.f06')
