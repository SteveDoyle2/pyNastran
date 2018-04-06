"""
SOL 145 plotter
"""
from __future__ import print_function
from six import iteritems

import matplotlib.pyplot as plt

from pyNastran.utils.log import get_logger2
from pyNastran.f06.flutter_response import FlutterResponse


def make_flutter_response(f06_filename, f06_units=None, out_units=None, log=None):
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
                #print('subcasei = %r' % new_subcase)
                if new_subcase > subcase:
                    log.debug('subcase=%s -> new_subcase=%s' % (subcase, new_subcase))
                    log.debug('modes1 = %s' % modes)
                    flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                              mach, density_ratio, method,
                                              modes, results,
                                              f06_units=f06_units, out_units=out_units)
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

                iline += 1
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    log.warning('breaking on nblank=100 a')
                    break
            if nblank == 100:
                log.warning('breaking on nblank=100 b')
                break

            # pulls the subcase id for the first subcase
            if last_line is not None:
                #log.debug('line%i_FSb = %r' % (iline, line))
                #log.debug('line%i_FSb = %r' % (iline-1, last_line.replace('     ', ' ')))
                sline = last_line.strip().split()
                isubcase = sline.index('SUBCASE')
                subcase = int(sline[isubcase + 1])
                log.debug('subcase = %s' % subcase)

            configuration_sline = f06_file.readline().split()
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
            else:
                raise NotImplementedError(point_sline)

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
                    'EIGENVALUE' not in sline)
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
                                  f06_units=f06_units, out_units=out_units)
        flutters[subcase] = flutter
    return flutters

def _get_units(units):
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
            raise NotImplementedError('units=%r' % units)
    else:
        assert isinstance(units, dict), 'units=%r' % (units)
    return units


def plot_flutter_f06(f06_filename, f06_units=None, out_units=None,
                     plot_type='tas', modes=None,
                     plot_vg=False, plot_vg_vf=False, plot_root_locus=False,
                     plot_kfreq_damping=False, show=True,
                     xlim=None, ylim_damping=None, ylim_freq=None, nopoints=False,
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

    Returns
    -------
    flutters : dict
        key : int
           subcase_id
        value : FlutterResponse()

    Supports:
    ---------
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
    ----------------
     o SOL 200
     o fixing mode switching problem
     o fixing unconverged points
    """
    flutters = make_flutter_response(
        f06_filename, f06_units=f06_units, out_units=out_units, log=log)

    make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq,
                       plot_type,
                       plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                       nopoints,
                       show=show)
    return flutters

def make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq,
                       plot_type,
                       plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                       nopoints,
                       show=True):
    """actually makes the flutter plots"""
    for unused_subcase, flutter in sorted(iteritems(flutters)):
        if plot_vg:
            flutter.plot_vg(modes=modes,
                            show=False,
                            xlim=xlim, ylim=ylim_damping)
        if plot_vg_vf:
            flutter.plot_vg_vf(modes=modes,
                               plot_type=plot_type,
                               show=False,
                               xlim=xlim,
                               ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                               nopoints=nopoints)
        if plot_root_locus:
            flutter.plot_root_locus(modes=modes, show=False)
        if plot_kfreq_damping:
            flutter.plot_kfreq_damping(modes=modes,
                                       ylim_damping=ylim_damping,
                                       ylim_kfreq=None,
                                       show=False)
    if show:
        plt.show()

if __name__ == '__main__':  # pragma: no cover
    plot_flutter_f06('bah_plane.f06')
