"""
SOL 145 plotter
"""
from __future__ import print_function
from six import iteritems

import matplotlib.pyplot as plt

from pyNastran.utils.log import get_logger2
from pyNastran.f06.flutter_response import FlutterResponse

def plot_flutter_f06(f06_filename, f06_units=None, out_units=None,
                     modes=None,
                     plot_vg=False, plot_vg_vf=False, plot_root_locus=False,
                     plot_kfreq_damping=False, show=True,
                     xlim=None, ylim_damping=None, ylim_freq=None):
    """
    Plots a flutter (SOL 145) deck

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
     o long tables (use LINE=500000)
     o SOL 200
     o fixing mode switching problem
     o fixing unconverged points
    """
    if f06_units is None:
        f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
    if out_units is None:
        out_units = {'velocity' : 'in/s', 'density' : 'slug/ft^3',
                     'altitude' : 'ft', 'dynamic_pressure' : 'psf'}

    log = get_logger2(log=None, debug=True, encoding='utf-8')
    flutters = {}
    iline = 0
    modes_to_plot = modes

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
                    log.info('subcase=%s -> new_subcase=%s' % (subcase, new_subcase))
                    log.info('modes1 = %s' % modes)
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
                log.info('subcase = %s' % subcase)

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
                if (sline
                    and 'PAGE' not in sline
                    and 'INFORMATION' not in sline
                    and 'EIGENVALUE' not in sline):
                    #print('sline = %s' % sline)
                    lines.append(sline)

            if found_existing_mode:
                results[mode-1].extend(lines)
            else:
                results.append(lines)
                modes.append(mode)
            #print('')

        log.info('modes = %s' % modes)
        flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                  mach, density_ratio, method,
                                  modes, results,
                                  f06_units=f06_units, out_units=out_units)
        flutters[subcase] = flutter

    make_flutter_plots(modes_to_plot, flutters, xlim, ylim_damping, ylim_freq,
                       plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                       show=show)
    return flutters

def make_flutter_plots(modes, flutters, xlim, ylim_damping, ylim_freq,
                       plot_vg, plot_vg_vf, plot_root_locus, plot_kfreq_damping,
                       show=True):
    """actually makes the flutter plots"""
    for subcase, flutter in sorted(iteritems(flutters)):
        if plot_vg:
            flutter.plot_vg(modes=modes,
                            show=False,
                            xlim=xlim, ylim=ylim_damping)
        if plot_vg_vf:
            flutter.plot_vg_vf(modes=modes,
                               show=False,
                               xlim=xlim,
                               ylim_damping=ylim_damping, ylim_freq=ylim_freq)
        if plot_root_locus:
            flutter.plot_root_locus(modes=modes, show=False)
        if plot_kfreq_damping:
            flutter.plot_kfreq_damping(modes=modes, show=False)
    if show:
        plt.show()

if __name__ == '__main__':  # pragma: no cover
    plot_flutter_f06('bah_plane.f06')
