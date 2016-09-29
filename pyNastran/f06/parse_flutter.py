"""
SOL 145 plotter
"""
from __future__ import print_function
from itertools import count
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pyNastran.utils.atmosphere import get_alt_for_density


class FlutterResponse(object):
    """storage object for single subcase SOL 145 results"""
    def __init__(self, subcase, configuration, xysym, xzsym, mach, density_ratio, method,
                 modes, results,
                 f06_units=None, out_units=None,
                ):
        """
        Parameters
        ----------
        subcase : int
            the subcase id
        method : str
            PK, PKNL, ???
        modes : List[int]; (default=None -> all)
            the modes; typically 1 to N
        results : varies
            method = PK
                List[List[float] * 7] * nmodes
                kfreq, 1/kfreq,                velocity, damping, freq, eigr, eigi
            method = PKNL
                List[List[float] * 9] * nmodes
                kfreq, 1/kfreq, density, mach, velocity, damping, freq, eigr, eigi

        Units
        -----
        Units are only applicable for quantities that have units
        (e.g. units on Mach don't do anything).
        All allowable fields are shown below.

        f06_units : dict[str] = str (default=None -> no units conversion)
            The units to to read from the F06.
            PK method:
                f06_units = {'velocity' : 'in/s'}
                The velocity units are the units for the FLFACT card in the BDF
            PKNL method:
                f06_units = {'velocity' : 'in/s', 'density' : 'slug/in^3'}
                The velocity/density units are the units for the FLFACT card in the BDF

        out_units dict[str] = str (default=None -> no units conversion)
            The units to store results in.
            PK method:
                out_units = {'velocity' : 'ft/s'}
            PKNL method:
                out_units = {
                    'velocity' : 'ft/s',
                    'eas' : 'knots',
                    'density' : 'slinch/in^3',
                    'altitude' : 'ft'
                }

        Unused Parameters
        -----------------
        configuration : str
            AEROSG2D...what is this???
        XY-SYMMETRY : str
            ASYMMETRIC, SYMMETRIC
            unused
        XZ-SYMMETRY : str
            ASYMMETRIC, SYMMETRIC
            unused
        """
        self.f06_units = f06_units
        self.out_units = out_units

        self.subcase = subcase
        self.configuration = configuration
        if method == 'PK':
            self.mach = mach
            self.xysym = xysym
            self.xzsym = xzsym
            self.density_ratio = density_ratio
            print('mach=%s' % mach)

        self.method = method
        self.modes = np.asarray(modes, dtype='int32')
        rho_ref = 1.

        self.ikfreq = 0
        self.ikfreq_inv = 1

        #print(results)
        results = np.asarray(results, dtype='float64')
        if self.method == 'PK':
            self.ivelocity = 2
            self.idamping = 3
            self.ifreq = 4
            self.ieigr = 5
            self.ieigi = 6
            self.results = results

            kvel = self._get_unit_factor('velocity')[0]
            results[:, :, self.ivelocity] *= kvel

        elif self.method == 'PKNL':
            self.idensity = 2
            self.imach = 3
            self.ivelocity = 4
            self.idamping = 5
            self.ifreq = 6
            self.ieigr = 7
            self.ieigi = 8
            self.ieas = 9
            self.iq = 10
            self.ialt = 11

            # eas = V * sqrt(rho / rhoSL)
            vel = results[:, :, self.ivelocity]#.ravel()
            rho = results[:, :, self.idensity]#.ravel()

            q = 0.5 * rho * vel**2
            #eas  = (2 * q / rho_ref)**0.5
            eas = vel * np.sqrt(rho / rho_ref)

            density_units1 = self.f06_units['density']
            altitude_units = self.out_units['altitude']
            #density_units1 = self.out_units['density']

            kdensity = self._get_density_unit_factor(density_units1, 'slug/ft^3')
            #kdensity = self._get_altitude_unit_factor(density_units2, 'slug/ft^3')
            kalt = self._get_altitude_unit_factor('ft', altitude_units)
            kvel = self._get_unit_factor('velocity')[0]
            #kpressure = self._get_unit_factor('dynamic_pressure')[0]
            kpressure = kdensity * kvel ** 2

            make_alt = False
            if make_alt:
                alt = np.array(
                    [get_alt_for_density(densityi, nmax=20) * kalt
                     for densityi in rho.ravel() * kdensity], dtype='float64').reshape(vel.shape)

                self.results = np.dstack([results, eas, q * kpressure, alt])
            else:
                self.results = np.dstack([results, eas, q * kpressure])
        else:
            raise NotImplementedError(method)
        #print(self.results.shape)

        colors = ['r', 'g', 'b', 'k']
        symbols = ['o', '*', 'x', 'v', '>', '<', '^']
        self._symbols = []
        for symbol in symbols:
            for color in colors:
                self._symbols.append(color + symbol)
        self.noline = False

    def set_plot_options(self, noline=False):
        self.noline = noline

    def _get_unit_factor(self, name):
        if not self.f06_units or not self.out_units:
            return 1.
        unit_f06 = self.f06_units[name]
        unit_out = self.out_units[name]

        #print('unit_f06=%r unit_out=%r' % (unit_f06, unit_out))
        if name == 'velocity':
            factor = self._get_velocity_unit_factor(unit_f06, unit_out)
        elif name == 'altitude':
            factor = self._get_altitude_unit_factor(unit_f06, unit_out)
        elif name == 'density':
            factor = self._get_altitude_unit_factor(unit_f06, unit_out)
        elif name in ['pressure', 'dynamic_pressure']:
            factor = self._get_pressure_unit_factor(unit_f06, unit_out)
        else:
            raise NotImplementedError(name)

        if self.out_units is not None:
            units = self.out_units[name]
        else:
            units = 'units'
        return factor, units

    @staticmethod
    def _get_velocity_unit_factor(unit_in, unit_out):
        """TODO: simplify this..."""
        if unit_in not in ['in/s', 'ft/s', 'knots', 'm/s']:
            msg = 'unit_in=%r not in [in/s, ft/s, knots, m/s]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['in/s', 'ft/s', 'knots', 'm/s']:
            msg = 'unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_out
            raise NotImplementedError(msg)

        if unit_in == unit_out:
            factor = 1.
        elif unit_in == 'in/s':
            if unit_out == 'ft/s':
                factor = 1./12.
            elif unit_out == 'knots':
                factor = 1./20.2537
            elif unit_out == 'm/s':
                factor = 1./39.3700432322835
            else:
                raise NotImplementedError('unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_out)
        elif unit_in == 'ft/s':
            if unit_out == 'in/s':
                factor = 12.
            elif unit_out == 'knots':
                factor = 1./1.687808333407337269
            elif unit_out == 'm/s':
                factor = 1./3.2808369360236246948
            else:
                raise NotImplementedError('unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_out)
        elif unit_in == 'knots':
            if unit_out == 'in/s':
                factor = 20.253700000888049004
            elif unit_out == 'ft/s':
                factor = 1.687808333407337269
            elif unit_out == 'm/s':
                factor = 1/1.9438427409666045875
            else:
                raise NotImplementedError('unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_out)
        elif unit_in == 'm/s':
            if unit_out == 'in/s':
                factor = 39.37004323228349989
            elif unit_out == 'ft/s':
                factor = 3.2808369360236246948
            elif unit_out == 'knots':
                factor = 1.9438427409666045875
            else:
                raise NotImplementedError('unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_out)
        else:
            raise NotImplementedError('unit_in=%r not in [in/s, ft/s, knots, m/s]' % unit_in)
        return factor

    @staticmethod
    def _get_altitude_unit_factor(unit_in, unit_out):
        """TODO: simplify this..."""
        if unit_in not in ['ft', 'm']:
            msg = 'unit_in=%r not in [ft, m]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['ft', 'm']:
            msg = 'unit_out=%r not in [ft, m]' % unit_out
            raise NotImplementedError(msg)

        if unit_in == unit_out:
            factor = 1.
        elif unit_in == 'm':
            if unit_out == 'ft':
                factor = 0.3048
            else:
                raise NotImplementedError('unit_out=%r not in [ft, m]' % unit_out)
        elif unit_in == 'ft':
            if unit_out == 'm':
                factor = 1./0.3048
            else:
                raise NotImplementedError('unit_out=%r not in [ft, m]' % unit_out)
        else:
            raise NotImplementedError('unit_in=%r not in [ft, m]' % unit_in)
        return factor

    @staticmethod
    def _get_pressure_unit_factor(unit_in, unit_out):
        """TODO: simplify this..."""
        if unit_in not in ['psi', 'psf']:
            msg = 'unit_in=%r not in [psi, psf]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['psi', 'psf']:
            msg = 'unit_out=%r not in [psi, psf]' % unit_out
            raise NotImplementedError(msg)

        if unit_in == unit_out:
            factor = 1.
        elif unit_in == 'psi':
            if unit_out == 'psf':
                factor = 1./144.
            else:
                raise NotImplementedError('unit_out=%r not in [psi, psf]' % unit_out)
        elif unit_in == 'psi':
            if unit_out == 'psf':
                factor = 144.
            else:
                raise NotImplementedError('unit_out=%r not in [psi, psf]' % unit_out)
        else:
            raise NotImplementedError('unit_in=%r not in [psi, psf]' % unit_in)
        return factor

    @staticmethod
    def _get_density_unit_factor(unit_in, unit_out):
        """TODO: simplify this..."""
        if unit_in not in ['slinch/in^3', 'slug/ft^3']:
            msg = 'unit_in=%r not in [slinch/in^3, slug/ft^3]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['slinch/in^3', 'slug/ft^3']:
            msg = 'unit_out=%r not in [slinch/in^3, slug/ft^3]' % unit_out
            raise NotImplementedError(msg)

        if unit_in == unit_out:
            factor = 1.
        elif unit_in == 'slinch/in^3':
            if unit_out == 'slug/ft^3':
                factor = 12.**4.
            else:
                msg = 'unit_out=%r not in [slinch/in^3, slug/ft^3]' % unit_out
                raise NotImplementedError(msg)
        elif unit_in == 'slug/ft^3':
            if unit_out == 'slinch/in^3':
                factor = 1./12.**4.
            else:
                msg = 'unit_out=%r not in [slinch/in^3, slug/ft^3]' % unit_out
                raise NotImplementedError(msg)
        else:
            msg = 'unit_in=%r not in [slinch/in^3, slug/ft^3]' % unit_in
            raise NotImplementedError(msg)
        return factor

    @property
    def symbols(self):
        """gets the symbols for the lines"""
        if not self.noline:
            symbols = [symbol + '-' for symbol in self._symbols]
        else:
            symbols = self._symbols
        return symbols

    @staticmethod
    def _set_xy_limits(xlim, ylim):
        """helper method for ``plot_vg``"""
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)

    def plot_vg(self, modes=None,
                fig=None,
                xlim=None, ylim=None,
                show=True, clear=False, legend=True,
                png_filename=None):
        """
        Make a V-g plot

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        """
        if modes is None:
            modes = self.modes
        else:
            modes = np.asarray(modes)
        if fig is None:
            plt.figure(self.subcase)

        self._set_xy_limits(xlim, ylim)
        imodes = np.searchsorted(self.modes, modes)
        symbols = self.symbols

        _kvelocity, velocity_units = self._get_unit_factor('velocity')
        for i, imode, mode in zip(count(), imodes, modes):
            vel = self.results[imode, :, self.ivelocity].ravel()
            damping = self.results[imode, :, self.idamping].ravel()
            plt.plot(vel, damping, symbols[i], label='Mode %i' % mode)

        plt.grid(True)
        plt.xlabel('Velocity [%s]' % velocity_units)

        plt.ylabel('Damping')

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        plt.suptitle(title)
        self._finalize_plot(legend, show, png_filename, clear)

    def plot_root_locus(self, modes=None,
                        fig=None,
                        xlim=None, ylim=None,
                        show=True, clear=False, legend=True,
                        png_filename=None):
        """
        Plots a root locus

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        """
        if modes is None:
            modes = self.modes
        else:
            modes = np.asarray(modes)

        if fig is None:
            fig = plt.figure()
        axes = fig.add_subplot(111)

        self._set_xy_limits(xlim, ylim)
        imodes = np.searchsorted(self.modes, modes)
        symbols = self.symbols

        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            real = self.results[imode, :, self.ieigr].ravel()
            imag = self.results[imode, :, self.ieigi].ravel()
            axes.plot(real, imag, symbol, label='Mode %i' % mode, markersize=0)

            s = np.linspace(.75, 50., len(real))
            #assert symbol[2] == '-', symbol
            axes.scatter(real, imag, s=s, color=symbol[0], marker=symbol[1])

        axes.grid(True)
        axes.set_xlabel('Eigenvalue (Real)')
        axes.set_ylabel('Eigenvalue (Imaginary)')

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        #self._finalize_plot(legend, show, png_filename, clear)
        if legend:
            axes.legend()

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            plt.clear()

    def plot_vg_vf(self, fig=None, modes=None, show=None, png_filename=None,
                   clear=False, legend=None,
                   xlim=None, ylim_damping=None, ylim_freq=None):
        """
        Make a V-g and V-f plot

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        """
        if modes is None:
            modes = self.modes
        else:
            modes = np.asarray(modes)
        if fig is None:
            fig = plt.figure() # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        imodes = np.searchsorted(self.modes, modes)
        symbols = self.symbols

        _kvelocity, velocity_units = self._get_unit_factor('velocity')

        legend_items = ['Mode %i' % mode for mode in modes]
        for i, imode, mode in zip(count(), imodes, modes):
            vel = self.results[imode, :, self.ivelocity].ravel()
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()
            damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            freq_axes.plot(vel, freq, symbols[i])

        damp_axes.set_xlabel('Velocity [%s]' % velocity_units)
        damp_axes.set_ylabel('Damping')
        damp_axes.grid(True)
        if xlim is not None:
            damp_axes.set_xlim(xlim)
        if ylim_damping is not None:
            damp_axes.set_ylim(ylim_damping)

        freq_axes.set_xlabel('Velocity [%s]' % velocity_units)
        freq_axes.set_ylabel('Frequency [Hz]')
        freq_axes.grid(True)

        if xlim is not None:
            freq_axes.set_xlim(xlim)
        if ylim_freq is not None:
            freq_axes.set_ylim(ylim_freq)

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename

        damp_axes.set_title(title)
        #plt.suptitle(title)

#self._finalize_plot(legend, show, png_filename, clear)
        if legend:
            damp_axes.legend(legend_items, fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2)
            #fig.subplots_adjust(hspace=0.25)
            #fig.subplots_adjust(hspace=.5)
            #plt.legend()
            #damp_axes.legend(legend_items, bbox_to_anchor=anchor, ncol=2)
            #fig.subplots_adjust(hspace=0.25)
            #fig.subplots_adjust(hspace=.5)

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            plt.clear()

    @staticmethod
    def _finalize_plot(legend, show, png_filename, clear):
        """common helper method"""
        if legend:
            plt.legend()

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            plt.clear()
        #return fig


def plot_flutter_f06(f06_filename, f06_units=None, out_units=None,
                     modes=None,
                     plot_vg=False, plot_vg_vf=True, plot_root_locus=True,
                     show=True,
                     xlim=None, ylim_damping=None, ylim_freq=None):
    """
    Plots a flutter (SOL 145) deck
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

    flutters = []
    iline = 0
    with open(f06_filename, 'r') as f06_file:
        subcase = 1
        results = []
        modes = []

        while 1:
            nblank = 0
            line = f06_file.readline()
            iline += 1
            #print('line%ia = %r' % (iline, line))
            while 'SUBCASE ' not in line and 'FLUTTER  SUMMARY' not in line:
                line = f06_file.readline()
                iline += 1
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    break
            if nblank == 100:
                break

            #print('line%ib = %r' % (iline, line))
            if 'SUBCASE' in line[109:]:
                sline = line.strip().split()
                isubcase = sline.index('SUBCASE')
                new_subcase = int(sline[isubcase + 1])
                #print('subcasei = %r' % new_subcase)
                if new_subcase > subcase:
                    print('\nsubcase=%s -> new_subcase=%s' % (subcase, new_subcase))
                    print('modes1 =', modes)
                    flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                              mach, density_ratio, method,
                                              modes, results,
                                              f06_units=f06_units, out_units=out_units)
                    flutters.append(flutter)
                    modes = []
                    results = []

                    subcase = new_subcase
                    #break
                continue

            #print('line%i_FS = %r' % (iline, line))
            while 'FLUTTER  SUMMARY' not in line:
                line = f06_file.readline()
                iline += 1
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    break
            if nblank == 100:
                break

            configuration_sline = f06_file.readline().split()
            iline += 1
            configuration = configuration_sline[2]
            xysym = configuration_sline[5]
            xzsym = configuration_sline[8]
            #print(configuration, xysym, xzsym)

            point_sline = f06_file.readline().split()
            iline += 1
            mode = int(point_sline[2])
            method = point_sline[-1]  # 13 for PN, 5 for PK
            #print(point_sline)
            if method == 'PK':
                mach = float(point_sline[6])
                density_ratio = float(point_sline[10])
                #method = point_sline[13]
                print(iline, mode, mach, density_ratio, method)
            elif method == 'PKNL':
                mach = None
                density_ratio = None
                print(iline, mode, method)
                f06_file.readline()
                iline += 1
            else:
                raise NotImplementedError(point_sline)
            if 1:
                if mode in modes:
                    print('found existing mode...')
                    continue
            modes.append(mode)

            # blanks
            f06_file.readline()
            f06_file.readline()
            iline += 2

            lines = []

            #KFREQ           1./KFREQ                                 VELOCITY       DAMPING     FREQUENCY      COMPLEX   EIGENVALUE - PK
            # KFREQ          1./KFREQ       DENSITY     MACH NO.      VELOCITY       DAMPING     FREQUENCY      COMPLEX   EIGENVALUE - PKNL
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

            results.append(lines)
            #print('')

        #print(len(results))
        print('modes =', modes)
        flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                  mach, density_ratio, method,
                                  modes, results,
                                  f06_units=f06_units, out_units=out_units)
        flutters.append(flutter)

        for flutter in flutters:
            if plot_vg:
                flutter.plot_vg(show=False,
                                xlim=xlim, ylim=ylim_damping)
            if plot_vg_vf:
                flutter.plot_vg_vf(show=False,
                                   xlim=xlim,
                                   ylim_damping=ylim_damping, ylim_freq=ylim_freq)
            if plot_root_locus:
                flutter.plot_root_locus(show=False)
        if show:
            plt.show()
    return flutters

if __name__ == '__main__':
    plot_flutter_f06('bah_plane.f06')
