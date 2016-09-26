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
                 f06_units=None, out_units=None, plot_units=None,
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

        plot_units dict[str] = str (default=None -> no units conversion)
            The units to use when plotting.
            TODO: Remove this???

            PK method:
                plot_units = {'velocity' : 'ft/s'}
            PKNL method:
                plot_units = {'velocity' : 'ft/s', 'eas' : 'knots', 'density': 'kg/m^3'}

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
        self.plot_units = plot_units

        self.subcase = subcase
        self.configuration = configuration
        self.xysym = xysym
        self.xzsym = xzsym
        self.method = method
        self.modes = np.asarray(modes, dtype='int32')
        self.density_ratio = density_ratio
        rho_ref = 1.

        self.ikfreq = 0
        self.iinvkfreq = 1

        #print(results)
        results = np.asarray(results, dtype='float64')
        if self.method == 'PK':
            self.ivelocity = 2
            self.idamping = 3
            self.ifreq = 4
            self.ieigr = 5
            self.ieigi = 6
            self.results = results

        elif self.method == 'PKNL':
            raise NotImplementedError(method)
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

            density_units = self.out_units['density']
            altitude_units = self.out_units['altitude']

            kdensity = self._get_altitude_unit_factor(density_units, 'slug/ft^3')
            kalt = self._get_altitude_unit_factor('ft', altitude_units)
            kpressure = self._get_unit_factor('dynamic_pressure')[0]

            alt = np.array(
                [get_alt_for_density(densityi, nmax=20) * kalt
                 for densityi in rho.ravel() * kdensity], dtype='float64').reshape(vel.shape)

            # TODO: this is wrong...
            self.results = np.hstack([results, eas, q * kpressure, alt])
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
        unit_out = self.out_units[name]
        unit_plot = self.plot_units[name]

        #print('unit_out=%r unit_out=%r' % (unit_out, unit_plot))
        if name == 'velocity':
            factor = self._get_velocity_unit_factor(unit_out, unit_plot)
        elif name == 'altitude':
            factor = self._get_altitude_unit_factor(unit_out, unit_plot)
        elif name == 'density':
            factor = self._get_altitude_unit_factor(unit_out, unit_plot)

        elif name in ['pressure', 'dynamic_pressure']:
            factor = self._get_pressure_unit_factor(unit_out, unit_plot)
        else:
            raise NotImplementedError(name)

        if self.plot_units is not None:
            units = self.plot_units[name]
        else:
            units = 'units'

        return factor, units

    @staticmethod
    def _get_velocity_unit_factor(unit_in, unit_out):
        """TODO: simplify this..."""
        if unit_in not in ['in/s', 'ft/s', 'knots', 'm/s']:
            msg = 'unit_out=%r not in [in/s, ft/s, knots, m/s]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['in/s', 'ft/s', 'knots', 'm/s']:
            msg = 'unit_plot=%r not in [in/s, ft/s, knots, m/s]' % unit_out
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
        if unit_out not in ['ft', 'm']:
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
        if unit_out not in ['psi', 'psf']:
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
        if unit_out not in ['slinch/in^3', 'slug/ft^3']:
            msg = 'unit_in=%r not in [slinch/in^3, slug/ft^3]' % unit_in
            raise NotImplementedError(msg)
        if unit_out not in ['ft', 'm']:
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

        kvelocity, velocity_units = self._get_unit_factor('velocity')
        for i, imode, mode in zip(count(), imodes, modes):
            vel = self.results[imode, :, self.ivelocity].ravel() * kvelocity
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
            plt.figure()

        self._set_xy_limits(xlim, ylim)
        imodes = np.searchsorted(self.modes, modes)
        symbols = self.symbols

        for i, imode, mode in zip(count(), imodes, modes):
            real = self.results[imode, :, self.ieigr].ravel()
            imag = self.results[imode, :, self.ieigi].ravel()
            plt.plot(real, imag, symbols[i], label='Mode %i' % mode)

        plt.grid(True)
        plt.xlabel('Eigenvalue (Real)')
        plt.ylabel('Eigenvalue (Imaginary)')

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        plt.suptitle(title)
        self._finalize_plot(legend, show, png_filename, clear)

    def plot_vg_vf(self, fig=None, modes=None, show=None, png_filename=None,
                   clear=False, legend=None):
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
            fig = plt.figure(self.subcase) # figsize=(12,9)
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        imodes = np.searchsorted(self.modes, modes)
        symbols = self.symbols

        kvelocity, velocity_units = self._get_unit_factor('velocity')
        for i, imode, mode in zip(count(), imodes, modes):
            vel = self.results[imode, :, self.ivelocity].ravel() * kvelocity
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()
            damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            freq_axes.plot(vel, freq, symbols[i])

        damp_axes.set_xlabel('Velocity [%s]' % velocity_units)
        damp_axes.set_ylabel('Damping')
        damp_axes.grid(True)

        freq_axes.set_xlabel('Velocity [%s]' % velocity_units)
        freq_axes.set_ylabel('Frequency [Hz]')
        freq_axes.grid(True)

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename

        damp_axes.set_title(title)
        #plt.suptitle(title)

        xlim_damp = None
        ylim_damp = None
        ylim_freq = None
        if xlim_damp:
            damp_axes.set_xlim(xlim_damp)
        if ylim_damp:
            damp_axes.set_ylim(ylim_damp)
        if ylim_freq:
            freq_axes.set_ylim(ylim_freq)

        #self._finalize_plot(legend, show, png_filename, clear)
        if legend:
            #damp_axes.legend(legend_items, fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2)
            #fig.subplots_adjust(hspace=0.25)
            #fig.subplots_adjust(hspace=.5)
            plt.legend()

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


def plot_flutter_f06(f06_filename, modes=None,
                     plot_vg=False, plot_vg_vf=True, plot_root_locus=True,
                     show=True):
    """
    TODO: support multiple subcases
    TODO: support long tables
    """
    f06_units = {'velocity' : 'in/s'}
    out_units = {'velocity' : 'in/s'}
    plot_units = {'velocity' : 'ft/s'}

    flutters = []
    with open(f06_filename, 'r') as f06_file:
        subcase = 1
        results = []
        modes = []
        while 1:
            nblank = 0
            line = f06_file.readline()
            #print('line = %s' % line.strip())
            while 'SUBCASE' not in line:
                line = f06_file.readline()
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    break
            if nblank == 100:
                break
            sline = line.strip().split()
            isubcase = sline.index('SUBCASE')
            new_subcase = int(sline[isubcase + 1])
            if new_subcase > subcase:
                print('new_subcase = ', subcase)
                break

            while 'FLUTTER  SUMMARY' not in line:
                line = f06_file.readline()
                if not line:
                    nblank += 1
                if nblank == 100:
                    print(line.strip())
                    break
            if nblank == 100:
                break

            configuration_sline = f06_file.readline().split()
            configuration = configuration_sline[2]
            xysym = configuration_sline[5]
            xzsym = configuration_sline[8]
            #print(configuration, xysym, xzsym)

            point_sline = f06_file.readline().split()
            mode = point_sline[2]
            mach = point_sline[6]
            density_ratio = point_sline[10]
            method = point_sline[13]
            print(mode, mach, density_ratio, method)
            if mode in modes:
                print('found existing mode...')
                continue
            modes.append(mode)

            # blanks
            f06_file.readline()
            f06_file.readline()

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
                if 'PAGE' not in sline and 'INFORMATION' not in sline and 'EIGENVALUE' not in sline:
                    lines.append(sline)
                    #print(sline)

            results.append(lines)
            #print('')

        #print(len(results))
        flutter = FlutterResponse(subcase, configuration, xysym, xzsym,
                                  mach, density_ratio, method,
                                  modes, results,
                                  f06_units=f06_units, out_units=out_units, plot_units=plot_units)
        flutters.append(flutter)

        if plot_vg:
            flutter.plot_vg(show=False)
        if plot_vg_vf:
            flutter.plot_vg_vf(show=False)
        if plot_root_locus:
            flutter.plot_root_locus(show=False)
        if show:
            plt.show()
    return flutters

if __name__ == '__main__':
    plot_flutter_f06('bah_plane.f06')
