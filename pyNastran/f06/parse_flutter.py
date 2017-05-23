"""
SOL 145 plotter
"""
from __future__ import print_function
from itertools import count
from six import iteritems
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pyNastran.utils.atmosphere import get_alt_for_density
from pyNastran.utils.log import get_logger2
from pyNastran.utils import object_attributes, object_methods

class FlutterResponse(object):
    """storage object for single subcase SOL 145 results"""

    def object_attributes(self, mode='public', keys_to_skip=None):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : List[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : List[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong
        """
        return object_methods(obj, mode=mode, keys_to_skip=keys_to_skip)

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
            #print('mach=%s' % mach)

        self.method = method
        self.modes = np.asarray(modes, dtype='int32')
        rho_ref = 1.

        self.ikfreq = 0
        self.ikfreq_inv = 1

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
            # velocity is the target
            self.names = ['kfreq', '1/kfreq', 'velocity', 'damping', 'freq', 'eigr', 'eigi']

        elif self.method == 'PKNL':
            # velocity is the target
            self.names = ['kfreq', '1/kfreq', 'density', 'velocity', 'damping', 'freq', 'eigr', 'eigi', 'eas', 'q', 'alt']
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


        # c - cyan
        # b - black
        # r - red
        # g - green
        # m - magenta
        # y - yellow
        #colors = ['b', 'c', 'g', 'k', 'm', 'r'] #, 'y']
        # D - wide diamond
        # h - hexagon
        # * - star
        # + - plus
        # 3 - 3 pointed star
        # o - circle
        # d - thin diamond
        # 1 - Y shape
        # s - square
        #shapes = ['D', 'h', '*', 's', 'd', '3', 'o', '1', '2', '4', 'x', '^', '<', '>'] # '+',
        #symbol_list = []
        #for shape in shapes:
            #for color in colors:
                #symbol_list.append('%s-%s' % (shape, color))
        self.noline = False
        self._symbols = []
        self.generate_symbols()

    def generate_symbols(self):
        """
        This symbol list is taken from a series of "good" colors (e.g. not yellow)
        and easily distinguishable shapes.  Far more combinations that is necessary
        is defined
        """
        colors = ['r', 'g', 'b', 'k']
        symbols = ['o', '*', 'x', 'v', '>', '<', '^']
        self._symbols = []
        for symbol in symbols:
            for color in colors:
                self._symbols.append(color + symbol)

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
                png_filename=None, **kwargs):
        """
        Make a V-g plot

        See ``plot_root_locus`` for arguments
        """
        _kvelocity, velocity_units = self._get_unit_factor('velocity')
        xlabel = 'Velocity [%s]' % velocity_units
        ylabel = 'Damping'
        ix = self.ivelocity
        iy = self.idamping

        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim,
                       show=show, clear=clear, legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    @property
    def flutter_speed(self, modes=None):
        """gets the flutter speed"""
        if modes is None:
            modes = self.modes
        else:
            modes = np.asarray(modes)

    def plot_root_locus(self, modes=None,
                        fig=None,
                        xlim=None, ylim=None,
                        show=True, clear=False, legend=True,
                        png_filename=None,
                        **kwargs):
        """
        Plots a root locus

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        fig : plt.figure
            figure object
        xlim : list[float/None, float/None]
            the x plot limits
        ylim : list[float/None, float/None]
            the y plot limits
        show : bool; default=True
            show the plot
        clear : bool; default=False
            clear the plot
        legend : bool; default=False
            show the legend

        kwargs : dict; default=None
           key : various matplotlib parameters
           value : depends

        Legend kwargs
        -------------
           loc : str
              'best'
           fancybox : bool; default=False
              makes the box look cool
           framealpha : float; 0.0 <= alpha <= 1.0
               0.0 - fully transparent
               1.0 - no transparency / opaque
        """
        xlabel = 'Eigenvalue (Real)'
        ylabel = 'Eigenvalue (Imaginary)'
        ix = self.ieigr
        iy = self.ieigi
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim,
                       show=show, clear=clear, legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    def _plot_x_y(self, ix, iy, xlabel, ylabel, scatter, modes=None,
                  fig=None,
                  xlim=None, ylim=None,
                  show=True, clear=False, legend=True,
                  png_filename=None,
                  **kwargs):
        """builds the plot"""
        self.fix()
        if kwargs is None:
            kwargs = {}

        modes, imodes = self._get_modes_imodes(modes)

        if fig is None:
            fig = plt.figure()
        axes = fig.add_subplot(111)

        self._set_xy_limits(xlim, ylim)
        symbols = self.symbols

        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            freq = self.results[imode, :, self.ifreq].ravel()
            real = self.results[imode, :, ix].ravel()
            imag = self.results[imode, :, iy].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            axes.plot(real[iplot], imag[iplot], symbol, label='Mode %i' % mode, markersize=0)

            if scatter:
                s = np.linspace(.75, 50., len(real))
                #assert symbol[2] == '-', symbol
                axes.scatter(real[iplot], imag[iplot], s=s, color=symbol[0], marker=symbol[1])

        axes.grid(True)
        axes.set_xlabel('Eigenvalue (Real)')
        axes.set_ylabel('Eigenvalue (Imaginary)')

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        if legend:
            axes.legend(**kwargs)

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            plt.clear()

    def plot_kfreq_damping(self, modes=None,
                           fig=None,
                           xlim=None, ylim=None,
                           show=True, clear=False, legend=True,
                           png_filename=None,
                           **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments
        """
        xlabel = 'KFreq'
        ylabel = 'Damping'
        ix = self.ikfreq
        iy = self.idamping
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim,
                       show=show,
                       clear=clear,
                       legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    def fix(self):
        """attempts to fix the mode switching"""
        print(self.names)

        # results[imode, ipoint, iresult]
        # 1. NaN all the invalid points
        freq = self.results[:, :, self.ifreq]
        iplot, jplot = np.where(freq == 0.0)
        self.results[iplot, jplot, :] = np.nan
        return

        #-----------------------------------------------------------------------
        # 2. sort the results based on velocity so we're going low to high
        #nmodes, npoints = self.results.shape[:2]

        #for imode in range(nmodes):
            #print(self.results[imode, :, self.ivelocity])
            #isort = np.argsort(self.results[imode, :, self.ivelocity])
            #self.results[imode, :, :] = self.results[imode, isort, :]

        #-----------------------------------------------------------------------
        # 3. sort the results based on damping, so we're going abs(high) to low
        #for ipoint in range(npoints):
            #isort = np.argsort(self.results[:, ipoint, self.idamping])
            #self.results[:, isort, :] = self.results[:, isort, :]

        # 4. find the critical mode
        # 5. ???

    def _get_modes_imodes(self, modes):
        """gets the index of the modes to plot"""
        if modes is None:
            modes = self.modes
        elif isinstance(modes, slice):
            start = modes.start
            if modes.stop is None:
                stop = len(self.modes) + 1
            stop = modes.stop
            step = modes.step
            modes = np.unique(range(start, stop, step))
        elif len(modes) == 0:
            raise RuntimeError('modes = %s' % modes)
        else:
            modes = np.unique(modes)
        assert 0 not in modes, modes

        if modes.max() > self.modes.max():
            imodes = np.where(modes <= self.modes.max())
            modes = modes[imodes]
        if len(modes) == 0:
            raise RuntimeError('No modes to plot...')
        imodes = np.searchsorted(self.modes, modes)
        return modes, imodes

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
        self.fix()
        if fig is None:
            fig = plt.figure() # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        modes, imodes = self._get_modes_imodes(modes)
        symbols = self.symbols

        _kvelocity, velocity_units = self._get_unit_factor('velocity')

        legend_items = ['Mode %i' % mode for mode in modes]
        for i, imode, mode in zip(count(), imodes, modes):
            vel = self.results[imode, :, self.ivelocity].ravel()
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()

            #iplot = np.where(freq > 0.0)
            #damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel, freq, symbols[i])

            iplot = np.where(freq != np.nan)
            damp_axes.plot(vel[iplot], damping[iplot], symbols[i], label='Mode %i' % mode)
            freq_axes.plot(vel[iplot], freq[iplot], symbols[i])

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
                    print(line.strip())
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
                    print()
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

            point_sline = f06_file.readline().split()
            iline += 1
            mode = int(point_sline[2])
            method = point_sline[-1]  # 13 for PN, 5 for PK

            #log.debug(point_sline)
            if method == 'PK':
                mach = float(point_sline[6])
                density_ratio = float(point_sline[10])
                #method = point_sline[13]
                if mode == 1:
                    print('# iline mode mach density_ratio method')
                print(iline, mode, mach, density_ratio, method)
            elif method == 'PKNL':
                mach = None
                density_ratio = None
                if mode == 1:
                    print('# iline mode method')
                print(iline, mode, method)
                f06_file.readline()
                iline += 1
            else:
                raise NotImplementedError(point_sline)

            if mode in modes:
                log.warning('found existing mode...')
                continue
            modes.append(mode)

            # blanks
            f06_file.readline()
            f06_file.readline()
            iline += 2

            lines = []

            # KFREQ  1./KFREQ                      VELOCITY  DAMPING  FREQUENCY   COMPLEX   EIGENVALUE - PK
            # KFREQ  1./KFREQ  DENSITY   MACH NO.  VELOCITY  DAMPING  FREQUENCY   COMPLEX   EIGENVALUE - PKNL
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

if __name__ == '__main__':
    plot_flutter_f06('bah_plane.f06')
