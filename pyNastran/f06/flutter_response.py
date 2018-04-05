from __future__ import print_function
from itertools import count
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pyNastran.utils.atmosphere import get_alt_for_density
from pyNastran.utils.atmosphere2 import (
    convert_altitude, convert_density, convert_pressure, convert_velocity,
    atm_density,
)
from pyNastran.utils import object_attributes, object_methods


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
        required_keys = ['altitude', 'velocity', 'eas', 'density', 'dynamic_pressure']
        for key in required_keys:
            assert key in f06_units, 'key=%r not in f06_units=%s' % (key, f06_units)
            assert key in out_units, 'key=%r not in out_units=%s' % (key, out_units)
        for key in f06_units.keys():
            assert key in required_keys, 'key=%r not in required_keys=%s' % (key, required_keys)


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

            # (imode, istep, iresult)
            results[:, :, self.ivelocity] *= kvel
            # velocity is the target
            self.names = ['kfreq', '1/kfreq', 'velocity', 'damping', 'freq', 'eigr', 'eigi']

        elif self.method == 'PKNL':
            # velocity is the target
            self.names = ['kfreq', '1/kfreq', 'density', 'velocity', 'damping',
                          'freq', 'eigr', 'eigi', 'eas', 'q', 'alt']

            #KFREQ  1./KFREQ  DENSITY  MACH  VELOCITY  DAMPING  FREQUENCY  COMPLEX EIGENVALUE
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
            self.set_pknl_results(results)
        else:
            raise NotImplementedError(method)

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

    def set_pknl_results(self, results):
        density_units_in = self.f06_units['density']

        # in/s
        vel = results[:, :, self.ivelocity]#.ravel()

        # slinch/in^3 - in_units
        rho = results[:, :, self.idensity]#.ravel()

        # good
        rho_ref = atm_density(0., R=1716., alt_units='ft',
                              density_units=density_units_in)

        q = 0.5 * rho * vel**2
        #eas  = (2 * q / rho_ref)**0.5

        # eas = V * sqrt(rho / rhoSL)
        keas = self._get_unit_factor('eas')[0]
        eas = vel * np.sqrt(rho / rho_ref) * keas
        #density_units2 = self.out_units['density']

        altitude_units = self.out_units['altitude']

        #print('density_units_in=%r density_units2=%r' % (density_units_in, density_units2))
        kdensityi = convert_density(1., density_units_in, 'slug/ft^3')
        kvel = self._get_unit_factor('velocity')[0]
        kdensity = self._get_unit_factor('density')[0]
        kpressure = kdensityi * kvel ** 2

        vel *= kvel
        make_alt = False
        if make_alt:
            rho_in_slug_ft3 = rho * kdensityi
            alt_ft = [get_alt_for_density(densityi, nmax=20)
                      for densityi in rho_in_slug_ft3.ravel()]

            ft_to_alt_unit = convert_altitude(1., 'ft', altitude_units)
            alt = np.array(alt_ft, dtype='float64').reshape(vel.shape) * ft_to_alt_unit

            rho *= kdensity
            results2 = np.dstack([results, eas, q * kpressure, alt])
        else:
            #kpressure = 1.
            rho *= kdensity
            results2 = np.dstack([results, eas, q * kpressure])

        results2[:, :, self.idensity] = rho
        results2[:, :, self.ivelocity] = vel
        self.results = results2

    def generate_symbols(self, colors=None, symbols=None):
        """
        This symbol list is taken from a series of "good" colors (e.g. not yellow)
        and easily distinguishable shapes.  Far more combinations that is necessary
        is defined
        """
        # max of 35 combinations
        if colors is None:
            colors = ['r', 'g', 'b', 'k', 'm'] # 5
        if symbols is None:
            symbols = ['o', '*', 'x', 'v', '>', '<', '^'] # 7
        self._symbols = []
        for symbol in symbols:
            for color in colors:
                self._symbols.append(color + symbol)

    def set_plot_options(self, noline=False):
        self.noline = noline

    def _get_unit_factor(self, name):
        if not self.f06_units or not self.out_units:
            msg = 'name=%r f06_units=%s out_units=%s' % (name, self.f06_units, self.out_units)
            raise RuntimeError(msg)
        unit_f06 = self.f06_units[name]
        unit_out = self.out_units[name]

        print('name=%s unit_f06=%r unit_out=%r' % (name, unit_f06, unit_out))
        if name in ['velocity', 'eas']:
            factor = convert_velocity(1., unit_f06, unit_out)
        elif name == 'altitude':
            factor = convert_altitude(1., unit_f06, unit_out)
        elif name == 'density':
            factor = convert_density(1., unit_f06, unit_out)
        elif name in ['pressure', 'dynamic_pressure']:
            factor = convert_pressure(1., unit_f06, unit_out)
        else:
            raise NotImplementedError(name)

        if self.out_units is not None:
            units = self.out_units[name]
        else:
            units = 'units'
        return factor, units

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
            xs = self.results[imode, :, ix].ravel()
            ys = self.results[imode, :, iy].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            axes.plot(xs[iplot], ys[iplot], symbol, label='Mode %i' % mode, markersize=0)

            if scatter:
                s = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                axes.scatter(xs[iplot], ys[iplot], s=s, color=symbol[0], marker=symbol[1])

        axes.grid(True)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)

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

    def _plot_x_y2(self, ix, iy, iy2, xlabel, ylabel1, ylabel2, scatter, modes=None,
                   fig=None,
                   xlim=None, ylim1=None, ylim2=None,
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
        gridspeci = gridspec.GridSpec(2, 4)

        damp_axes = fig.add_subplot(gridspeci[0, :3])
        freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        if xlim:
            damp_axes.set_xlim(xlim)
        if ylim1:
            damp_axes.set_ylim(ylim1)
        if ylim2:
            damp_axes.set_ylim(ylim2)
        symbols = self.symbols

        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            y1s = self.results[imode, :, iy].ravel()
            y2s = self.results[imode, :, iy2].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            damp_axes.plot(xs[iplot], y1s[iplot], symbol, label='Mode %i' % mode, markersize=0)
            freq_axes.plot(xs[iplot], y2s[iplot], symbol, label='Mode %i' % mode, markersize=0)

            if scatter:
                s = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                damp_axes.scatter(xs[iplot], y1s[iplot], s=s, color=symbol[0], marker=symbol[1])
                freq_axes.scatter(xs[iplot], y2s[iplot], s=s, color=symbol[0], marker=symbol[1])

        damp_axes.grid(True)
        damp_axes.set_xlabel(xlabel)
        damp_axes.set_ylabel(ylabel1)

        freq_axes.grid(True)
        freq_axes.set_xlabel(xlabel)
        freq_axes.set_ylabel(ylabel2)

        velocity_units = self.out_units['eas']
        _set_eas_tas_xlabel('eas', damp_axes, freq_axes, velocity_units)

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        if legend:
            damp_axes.legend(**kwargs)

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            plt.clear()


    def plot_kfreq_damping(self, modes=None,
                           fig=None,
                           xlim=None,
                           show=True, clear=False, legend=True,
                           png_filename=None,
                           ylim_damping=None,
                           ylim_kfreq=None,
                           **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments
        """
        xlabel = 'Equivalent Airspeed'
        ylabel = 'Damping'
        ylabel2 = 'KFreq'
        ix = self.ieas
        iy = self.idamping
        iy2 = self.ikfreq
        scatter = True
        self._plot_x_y2(ix, iy, iy2, xlabel, ylabel, ylabel2, scatter,
                        modes=modes, fig=fig,
                        xlim=xlim, ylim1=ylim_damping, ylim2=ylim_kfreq,
                        show=show,
                        clear=clear,
                        legend=legend,
                        png_filename=png_filename,
                        **kwargs)

    def plot_kfreq_damping2(self, modes=None,
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
        ylabel2 = 'Frequency'
        ix = self.ikfreq
        iy = self.idamping
        iy2 = self.ifreq
        scatter = True
        self._plot_x_y2(ix, iy, iy2, xlabel, ylabel, ylabel2, scatter,
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

    def plot_vg_vf(self, fig=None, damp_axes=None, freq_axes=None, modes=None, show=None,
                   plot_type='tas',
                   png_filename=None, clear=False, legend=None,
                   xlim=None, ylim_damping=None, ylim_freq=None, nopoints=False):
        """
        Make a V-g and V-f plot

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        """
        #self.fix()
        if fig is None:
            fig = plt.figure() # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        modes, imodes = self._get_modes_imodes(modes)
        symbols = self.symbols

        legend_items = ['Mode %i' % mode for mode in modes]
        for i, imode, mode in zip(count(), imodes, modes):
            if plot_type == 'tas':
                vel = self.results[imode, :, self.ivelocity].ravel()
                velocity_units = self.out_units['velocity']
            elif plot_type == 'eas':
                vel = self.results[imode, :, self.ieas].ravel()
                velocity_units = self.out_units['eas']
            else:
                raise NotImplementedError("plot_type=%r not in ['tas', 'eas']")
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()

            #iplot = np.where(freq > 0.0)
            #damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel, freq, symbols[i])

            #iplot = np.where(freq != np.nan)
            #damp_axes.plot(vel[iplot], damping[iplot], symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel[iplot], freq[iplot], symbols[i])
            if symbols and not nopoints:
                symbol = symbols[i]
                damp_axes.plot(vel, damping, symbol, label='Mode %i' % mode)
                freq_axes.plot(vel, freq, symbol)
            else:
                damp_axes.plot(vel, damping, label='Mode %i' % mode)
                freq_axes.plot(vel, freq)

        _set_eas_tas_xlabel(plot_type, damp_axes, freq_axes, velocity_units)

        damp_axes.set_ylabel('Damping')
        damp_axes.grid(True)
        if xlim is not None:
            damp_axes.set_xlim(xlim)
        if ylim_damping is not None:
            damp_axes.set_ylim(ylim_damping)

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
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip)

def _set_eas_tas_xlabel(plot_type, damp_axes, freq_axes, velocity_units):
    if plot_type == 'tas':
        damp_axes.set_xlabel('Velocity [%s]' % velocity_units)
        freq_axes.set_xlabel('Velocity [%s]' % velocity_units)
    elif plot_type == 'eas':
        damp_axes.set_xlabel('Equivalent Airspeed [%s]' % velocity_units)
        freq_axes.set_xlabel('Equivalent Airspeed [%s]' % velocity_units)
    else:
        raise NotImplementedError("plot_type=%r not in ['tas', 'eas']")
