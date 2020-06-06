from itertools import count
from typing import  List, Dict, Union, Tuple, Iterable, Optional, Any

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    IS_MATPLOTLIB = True
except ImportError:  # pragma: no cover
    IS_MATPLOTLIB = False

from pyNastran.utils.atmosphere import (
    get_alt_for_density, atm_density,
    convert_altitude, convert_velocity, convert_density, convert_pressure,
)
from pyNastran.utils import object_attributes, object_methods


class FlutterResponse:
    """storage object for single subcase SOL 145 results"""

    def __init__(self, subcase: int, configuration: str,
                 xysym: str, xzsym: str,
                 mach: float, density_ratio: float, method: str,
                 modes: List[int], results: Any,
                 f06_units: Union[None, str, Dict[str, str]]=None,
                 out_units: Union[None, str, Dict[str, str]]=None,
                 make_alt: bool=False) -> None:
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
                f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3'}
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
        self.make_alt = make_alt
        self.f06_units = f06_units
        self.out_units = out_units
        required_keys = ['altitude', 'velocity', 'eas', 'density', 'dynamic_pressure']
        for key in required_keys:
            assert key in f06_units, 'key=%r not in f06_units=%s' % (key, f06_units)
            assert key in out_units, 'key=%r not in out_units=%s' % (key, out_units)
        for key in f06_units:
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

        results = _asarray(results)
        if self.method in ['PK', 'KE']:
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
        self._symbols = []  # type: List[str]
        self._colors = []  # type: List[str]
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
        if self.make_alt:
            rho_in_slug_ft3 = rho * kdensityi
            alt_ft = [get_alt_for_density(densityi, density_units='slug/ft^3',
                                          alt_units='ft', nmax=20)
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
        is defined.

        """
        # rgbkm - max of 35 combinations
        # C0-10 - max of 70 combinations
        if colors is None:
            #colors = ['r', 'g', 'b', 'k', 'm'] # 5
            colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'] # 10
        if symbols is None:
            symbols = ['o', '*', 'x', 'v', '>', '<', '^'] # 7
        for symbol in symbols:
            for color in colors:
                self._symbols.append(symbol)
                self._colors.append(color)

    def set_plot_options(self, noline: bool=False) -> None:
        self.noline = noline

    def _get_unit_factor(self, name: str) -> Tuple[float, str]:
        if not self.f06_units or not self.out_units:
            msg = 'name=%r f06_units=%s out_units=%s' % (name, self.f06_units, self.out_units)
            raise RuntimeError(msg)
        unit_f06 = self.f06_units[name]
        unit_out = self.out_units[name]

        #print('name=%s unit_f06=%r unit_out=%r' % (name, unit_f06, unit_out))
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

    def plot_vg(self, fig=None, modes=None,
                plot_type='tas',
                xlim=None, ylim_damping=None,
                clear=False, legend=True,
                png_filename=None, show=True, **kwargs):
        """
        Make a V-g plot

        See ``plot_root_locus`` for arguments
        """
        ix, xlabel = self._plot_type_to_ix_xlabel(plot_type)
        ylabel = 'Damping'
        iy = self.idamping
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim_damping,
                       show=show, clear=clear, legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    #@property
    #def flutter_speed(self, modes=None):
        #"""gets the flutter speed"""
        #if modes is None:
            #modes = self.modes
        #else:
            #modes = np.asarray(modes)

    def plot_root_locus(self, modes=None,
                        fig=None, axes=None,
                        xlim=None, ylim=None,
                        show=True, clear=False, close=False, legend=True,
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
        axes : plt.axes
            axes object
        xlim : list[float/None, float/None]
            the x plot limits
        ylim : list[float/None, float/None]
            the y plot limits
        show : bool; default=True
            show the plot
        clear : bool; default=False
            clear the plot
        legend : bool; default=True
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
            1.0 - no transparency / opaque
            0.0 - fully transparent

        """
        xlabel = 'Eigenvalue (Real)'
        ylabel = 'Eigenvalue (Imaginary)'
        ix = self.ieigr
        iy = self.ieigi
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, axes=axes, xlim=xlim, ylim=ylim,
                       show=show, clear=clear, close=close, legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    def _plot_x_y(self, ix, iy, xlabel, ylabel, scatter, modes=None,
                  fig=None, axes=None,
                  xlim=None, ylim=None,
                  show=True, clear=False, close=False, legend=True,
                  png_filename=None,
                  **kwargs):
        """builds the plot"""
        self.fix()
        if kwargs is None:
            kwargs = {}

        modes, imodes = _get_modes_imodes(self.modes, modes)

        if fig is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)
        symbols = self._symbols
        colors = self._colors
        linestyle = 'None' if self.noline else '-'

        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            color = colors[i]
            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            ys = self.results[imode, :, iy].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            axes.plot(xs[iplot], ys[iplot], color=color, marker=symbol, label='Mode %i' % mode,
                      linestyle=linestyle, markersize=0)

            if scatter:
                scatteri = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                #axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=symbol[0], marker=symbol[1])
                axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=color, marker=symbol)

        axes.grid(True)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        if xlim:
            axes.set_xlim(xlim)
        if ylim:
            axes.set_ylim(ylim)

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
            fig.clear()
        if close:
            plt.close()
        return axes

    def _plot_x_y2(self, ix, iy1, iy2, xlabel, ylabel1, ylabel2, scatter, modes=None,
                   fig=None, axes1=None, axes2=None,
                   xlim=None, ylim1=None, ylim2=None,
                   nopoints=False, noline=False,
                   show=True, clear=False, close=False, legend=True,
                   png_filename=None,
                   **kwargs):
        """
        Builds the plot

        Parameters
        ----------
        scatter : bool
            draw the points with growth along the x-axis

        """
        self.fix()
        if kwargs is None:
            kwargs = {}

        modes, imodes = _get_modes_imodes(self.modes, modes)

        if fig is None:
            fig = plt.figure()
            gridspeci = gridspec.GridSpec(2, 4)
            axes1 = fig.add_subplot(gridspeci[0, :3])
            axes2 = fig.add_subplot(gridspeci[1, :3], sharex=axes1)

        if xlim:
            axes1.set_xlim(xlim)
        if ylim1:
            axes1.set_ylim(ylim1)
        if ylim2:
            axes2.set_ylim(ylim2)

        symbols = self._symbols
        colors = self._colors

        linestyle = 'None' if noline else '-'
        if nopoints: # and noline is False:
            scatter = False

        markersize = None
        if noline:
            markersize = 0

        #showline = not noline
        #showpoints = not nopoints

        legend_elements = []
        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            color = colors[i]

            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            y1s = self.results[imode, :, iy1].ravel()
            y2s = self.results[imode, :, iy2].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)

            # plot the line
            label = 'Mode %i' % mode
            legend_element = Line2D([0], [0], color=color, marker=symbol, label=label, linestyle=linestyle)
            if nopoints:
                symbol = 'None'
            if scatter:
                scatteri = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                axes1.scatter(xs[iplot], y1s[iplot],
                              s=scatteri, color=color, marker=symbol)
                axes2.scatter(xs[iplot], y2s[iplot],
                              s=scatteri, color=color, marker=symbol)

                # Draw the line
                axes1.plot(xs[iplot], y1s[iplot], marker=symbol, label=label,
                           color=color, markersize=markersize, linestyle=linestyle)
                axes2.plot(xs[iplot], y2s[iplot], marker=symbol,
                           color=color, markersize=markersize, linestyle=linestyle)
            else:
                axes1.plot(xs[iplot], y1s[iplot], marker=symbol, label=label,
                           color=color, markersize=markersize, linestyle=linestyle)
                axes2.plot(xs[iplot], y2s[iplot], marker=symbol,
                           color=color, markersize=markersize, linestyle=linestyle)
            legend_elements.append(legend_element)

        axes1.grid(True)
        axes1.set_xlabel(xlabel)
        axes1.set_ylabel(ylabel1)

        axes2.grid(True)
        axes2.set_xlabel(xlabel)
        axes2.set_ylabel(ylabel2)

        title = 'Subcase %i' % self.subcase
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        if legend:
            axes1.legend(handles=legend_elements, **kwargs)

        if show:
            plt.show()
        if png_filename:
            plt.savefig(png_filename)
        if clear:
            fig.clear()
        if close:
            plt.close()

    def plot_kfreq_damping(self, modes=None,
                           plot_type='tas',
                           fig=None, damp_axes=None, freq_axes=None,
                           xlim=None,
                           show=True, clear=False, close=False, legend=True,
                           png_filename=None,
                           ylim_damping=None,
                           ylim_kfreq=None,
                           nopoints=False, noline=False,
                           **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments
        """
        ylabel1 = 'Damping'
        ylabel2 = 'KFreq [rad]'

        ix, xlabel = self._plot_type_to_ix_xlabel(plot_type)
        iy1 = self.idamping
        iy2 = self.ikfreq
        scatter = True
        self._plot_x_y2(ix, iy1, iy2, xlabel, ylabel1, ylabel2, scatter,
                        modes=modes, fig=fig, axes1=damp_axes, axes2=freq_axes,
                        xlim=xlim, ylim1=ylim_damping, ylim2=ylim_kfreq,
                        nopoints=nopoints, noline=noline,
                        show=show, clear=clear, close=close,
                        legend=legend,
                        png_filename=png_filename,
                        **kwargs)

    def plot_kfreq_damping2(self, modes=None,
                            fig=None,
                            xlim=None, ylim=None,
                            show=True, clear=False, close=False, legend=True,
                            png_filename=None,
                            **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments

        """
        xlabel = 'KFreq [rad]'
        ylabel1 = 'Damping'
        ylabel2 = 'Frequency'
        ix = self.ikfreq
        iy1 = self.idamping
        iy2 = self.ifreq
        scatter = True
        self._plot_x_y2(ix, iy1, iy2, xlabel, ylabel1, ylabel2, scatter,
                        modes=modes, fig=fig, xlim=xlim, ylim=ylim,
                        show=show, clear=clear, close=close,
                        legend=legend,
                        png_filename=png_filename,
                        **kwargs)

    def fix(self):
        """attempts to fix the mode switching"""
        #print(self.names)

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

    def plot_vg_vf(self, fig=None, damp_axes=None, freq_axes=None, modes=None,
                   plot_type='tas',
                   clear=False, close=False, legend=True,
                   xlim=None, ylim_damping=None, ylim_freq=None,
                   nopoints=False, noline=False, png_filename=None, show=None):
        """
        Make a V-g and V-f plot

        Parameters
        ----------
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        legend : bool; default=True
            should the legend be shown
        """
        #self.fix()
        if fig is None:
            fig = plt.figure() # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        modes, imodes = _get_modes_imodes(self.modes, modes)
        symbols = self._symbols
        colors = self._colors

        if nopoints:
            symbols = ['None'] * len(symbols)
        linestyle = 'None' if noline else '-'

        ix, xlabel = self._plot_type_to_ix_xlabel(plot_type)
        for i, imode, mode in zip(count(), imodes, modes):
            color = colors[i]
            symbol = symbols[i]

            vel = self.results[imode, :, ix].ravel()
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()

            #iplot = np.where(freq > 0.0)
            #damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel, freq, symbols[i])

            #iplot = np.where(freq != np.nan)
            #damp_axes.plot(vel[iplot], damping[iplot], symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel[iplot], freq[iplot], symbols[i])
            #print(color, symbol, linestyle)
            damp_axes.plot(vel, damping, color=color, marker=symbol, linestyle=linestyle, label='Mode %i' % mode)
            freq_axes.plot(vel, freq, color=color, marker=symbol, linestyle=linestyle)

        damp_axes.set_xlabel(xlabel)
        freq_axes.set_xlabel(xlabel)
        damp_axes.set_ylabel('Damping')
        freq_axes.set_ybound(lower=0.)

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
            damp_axes.legend(fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2)
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
            fig.clear()
        if close:
            plt.close()

    def export_to_veas(self, veas_filename: str, modes: Optional[List[int]]=None) -> None:
        """
        Exports a ZONA .veas file

        Parameters
        ----------
        veas_filename : str
            the filename to write
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N

        *.VEAS

        ' DAMPING & FREQUENCY X-Y PLOT FILE OF PLTVG SETID=     714 FOR FLUTTER/ASE ID=     714 NMODE=   12'
        'EQUIVALENT V   G,MODE--1   G,MODE--2   G,MODE--3   G,MODE--4   G,MODE--5   G,MODE--6   G,MODE--7   G,MODE--8   G,MODE--9   G,MODE-10   G,MODE-11   G,MODE-12EQUIVALENT V WHZ,MODE--1 WHZ,MODE--2 WHZ,MODE--3 WHZ,MODE--4 WHZ,MODE--5 WHZ,MODE--6 WHZ,MODE--7 WHZ,MODE--8 WHZ,MODE--9 WHZ,MODE-10 WHZ,MODE-11 WHZ,MODE-12'
        '  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.4277E+01  7.7705E+01  8.8016E+01  2.2346E+02  2.4656E+02  3.2366E+02  4.5518E+02  4.9613E+02  6.7245E+02  7.3197E+02  8.0646E+02  0.0000E+00'
        '  4.9374E-01 -1.6398E-03 -5.4768E-04 -2.3136E-04 -3.5776E-04 -6.2358E-05 -3.0348E-04 -7.7492E-05 -2.0301E-05 -1.7733E-04 -9.2383E-05 -1.2849E-05 -2.8854E-05  4.9374E-01  1.4272E+01  7.7726E+01  8.8010E+01  2.2347E+02  2.4655E+02  3.2364E+02  4.5516E+02  4.9612E+02  6.7245E+02  7.3195E+02  8.0646E+02  0.0000E+00'
        """
        modes, nmodes = self._modes_nmodes(modes)

        damping_modes = []
        omega_modes = []
        for mode in modes:
            if mode < 10:
                gmode = '   G,MODE--%i' % mode
                wmode = ' WHZ,MODE--%i' % mode
            elif mode < 100:
                gmode = '   G,MODE-%2i' % mode
                wmode = ' WHZ,MODE-%2i' % mode
            else:
                gmode = '   G,MODE%3i' % mode
                wmode = ' WHZ,MODE%3i' % mode
            damping_modes.append(gmode)
            omega_modes.append(wmode)

        headers = ['EQUIVALENT V'] + damping_modes + ['EQUIVALENT V'] + omega_modes
        with open(veas_filename, 'w') as veas_file:
            veas_file.write(' DAMPING & FREQUENCY X-Y PLOT FILE OF PLTVG '
                            'SETID=       1 FOR FLUTTER/ASE ID=       1 NMODE=  %3i\n' % nmodes)
            veas_file.write(''.join(headers) + '\n')
            nspeeds = self.results.shape[1]
            for i in range(nspeeds):
                damping = self.results[:, i, self.idamping]
                eas = self.results[0, i, self.ieas]
                omega = self.results[:, i, self.ifreq] # in Hz (what WHZ means)

                values = [eas] + damping.tolist() + [eas] + omega.tolist()
                str_values = (' %11.4E' % value for value in values)
                veas_file.write(''.join(str_values) + '\n')

    def _imodes(self, modes):
        """gets the imodes from the modes"""
        if modes is None:
            nmodes = self.results.shape[0]
            imodes = range(nmodes)
        else:
            imodes = modes - 1
        return imodes

    def _modes_nmodes(self, modes: Optional[Iterable[int]]) -> Tuple[Iterable[int], int]:
        """gets the modes and nmodes"""
        if modes is None:
            nmodes = self.results.shape[0]
            modes = range(1, nmodes + 1)
        else:
            nmodes = max(modes)
        return modes, nmodes

    def export_to_f06(self, f06_filename: str,
                      modes: Optional[List[int]]=None,
                      page_stamp: Optional[str]=None,
                      page_num: int=1) -> int:
        if page_stamp is None:
            page_stamp = 'PAGE %i'
        # nmodes, vel, res
        imodes = self._imodes(modes)

        with open(f06_filename, 'w') as f06_file:
            for imode in imodes:
                #'      MACH 0.0                                                                                                                      '
                f06_file.write('0                                                                                                            SUBCASE %i\n' % self.subcase)
                f06_file.write('0                                                       FLUTTER  SUMMARY\n')
                f06_file.write('                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = ASYMMETRIC\n')
                f06_file.write('       POINT = %4i     METHOD = %s\n' % (imode + 1, self.method))
                f06_file.write('\n')
                f06_file.write('\n')

                f06_file.write('    KFREQ          1./KFREQ       DENSITY     MACH NO.      VELOCITY       DAMPING     FREQUENCY      COMPLEX   EIGENVALUE\n')
                for res in self.results[imode, :, :9]:
                    #print(res)
                    kfreq, kfreqi, rho, mach, vel, damp, freq, eigr, eigi = res
                    #                 kfreq      ikfreq  rho      mach     vel    damp   freq    eigr       eigi
                    f06_file.write(' %8.4f      %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n' % (
                        kfreq, kfreqi, rho, mach, vel, damp, freq, eigr, eigi,
                    ))
                #'1                                                                          DECEMBER  14, 2018  MSC.NASTRAN  6/17/05   PAGE    12\n'
                f06_file.write(page_stamp % page_num)
                page_num += 1
        return page_num

    def export_to_zona(self, zona_filename: str,
                       modes: Optional[List[int]]=None,
                       xlim: Optional[List[float]]=None,
                       plot_type: str='tas',
                       damping_ratios: Optional[List[float]]=None) -> str:
        """
        Writes a custom ZONA flutter file

        Parameters
        ----------
        zona_filename : str
            the filename to write
        modes : List[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N

        TODO: not done
        """
        if damping_ratios is None:
            damping_ratios = [0., 0.01, 0.02, 0.05, 0.1, 0.15]
        if xlim is None:
            xlim = [None, None]

        modes, imodes = _get_modes_imodes(self.modes, modes)
        #unused_legend_items = ['Mode %i' % mode for mode in modes]
        ix, unused_xlabel = self._plot_type_to_ix_xlabel(plot_type)

        # these are the required damping levels to plot
        msg = ''
        for damping_ratio in damping_ratios:
            msgi = ''
            xlim_min = xlim[0]
            xlim_max = xlim[1]
            for unused_i, imode, mode in zip(count(), imodes, modes):
                vel = self.results[imode, :, ix].ravel()
                damping = self.results[imode, :, self.idamping].ravel()
                freq = self.results[imode, :, self.ifreq].ravel()

                # consider flutter to be at 0, -0.01, ... damping ratio
                # (so we have a damping margin)
                inot_nan = np.isfinite(damping)
                idamp = np.where(damping[inot_nan] > -damping_ratio)[0]

                if len(idamp) == 0:
                    continue

                # limit the plot based on the xlimits
                veli = vel[inot_nan][idamp]
                if xlim_min is None and xlim_max is None:
                    jvel = None
                elif xlim_min is not None and xlim_max is not None:
                    jvel = np.where((xlim_min <= veli) & (veli <= xlim_max))[0]
                elif xlim_min is not None:
                    jvel = np.where(xlim_min <= veli)[0]
                elif xlim_max is not None:
                    jvel = np.where(veli <= xlim_max)[0]
                else:  # pragma: no cover
                    raise RuntimeError('xlim min/max are unclear')

                if jvel is None or len(jvel) == 0:
                    continue

                if jvel is None:
                    # no xlimits
                    idampi = idamp[0]
                    veli = vel[inot_nan][idampi]
                    dampi = damping[inot_nan][idampi]
                    freqi = freq[inot_nan][idampi]
                else:
                    # xlimits
                    jveli = jvel[0]
                    veli = vel[inot_nan][idamp][jveli]
                    dampi = damping[inot_nan][idamp][jveli]
                    freqi = freq[inot_nan][idamp][jveli]
                msgi += '%s %s %s %s\n' % (mode, veli, dampi, freqi)
            if msgi:
                msg += 'mode, V, damp, freq: (damping ratio=%s)\n%s\n' % (damping_ratio, msgi)

        ## TODO: doesn't always have data...
        with open(zona_filename, 'w') as zona_file:
            zona_file.write(msg)
        return msg

    def _plot_type_to_ix_xlabel(self, plot_type):
        """helper method for ``plot_vg_vf``"""
        plot_type = plot_type.lower()
        if plot_type == 'tas':
            ix = self.ivelocity
            velocity_units = self.out_units['velocity']
            xlabel = 'Velocity [%s]' % velocity_units
        elif plot_type == 'eas':
            ix = self.ieas
            velocity_units = self.out_units['eas']
            xlabel = 'Equivalent Airspeed [%s]' % velocity_units
        elif plot_type == 'alt':
            ix = self.ialt
            alt_units = self.out_units['altitude']
            xlabel = 'Altitude [%s]' % alt_units
        elif plot_type == 'kfreq':
            ix = self.ikfreq
            xlabel = 'Reduced Frequency [rad]'
        elif plot_type == 'rho':
            ix = self.idensity
            density_units = self.out_units['density']
            xlabel = 'Density [%s]' % density_units
        elif plot_type == 'freq':
            ix = self.ifreq
            xlabel = 'Frequency [Hz]'
        elif plot_type in ['1/kfreq', 'ikfreq']:
            ix = self.ikfreq_inv
            xlabel = '1/KFreq [1/rad]'
        elif plot_type == 'eigr':
            ix = self.ieigr
            xlabel = 'Eigenvalue (Real)'
        elif plot_type == 'eigr':
            ix = self.ieigi
            xlabel = 'Eigenvalue (Imaginary)'
        elif plot_type in ['damp', 'damping']:
            ix = self.idamping
            xlabel = 'Damping'
        else:
            raise NotImplementedError("plot_type=%r not in ['tas', 'eas', 'alt', 'kfreq', "
                                      "'1/kfreq', 'freq', 'damp', 'eigr', 'eigi']")
        return ix, xlabel

    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
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
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip,
                                 filter_properties=filter_properties)

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


def _get_modes_imodes(all_modes, modes):
    """gets the index of the modes to plot"""
    if modes is None:
        modes = all_modes
    elif isinstance(modes, slice):
        start = modes.start
        if modes.stop is None:
            stop = len(all_modes) + 1
        stop = modes.stop
        step = modes.step
        if step is None:
            step = 1
        modes = np.unique(range(start, stop, step))
    elif len(modes) == 0:
        raise RuntimeError('modes = %s' % modes)
    else:
        assert 0 not in modes, modes
        modes = np.unique(modes)
    assert 0 not in modes, modes

    if modes.max() > all_modes.max():
        imodes = np.where(modes <= all_modes.max())
        modes = modes[imodes]
    if len(modes) == 0:
        raise RuntimeError('No modes to plot...')
    imodes = np.searchsorted(all_modes, modes)
    return modes, imodes

def _asarray(results):
    """casts the results array"""
    try:
        results = np.asarray(results, dtype='float64')
    except:
        results2 = []
        fix_kfreq = False
        for mode_result in results:
            mode_result2 = []
            for row in mode_result:
                #print(row)
                row2 = []
                for i, row_entry in enumerate(row):
                    if i == 0 and '****' in row_entry:
                        row_entry2 = np.nan
                        fix_kfreq = True
                    elif i == 1 and 'INF' in row_entry:
                        row_entry2 = np.inf
                    elif i in [2, 3, 4] and row_entry in ['UNSTABL-SYSTEM', 'STABLE-SYSTEM-']:
                        row_entry2 = np.nan
                    else:
                        try:
                            row_entry2 = float(row_entry)
                        except:
                            raise ValueError(f'i={i} row_entry={row_entry!r}')
                    row2.append(row_entry2)
                mode_result2.append(row2)
            results2.append(mode_result2)
        results = np.array(results2, dtype='float64')

        if fix_kfreq:
            #inan = np.isnan(results[:, :, :])
            #kfreq_inv = results[:, :, 1]
            #print('kfreq.shape', kfreq_inv.shape)
            #print('inan.shape', inan.shape)
            #results[inan] = 1 / kfreq_inv[inan]
            nmodes = results.shape[0]
            for imode in range(nmodes):
                inan = np.isnan(results[imode, :, 0])
                if len(inan) == 0:
                    continue
                #print(inan)
                results[imode, inan, 0] = 1. / results[imode, inan, 1]
    return results
