from __future__ import annotations
from itertools import count
from typing import Iterable, Optional, Any, TYPE_CHECKING

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

from pyNastran.utils.atmosphere import (
    get_alt_for_density, atm_density,
    convert_altitude, convert_velocity, convert_density, convert_pressure,
)
from pyNastran.utils import object_attributes, object_methods
from pyNastran.utils.numpy_utils import float_types

if TYPE_CHECKING and IS_MATPLOTLIB:  # pragma: no cover
    from matplotlib.axes import Axes


class FlutterResponse:
    """storage object for single subcase SOL 145 results"""

    def __repr__(self) -> str:
        #from pyNastran.utils import object_stats
        #print(object_stats(self))
        #configuration : 'AEROSG2D'
        #f06_units : {'velocity': 'm/s', 'density': 'kg/m^3', 'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'm/s'}
        #ialt   : 11
        #idamping : 5
        #idensity : 2
        #ieas   : 9
        #ieigi  : 8
        #ieigr  : 7
        #ifreq  : 6
        #ikfreq : 0
        #ikfreq_inv : 1
        #imach  : 3
        #iq     : 10
        #ivelocity : 4
        #make_alt : False
        #method : 'PKNL'
        #modes  : array([1, 2])
        #names  : ['kfreq', '1/kfreq', 'density', 'velocity', 'damping', 'freq', 'eigr', 'eigi', 'eas', 'q', 'alt']
        #noline : False
        #out_units : {'velocity': 'in/s', 'density': 'slinch/in^3', 'altitude': 'ft', 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        xyz_sym = ''
        if hasattr(self, 'xysym'):
            xyz_sym += f'xysym  = {self.xysym!r}\n'
            xyz_sym += f'xzsym  = {self.xzsym!r}\n'

        msg = (
            'FlutterResponse:\n'
            f'subcase= {self.subcase:d}\n'
            f'{xyz_sym}'
            f'f06_units  = {self.f06_units}\n'
            f'out_units  = {self.out_units}\n'
            f'names  = {self.names}; n={len(self.names)}\n\n'
            f'method  = {self.method!r}\n'
            f'modes  = {self.modes}; n={len(self.modes)}\n'
            f'results.shape = {self.results.shape}; (nmodes, npoint, nresults)\n'
            #configuration : 'AEROSG2D'
            #density_ratio : 1.0
            #ikfreq : 0
            #ikfreq_inv : 1
            #ivelocity : 2
            #idamping : 3
            #ieigi  : 6
            #ieigr  : 5
            #ifreq  : 4
            #mach   : 0.0
            #make_alt : False
        )
        return msg

    @classmethod
    def from_nx(cls, method: str, fdata: np.ndarray,
                subcase_id: int=1, cref: float=1.0,
                is_xysym: bool=False, is_xzsym: bool=False):
        """
        Parameters
        ----------
        method : str
            expected PKNL
        subcase_id : int
            self explanatory
        cref : chord; default=1.0
            Reference length for reduced frequency
            used to compute b (reference semi-chord); found per the AERO card
        is_xysym : bool; default=False
            is the model symmetric about the xy plane?
        is_xzsym : bool; default=False
            is the model symmetric about the xz plane?

        """
        b = cref / 2.0
        configuration = 'AEROSG2D' #  TODO: what is this?
        xysym = '???'
        xzsym = '???'
        #[mach, rho, velocity, damping, freq]
        modes = None
        mach = 0.0
        density_ratio = 1.0
        assert method == 'PKNL', method
        nmodes, npoints, five = fdata.shape
        assert five == 5, fdata.shape

        #fdata   [                rho,     mach, velocity, damping, freq]
        #results [kfreq, 1/kfreq, density, mach, velocity, damping, freq, eigr, eigi]
        #
        # kfreq = omega*b / V
        results = np.zeros((nmodes, npoints, 9), dtype=fdata.dtype)
        #rhos = fdata[:, :, 0]
        #machs = fdata[:, :, 1]
        vels = fdata[:, :, 2]
        damps = fdata[:, :, 3]
        freqs = fdata[:, :, 4]
        omegas = 2 * np.pi * freqs
        kfreqs = omegas * b / vels
        #print(f'b = {b}')
        #print('omega:\n', omegas)
        #print('V:\n', vels)
        #print('kfreq:\n', kfreqs)

        # lambda = omega*(zeta + 1j)
        eigr = omegas * damps / 2
        eigc = omegas

        one_over_kfreq = np.full(kfreqs.shape, np.nan, dtype=kfreqs.dtype)
        ikfreq = np.where(kfreqs != 0)
        one_over_kfreq[ikfreq] = 1 / kfreqs[ikfreq]
        results[:, :, 0] = kfreqs
        results[:, :, 1] = one_over_kfreq
        results[:, :, 7] = eigr  # real part of eigenvalue
        results[:, :, 8] = eigc  # complex part of eigenvalue
        results[:, :, [2, 3, 4, 5, 6]] = fdata[:, :, [0, 1, 2, 3, 4]]
        modes = np.arange(nmodes, dtype='int32') + 1

        f06_units = get_flutter_units('english_in')
        #{
            #'altitude': 'ft',
            #'density': 'slinch/in^3',
            #'velocity': 'in/s',
            #'eas': 'in/s',
            #'dynamic_pressure': 'psi',
        #}
        #out_units = get_flutter_units('english_kt')
        out_units = get_flutter_units('english_in')
        #{
            #'altitude': 'ft',
            #'density': 'slug/ft^3',
            #'velocity': 'knots',
            #'eas': 'knots',
            #'dynamic_pressure': 'psf',
        #}
        resp = FlutterResponse(
            subcase_id, configuration, xysym, xzsym, mach, density_ratio,
            method, modes, results, f06_units=f06_units, out_units=out_units)
        if 0:  # pragma: no cover
            resp.plot_root_locus(modes=None, fig=None, axes=None, xlim=None, ylim=None,
                                 ncol=ncol, show=False,
                                 clear=False, close=False, legend=True, png_filename=None)
            resp.plot_vg_vf(fig=None, damp_axes=None, freq_axes=None, modes=None, plot_type='eas',
                             legend=True, xlim=None,
                            ylim_damping=None, ylim_freq=None, vd_limit=None,
                            damping_limit=None, nopoints=False, noline=False,
                            clear=False, close=False,
                            ncol=ncol, png_filename=None, show=False)
            resp.plot_kfreq_damping(modes=None, plot_type='tas', fig=None, damp_axes=None, freq_axes=None,
                                    xlim=None, ncol=ncol, legend=True,
                                    show=False, clear=False, close=False, png_filename=None,
                                    ylim_damping=None, ylim_kfreq=None, vd_limit=None, damping_limit=None,
                                    nopoints=False, noline=False)
            resp.plot_kfreq_damping2(modes=None, fig=None, xlim=None, ylim=None,
                                     show=True, clear=False, close=False,
                                     ncol=ncol, legend=True, png_filename=None)
        return resp

    def __init__(self, subcase: int, configuration: str,
                 xysym: str, xzsym: str,
                 mach: float, density_ratio: float, method: str,
                 modes: list[int], results: Any,
                 f06_units: None | str | dict[str, str]=None,
                 out_units: None | str | dict[str, str]=None,
                 make_alt: bool=False) -> None:
        """
        Parameters
        ----------
        subcase : int
            the subcase id
        method : str
            PK, PKNL, ???
        modes : list[int]; (default=None -> all)
            the modes; typically 1 to N
        results : varies
            method = PK
                list[list[float] * 7] * nmodes
                kfreq, 1/kfreq,                velocity, damping, freq, eigr, eigi
            method = PKNL
                list[list[float] * 9] * nmodes
                kfreq, 1/kfreq, density, mach, velocity, damping, freq, eigr, eigi

        Units
        -----
        Units are only applicable for quantities that have units
        (e.g. units on Mach don't do anything).
        All allowable fields are shown below.

        f06_units : dict[str] = str (default=None -> no units conversion)
            The units to read from the F06.
            PK method:
                f06_units = {'velocity' : 'in/s'}
                The velocity units are the units for the FLFACT card in the BDF
            PKNL method:
                f06_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3', 'altitude' : 'ft', 'dynamic_pressure': 'psi'}
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
                    'altitude' : 'ft',
                    'dynamic_pressure' : 'psf',
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
        if f06_units is None and out_units is None:
            pass
        else:
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
        else:  # pragma: no cover
            raise NotImplementedError(method)

        # c - cyan
        # b - black
        # r - red
        # g - green
        # m - magenta
        # y - yellow
        #colors = ['b', 'c', 'g', 'k', 'm', 'r'] #, 'y'
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
        self._symbols: list[str] = []
        self._colors: list[str] = []
        self.generate_symbols()

    def set_pknl_results(self, results: np.ndarray):
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
        kpressure = self._get_unit_factor('dynamic_pressure')[0]

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

    def generate_symbols(self, colors=None, symbols=None, imethod: int=0):
        """
        This symbol list is taken from a series of "good" colors (e.g. not yellow)
        and easily distinguishable shapes.  Far more combinations that is necessary
        is defined.

        Parameters
        ----------
        imethod: int; default=0
            imethod=0: loop over symbols and then colors
            imethod=1: loop over colors and then symbols

        """
        # rgbkm - max of 35 combinations
        # C0-10 - max of 70 combinations
        if colors is None:
            #colors = ['r', 'g', 'b', 'k', 'm'] # 5
            colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'] # 10

        if symbols is None:
            symbols = ['o', '*', 'x', 'v', '>', '<', '^'] # 7
        if imethod == 0:
            for symbol in symbols:
                for color in colors:
                    self._symbols.append(symbol)
                    self._colors.append(color)
        else:
            xmin = 0.0
            xmax = 0.8
            nmodes = len(self.modes)
            colormap_name = 'viridis'
            x = np.linspace(xmin, xmax, num=nmodes)[::-1]
            colormap = plt.get_cmap(colormap_name)
            colors = colormap(x)
            for color in colors:
                for symbol in symbols:
                    self._symbols.append(symbol)
                    self._colors.append(color)

    def set_plot_options(self, noline: bool=False) -> None:
        self.noline = noline

    def _get_unit_factor(self, name: str) -> tuple[float, str]:
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
        else:  # pragma: no cover
            raise NotImplementedError(name)

        if self.out_units is not None:
            units = self.out_units[name]
        else:
            units = 'units'
        return factor, units

    def plot_vg(self, fig=None, modes=None,
                plot_type='tas',
                xlim=None, ylim_damping=None,
                ncol: int=0,
                clear=False, legend=True,
                png_filename=None, show=True, **kwargs):
        """
        Make a V-g plot

        See ``plot_root_locus`` for arguments
        """
        ix, xlabel = self._plot_type_to_ix_xlabel(plot_type)
        ylabel = 'Viscous Damping'
        iy = self.idamping
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim_damping,
                       ncol=ncol,
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
                        ncol: int=0,
                        show: bool=True, clear: bool=False,
                        close: bool=False, legend: bool=True,
                        png_filename=None,
                        **kwargs):
        """
        Plots a root locus

        Parameters
        ----------
        modes : list[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        fig : plt.Figure
            the figure object
        axes : plt.Axes
            the axes object
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
        xlabel = r'Eigenvalue (Real); $\omega \zeta$'
        ylabel = r'Eigenvalue (Imaginary); $\omega$'
        ix = self.ieigr
        iy = self.ieigi
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, axes=axes, xlim=xlim, ylim=ylim,
                       ncol=ncol,
                       show=show, clear=clear, close=close, legend=legend,
                       png_filename=png_filename,
                       **kwargs)

    def _plot_x_y(self, ix: int, iy: int,
                  xlabel: str, ylabel: str,
                  scatter: bool,
                  modes=None,
                  fig=None, axes=None,
                  xlim=None, ylim=None,
                  ncol: int=0,
                  show: bool=True, clear: bool=False,
                  close: bool=False, legend: bool=True,
                  png_filename=None,
                  **kwargs):
        """builds the plot"""
        self.fix()
        if kwargs is None:
            kwargs = {}

        modes, imodes = _get_modes_imodes(self.modes, modes)
        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)

        if fig is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)

        symbols, colors = self._get_symbols_colors_from_modes(modes)
        linestyle = 'None' if self.noline else '-'

        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[i]
            color = colors[i]
            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            ys = self.results[imode, :, iy].ravel()

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            line = axes.plot(xs[iplot], ys[iplot],
                             color=color, marker=symbol, label=f'Mode {mode:d}',
                             linestyle=linestyle, markersize=0)

            if scatter:
                scatteri = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                #axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=symbol[0], marker=symbol[1])
                axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=color, marker=symbol)

        axes.grid(True)
        #axes.set_xlabel(xlabel  + '; _plot_x_y')
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
            # bbox_to_anchor=(1.125, 1.), ncol=ncol,
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
                   nopoints: bool=False, noline: bool=False,
                   ncol: int=0,
                   show: bool=True, clear: bool=False,
                   close: bool=False, legend: bool=True,
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
        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)

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

        symbols, colors = self._get_symbols_colors_from_modes(modes)

        linestyle = 'None' if noline else '-'
        if nopoints:  # and noline is False:
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
        #axes1.set_xlabel(xlabel + '; _plot_x_y2')
        axes1.set_ylabel(ylabel1)

        axes2.grid(True)
        axes2.set_xlabel(xlabel)
        axes2.set_ylabel(ylabel2)

        title = f'Subcase {self.subcase:d}'
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        if legend:
            legend_kwargs = {}
            for key, value in kwargs.items():
                if key in {'ylim'}:
                    continue
                legend_kwargs[key] = value

            axes1.legend(handles=legend_elements, **legend_kwargs)              # TODO: switch to figure...
            #axes1.legend(handles=legend_elements, ncol=ncol, **legend_kwargs)  # TODO: switch to figure...
            #fig.legend(handles=legend_elements, ncol=ncol, **legend_kwargs)

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
                           vd_limit=None, damping_limit=None,
                           nopoints=False, noline=False,
                           **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments
        """
        assert vd_limit is None or isinstance(vd_limit, float_types), vd_limit
        assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit

        ylabel1 = r'Viscous Damping; $g = 2 \gamma $'
        ylabel2 = r'KFreq [rad]; $ \omega c / (2 V)$'

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
                            show: bool=True, clear: bool=False,
                            close: bool=False, legend: bool=True,
                            png_filename=None,
                            **kwargs):
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments

        """
        xlabel = r'KFreq [rad]; $ \omega c / (2 V)$'
        ylabel1 = r'Viscous Damping; $g = 2 \gamma $'
        ylabel2 = 'Frequency [Hz]'
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
        return
        #print(self.names)

        # results[imode, ipoint, iresult]
        # 1. NaN all the invalid points
        freq = self.results[:, :, self.ifreq]
        iplot, jplot = np.where(freq == 0.0)
        self.results[iplot, jplot, :] = np.nan
        return

        #-----------------------------------------------------------------------
        # 2. sort the results based on velocity, so we're going low to high
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

    def _get_symbols_colors_from_modes(self, modes: np.ndarray) -> tuple[list[str], list[str]]:
        """
        We need to make sure we have a symbol and color for each mode,
        even if we repeat them.

        For the colors, calculate how many more we need N = ceil(nmodes/ncolors)
        and just duplicate colors N times.
        """
        nmodes = len(modes)
        symbols, colors = _symbols_colors_from_nlines(self._colors, self._symbols, nmodes)
        return symbols, colors

    def plot_vg_vf(self, fig=None, damp_axes=None, freq_axes=None, modes=None,
                   plot_type: str='tas',
                   clear: bool=False, close: bool=False, legend: bool=True,
                   xlim=None, ylim_damping=None, ylim_freq=None,
                   vd_limit=None, damping_limit=None,
                   nopoints: bool=False, noline: bool=False,
                   ncol: int=0,
                   png_filename=None, show: bool=False):
        """
        Make a V-g and V-f plot

        Parameters
        ----------
        fig : plt.Figure; default=None
            the fig object
        damp_axes : plt.Axes; default=None
            the ax object
        freq_axes : plt.Axes; default=None
            the ax object
        fig : plt.Figure; default=None
            the fig object
        modes : list[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        plot_type : str; default='tas'
           tas, eas, alt, kfreq, 1/kfreq, freq, damp, eigr, eigi, q, mach
        legend : bool; default=True
            should the legend be shown
        show : bool; default=None???
            show/don't show the plot

        """
        assert vd_limit is None or isinstance(vd_limit, float_types), vd_limit
        assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit
        #self.fix()
        if fig is None:
            fig = plt.figure() # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        modes, imodes = _get_modes_imodes(self.modes, modes)
        symbols, colors = self._get_symbols_colors_from_modes(modes)

        if nopoints:
            symbols = ['None'] * len(symbols)
        linestyle = 'None' if noline else '-'

        #plot_type = ['tas', 'eas', 'alt', 'kfreq', '1/kfreq', 'freq', 'damp', 'eigr', 'eigi', 'q', 'mach']
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
            label = f'Mode {mode:d}; freq={freq[0]:.3g}'
            damp_axes.plot(vel, damping, color=color, marker=symbol, linestyle=linestyle, label=label)
            freq_axes.plot(vel, freq, color=color, marker=symbol, linestyle=linestyle)

        damp_axes.set_xlabel(xlabel)
        freq_axes.set_xlabel(xlabel)
        damp_axes.set_ylabel(r'Viscous Damping; $g = 2 \gamma $')
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

        _add_damping_limit(plot_type, damp_axes, damping_limit)
        _add_vd_limit(plot_type, damp_axes, freq_axes, vd_limit)

        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)
        if legend:
            damp_axes.legend(fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2, ncol=ncol)
            #damp_axes.legend(fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2, ncol=ncol)
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

    def export_to_csv(self, csv_filename: str, modes: Optional[list[int]]=None) -> None:
        """
        Exports a ZONA .veas file

        Parameters
        ----------
        veas_filename : str
            the filename to write
        modes : list[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N

        *.VEAS

        ' DAMPING & FREQUENCY X-Y PLOT FILE OF PLTVG SETID=     714 FOR FLUTTER/ASE ID=     714 NMODE=   12'
        '  IPOINT,  VEL_1,       VEL_2,       RHO_1,       RHO_2, ...'
        '  1,       0.0000E+00,  0.0000E+00,  0.0000E+00,  0.0000E+00'
        '  2,       4.9374E-01, -1.6398E-03, -5.4768E-04, -2.3136E-04'
        """
        imodes = self._imodes(modes)
        assert len(imodes) > 0, imodes
        headers = ['ipoint']
        for name in self.names:
            headersi = [f'{name}_imode={imode}' for imode in imodes]
            headers += headersi
        #print(headers)

        nspeed = self.results.shape[1]
        ipoint = np.arange(1, nspeed+1, 1, dtype='int32')

        csv_filename2 = _apply_subcase_to_filename(csv_filename, self.subcase)
        with open(csv_filename2, 'w') as csv_file:
            csv_file.write(f'# subcase = {self.subcase}\n')
            for key, unit in self.out_units.items():
                csv_file.write(f'# {key} = {unit}\n')
            csv_file.write('# ' + ','.join(headers) + '\n')

            results = [ipoint]
            for iresult, name in enumerate(self.names):
                for imode in imodes:
                    result = self.results[imode, :, iresult]
                    results.append(result)
            results_array = np.column_stack(results)
            np.savetxt(csv_file, results_array, delimiter=',')

            #for i in range(nspeeds):
            #    damping = self.results[:, i, self.idamping]
            #    asdf

    def export_to_veas(self, veas_filename: str, modes: Optional[list[int]]=None) -> None:
        """
        Exports a ZONA .veas file

        Parameters
        ----------
        veas_filename : str
            the filename to write
        modes : list[int] / int ndarray; (default=None -> all)
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
            gmode = f'   G,MODE{mode:->3}'  # 3 characters; pad with -
            wmode = f' WHZ,MODE{mode:->3}'
            damping_modes.append(gmode)
            omega_modes.append(wmode)

        headers = ['EQUIVALENT V'] + damping_modes + ['EQUIVALENT V'] + omega_modes
        veas_filename2 = _apply_subcase_to_filename(veas_filename, self.subcase)
        with open(veas_filename2, 'w') as veas_file:
            veas_file.write(' DAMPING & FREQUENCY X-Y PLOT FILE OF PLTVG '
                            f'SETID=       1 FOR FLUTTER/ASE ID=       1 NMODE= {nmodes:4d}\n')
            veas_file.write(''.join(headers) + '\n')
            nspeeds = self.results.shape[1]
            for i in range(nspeeds):
                damping = self.results[:, i, self.idamping]
                eas = self.results[0, i, self.ieas]
                omega = self.results[:, i, self.ifreq]  # in Hz (what WHZ means)

                values = [eas] + damping.tolist() + [eas] + omega.tolist()
                str_values = (' %11.4E' % value for value in values)
                veas_file.write(''.join(str_values) + '\n')

    def _imodes(self, modes: Optional[np.ndarray, slice[int] |
                                      tuple[int] | list[int]]) -> np.ndarray:
        """gets the imodes from the modes"""
        if modes is None:
            nmodes = self.results.shape[0]
            imodes = range(nmodes)
        elif isinstance(modes, slice):
            nmodes = self.results.shape[0]
            all_modes = np.arange(0, nmodes, dtype='int32')
            imodes = all_modes[modes]
        elif isinstance(modes, (list, tuple)):
            modes = np.array(modes, dtype='int32')
            imodes = modes - 1
        else:
            imodes = modes - 1
        return imodes

    def _modes_nmodes(self, modes: Optional[Iterable[int]]) -> tuple[Iterable[int], int]:
        """gets the modes and nmodes"""
        if modes is None:
            nmodes = self.results.shape[0]
            modes = range(1, nmodes + 1)
        else:
            nmodes = max(modes)
        return modes, nmodes

    def export_to_f06(self, f06_filename: str,
                      modes: Optional[list[int]]=None,
                      page_stamp: Optional[str]=None,
                      page_num: int=1) -> int:
        f06_filename2 = _apply_subcase_to_filename(f06_filename, self.subcase)
        with open(f06_filename2, 'w') as f06_file:
            page_num = self.export_to_f06_file(
                f06_file, modes=modes,
                page_stamp=page_stamp, page_num=page_num)
        return page_num

    def export_to_f06_file(self, f06_file: str,
                           modes: Optional[list[int]]=None,
                           page_stamp: Optional[str]=None,
                           page_num: int=1) -> int:
        imodes = self._imodes(modes)
        if page_stamp is None:
            page_stamp = 'PAGE %i'
        for imode in imodes:
            #'      MACH 0.0                                                                                                                      '
            f06_file.write(f'0                                                                                                            SUBCASE {self.subcase:d}\n')
            f06_file.write('0                                                       FLUTTER  SUMMARY\n')
            f06_file.write('                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = ASYMMETRIC\n')
            f06_file.write('       POINT = %4i     METHOD = %s\n' % (imode + 1, self.method))
            f06_file.write('\n')
            f06_file.write('\n')

            f06_file.write('    KFREQ          1./KFREQ       DENSITY     MACH NO.      VELOCITY       DAMPING     FREQUENCY      COMPLEX   EIGENVALUE\n')
            for res in self.results[imode, :, :9]:
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
                       modes: Optional[list[int]]=None,
                       xlim: Optional[list[float]]=None,
                       plot_type: str='tas',
                       damping_ratios: Optional[list[float]]=None) -> str:
        """
        Writes a custom ZONA flutter file

        Parameters
        ----------
        zona_filename : str
            the filename to write
        modes : list[int] / int ndarray; (default=None -> all)
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
        zona_filename2 = _apply_subcase_to_filename(zona_filename, self.subcase)
        with open(zona_filename2, 'w') as zona_file:
            zona_file.write(msg)
        return msg

    def _plot_type_to_ix_xlabel(self, plot_type: str) -> tuple[int, str]:
        """helper method for ``plot_vg_vf``"""
        plot_type = plot_type.lower()
        if plot_type == 'tas':
            ix = self.ivelocity
            velocity_units = self.out_units['velocity']
            xlabel = f'Velocity [{velocity_units}]'
        elif plot_type == 'eas':
            ix = self.ieas
            velocity_units = self.out_units['eas']
            xlabel = f'Equivalent Airspeed [{velocity_units}]'
        elif plot_type == 'alt':
            ix = self.ialt
            alt_units = self.out_units['altitude']
            xlabel = f'Altitude [{alt_units}]'
        elif plot_type == 'kfreq':
            ix = self.ikfreq
            xlabel = r'Reduced Frequency [rad]; $\omega c / (2V) $'
        elif plot_type == 'rho':
            ix = self.idensity
            density_units = self.out_units['density']
            xlabel = f'Density [{density_units}]'
        elif plot_type == 'q':
            ix = self.iq
            pressure_unit = self.out_units['dynamic_pressure']
            xlabel = f'Dynamic Pressure [{pressure_unit}]'
        elif plot_type == 'mach':
            ix = self.imach
            xlabel = 'Mach'
        elif plot_type == 'freq':
            ix = self.ifreq
            xlabel = 'Frequency [Hz]'
        elif plot_type in ['1/kfreq', 'ikfreq']:
            ix = self.ikfreq_inv
            xlabel = r'1/KFreq [1/rad]; $2V / (\omega c) $'
        elif plot_type == 'eigr':
            ix = self.ieigr
            xlabel = r'Eigenvalue (Real); $\omega \gamma$'
        elif plot_type == 'eigi':
            ix = self.ieigi
            xlabel = r'Eigenvalue (Imaginary); $\omega$'
        elif plot_type in ['damp', 'damping']:
            ix = self.idamping
            xlabel = r'Viscous Damping; $g = 2 \gamma $'
        else:  # pramga: no cover
            raise NotImplementedError("plot_type=%r not in ['tas', 'eas', 'alt', 'kfreq', "
                                      "'1/kfreq', 'freq', 'damp', 'eigr', 'eigi', 'q', 'mach', 'alt']")
        return ix, xlabel

    def object_attributes(self, mode: str='public', keys_to_skip=None,
                          filter_properties: bool=False):
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
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : list[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public', keys_to_skip=None):
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
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : list[str]
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
    elif len(modes) == 0:  # pragma: no cover
        raise RuntimeError('modes = %s' % modes)
    else:
        assert 0 not in modes, modes
        modes = np.unique(modes)
    assert 0 not in modes, modes

    if modes.max() > all_modes.max():
        imodes = np.where(modes <= all_modes.max())
        modes = modes[imodes]
    if len(modes) == 0:  # pragma: no cover
        raise RuntimeError('No modes to plot...')
    imodes = np.searchsorted(all_modes, modes)
    return modes, imodes


def _asarray(results, allow_fix_kfreq: bool=True):
    """casts the results array"""
    allow_fix_kfreq = False
    if not allow_fix_kfreq:
        results = np.asarray(results, dtype='float64')
        return results

    try:
        results = np.asarray(results, dtype='float64')
    except Exception:
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
                        except Exception:  # pragma: no cover
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


def _add_damping_limit(plot_type: str,
                      damp_axes: Axes,
                      damping_limit: Optional[float],
                      linewidth: int=2) -> None:
    if damping_limit is None:
        return
    #damp_label = f'Damping={damping_limit*100:.1f}'
    #plt.axhline(y=1.0, color="black", linestyle="--")
    damp_axes.axhline(y=0., color='k', linestyle='--', linewidth=linewidth,
                      label=f'Viscous Damping=0%')
    damp_axes.axhline(y=damping_limit, color='k', linestyle='-', linewidth=linewidth,
                      label=f'Abs Viscous Damping={damping_limit*100:.1f}%')


def _add_vd_limit(plot_type: str,
                  damp_axes: Axes, freq_axes: Axes,
                  vd_limit: Optional[float],
                  #linestyle: str='--',
                  linewidth: int=2) -> None:
    if vd_limit is None or plot_type not in {'tas', 'eas'}:
        return
    # ax.text(0, vd_limit, 'Damping Limit')
    label1 = f'Vd={vd_limit:.1f}'
    damp_axes.axvline(x=vd_limit, color='k', linestyle='--',
                      linewidth=linewidth, label=label1)
    freq_axes.axvline(x=vd_limit, color='k', linestyle='--',
                      linewidth=linewidth)

    vd_limit_115 = 1.15 * vd_limit
    label2 = f'1.15 Vd={vd_limit_115:.1f}'
    damp_axes.axvline(x=vd_limit_115, color='k', linestyle='-',
                      linewidth=linewidth, label=label2)
    freq_axes.axvline(x=vd_limit_115, color='k', linestyle='-',
                      linewidth=linewidth)


def get_flutter_units(units: Optional[str | dict[str, str]]) -> Optional[str | dict[str, str]]:
    """gets the units"""
    if units is None:
        units = 'english_in'
        #units = {
            #'velocity' : 'in/s', 'density' : 'slug/ft^3',
            #'altitude' : 'ft', 'dynamic_pressure' : 'psf', 'eas':'ft/s'}

    if isinstance(units, str):
        units = units.lower()
        # https://www.dynasupport.com/howtos/general/consistent-units
        # mm, Mg, s / si_ton
        # mm, Mg, s
        #units = {
            #'velocity' : 'mm/s', 'density' : 'Mg/mm^3',
            #'altitude' : 'm', 'dynamic_pressure' : 'MPa', 'eas':'m/s'}

        # units should be consistent
        # what's going on with altitude having inconsistent units?
        if units == 'si':
            units = {'velocity': 'm/s', 'density': 'kg/m^3',
                     'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'm/s'}
        elif units == 'si_mm':
            units = {'velocity': 'mm/s', 'density': 'Mg/mm^3',
                     'altitude': 'm', 'dynamic_pressure': 'MPa', 'eas': 'mm/s'}
        #elif units == 'si_cmgs':
            #units = {'velocity': 'cm/s', 'density': 'g/cm^3',
                     #'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'cm/s'}
        elif units == 'english_in':
            units = {'velocity': 'in/s', 'density': 'slinch/in^3',
                     'altitude': 'ft', 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        elif units == 'english_ft':
            units = {'velocity': 'ft/s', 'density': 'slug/ft^3',
                     'altitude': 'ft', 'dynamic_pressure': 'psf', 'eas': 'ft/s'}
        elif units == 'english_kt':
            units = {'velocity': 'knots', 'density': 'slug/ft^3',
                     'altitude': 'ft', 'dynamic_pressure': 'psf', 'eas': 'knots'}
        else:  # pragma: no cover
            raise NotImplementedError(f'units={units!r} must be in [si, si_mm, '
                                      'english_in, english_ft, english_kt]')
    else:  # pragma: no cover
        assert isinstance(units, dict), f'units={units!r}'
    return units

def _apply_subcase_to_filename(filename: str, subcase: int) -> str:
    """helper for filename management"""
    filename_out = filename
    if '%' in filename:
        filename_out = filename % subcase
    return filename_out

def _update_ncol(nmodes: int, ncol: int=0) -> int:
    """Updates ncol to be a valid number"""
    if ncol > 0:
        return ncol
    nmodes_per_column = 40
    ncol = nmodes // nmodes_per_column
    if nmodes % nmodes_per_column > 0:
        ncol += 1
    ncol = min(ncol, nmodes)
    return ncol


def _symbols_colors_from_nlines(colors: list[str], symbols: list[str],
                                nlines: int) -> tuple[list[str], list[str]]:
    """
    Repeats colors/symbols if there are too many lines

    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
              'C7', 'C8', 'C9', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0', 'C1', 'C2', 'C3',
              'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0',
              'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7',
              'C8', 'C9']
    symbols = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*',
               'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v',
               '>', '>', '>', '>', '>', '>', '>', '>', '>', '>', '<', '<', '<', '<', '<', '<', '<', '<', '<', '<',
               '^', '^', '^', '^', '^', '^', '^', '^', '^', '^']
    """
    ncolors = len(colors)
    nsymbols = len(symbols)
    if ncolors < nlines:
        kcolor = int(np.ceil(nlines / ncolors))
        colors = colors * kcolor

    if nsymbols < nlines:
        ksymbol = int(np.ceil(nlines / nsymbols))
        symbols = symbols * ksymbol
    return symbols, colors
