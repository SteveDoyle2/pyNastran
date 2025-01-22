"""
TODO:
 - more control of font sizes
 - control over point size / text annotation size
 -
"""
from __future__ import annotations
from copy import deepcopy
import warnings
from itertools import count
from typing import Iterable, TextIO, Optional, Any, TYPE_CHECKING

import numpy as np
import scipy
import scipy.interpolate
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

from pyNastran.utils.atmosphere import (
    get_alt_for_density, atm_density,
    convert_altitude, convert_velocity, convert_density, convert_pressure,
)
from pyNastran.utils import object_attributes, object_methods, PathLike
from pyNastran.utils.numpy_utils import float_types

if TYPE_CHECKING and IS_MATPLOTLIB:  # pragma: no cover
    from matplotlib.axes import Axes

Color = str
LineData = tuple[str, float, Color, str]
LINESTYLES = ['-', '--', '-.', ':', 'None', ' ', '',
              'solid', 'dashed', 'dashdot', 'dotted']
Crossing = tuple[float, float, float]
Limit = tuple[Optional[float], Optional[float]] | None


class FlutterResponse:
    """storage object for single subcase SOL 145 results"""

    def __eq__(self, flutter_response: FlutterResponse) -> bool:
        return True
    def export_to_hdf5(self, h5_file, encoding: str):
        return #h5_file.create_dataset('data', data=np.ones(10))

    def __repr__(self) -> str:
        #from pyNastran.utils import object_stats
        #print(object_stats(self))
        #configuration : 'AEROSG2D'
        #in_units : {'velocity': 'm/s', 'density': 'kg/m^3', 'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'm/s'}
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
        #method : 'PKNL'
        #modes  : array([1, 2])
        #names  : ['kfreq', '1/kfreq', 'density', 'velocity', 'damping', 'freq', 'eigr', 'eigi', 'eas', 'q', 'alt']
        #noline : False
        #out_units : {'velocity': 'in/s', 'density': 'slinch/in^3', 'altitude': 'ft',
        # 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        xyz_sym = ''
        if hasattr(self, 'xysym'):
            xyz_sym += f'xysym  = {self.xysym!r}\n'
            xyz_sym += f'xzsym  = {self.xzsym!r}\n'

        msg = (
            'FlutterResponse:\n'
            f'subcase= {self.subcase:d}\n'
            f'{xyz_sym}'
            f'in_units   = {self.in_units}\n'
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
        )
        return msg

    def get_stats(self) -> str:
        return f'FlutterResponse(isubcase={self.subcase})'

    @classmethod
    def from_nx(cls, method: str, fdata: np.ndarray,
                subcase_id: int=1, subtitle: str='', label: str='',
                cref: float=1.0,
                is_xysym: bool=False, is_xzsym: bool=False,
                in_units: dict[str, str]=None):
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
        if in_units is None:
            in_units = get_flutter_units('english_in')
        b = cref / 2.0
        configuration = 'AEROSG2D'  #  TODO: what is this?
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

        resp = FlutterResponse(
            subcase_id, configuration, xysym, xzsym, mach, density_ratio,
            method, modes, results, in_units=in_units,
            subtitle=subtitle,
            label=label,
            use_rhoref=False,
            make_alt=False)
        if 0:  # pragma: no cover
            resp.plot_root_locus(modes=None, fig=None, axes=None, xlim=None, ylim=None,
                                 ncol=ncol, show=False,
                                 clear=False, close=False, legend=True, png_filename=None)
            resp.plot_vg_vf(fig=None, damp_axes=None, freq_axes=None, modes=None, plot_type='eas',
                            legend=True, xlim=None,
                            ylim_damping=None, ylim_freq=None,
                            damping_limit=None,
                            clear=False, close=False,
                            ncol=ncol, png_filename=None, show=False)
            resp.plot_kfreq_damping(modes=None, plot_type='tas', fig=None, damp_axes=None, freq_axes=None,
                                    xlim=None, ncol=ncol, legend=True,
                                    show=False, clear=False, close=False, png_filename=None,
                                    ylim_damping=None, ylim_kfreq=None, damping_limit=None)
            resp.plot_kfreq_damping2(modes=None, fig=None, xlim=None, ylim=None,
                                     show=True, clear=False, close=False,
                                     ncol=ncol, legend=True, png_filename=None)
        return resp

    def __init__(self, subcase: int, configuration: str,
                 xysym: str, xzsym: str,
                 mach: float, density_ratio: float, method: str,
                 modes: list[int], results: Any,
                 in_units: None | str | dict[str, str]=None,
                 subtitle: str='', label: str='',
                 use_rhoref: bool | float=False,
                 make_alt: bool=True,
                 eigenvector: Optional[np.ndarray]=None,
                 eigr_eigi_velocity: Optional[np.ndarray]=None) -> None:
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
        use_rhoref: bool; default=False
            False: assume the density in the table is absolute density
            True: assume the density should be defined by sea level density,
                  so density is a density ratio
            float: use specified rho_ref
        in_units : dict[str] = str (default=None -> no units conversion)
            The units to read from the F06.  Units are only applicable for
            quantities that have units (e.g. units on Mach don't do anything).
            All allowable fields are shown below.

            PK method:
                in_units = {'velocity' : 'in/s'}
                The velocity units are the units for the FLFACT card in the BDF
            PKNL method:
                in_units = {'velocity' : 'in/s', 'density' : 'slinch/in^3', 'altitude' : 'ft', 'dynamic_pressure': 'psi'}
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
        if eigr_eigi_velocity is None:
            eigr_eigi_velocity = np.zeros((0, 3), dtype='float64')
        if eigenvector is None:
            eigenvector = np.array([], dtype='complex128')
        else:
            assert eigenvector.ndim == 2, eigenvector.shape

        assert eigr_eigi_velocity.ndim == 2, eigr_eigi_velocity
        assert eigr_eigi_velocity.shape[1] == 3, eigr_eigi_velocity
        self.eigenvector = eigenvector
        self.eigr_eigi_velocity = eigr_eigi_velocity

        self.in_units = in_units
        self.out_units = ''
        self.make_alt = make_alt
        in_units_dict = get_flutter_units(in_units)

        self.is_complex = True
        self.is_real = True
        self.nonlinear_factor = None
        self.subcase = subcase
        self.subtitle = subtitle
        self.label = label
        self.configuration = configuration
        if method == 'PK':
            self.mach = mach
            self.xysym = xysym
            self.xzsym = xzsym
            self.density_ratio = density_ratio

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
            #flutter.results[:, :, flutter.idensity] *= 1.146e-7
            #print(f'use_rhoref = {use_rhoref}')
            if isinstance(use_rhoref, float_types):
                rhoref = use_rhoref
                results[:, :, self.idensity] *= rhoref
            elif use_rhoref:
                rhoref = atm_density(alt=0, density_units=in_units_dict['density'])
                #print(f'applying rhoref={rhoref}')
                results[:, :, self.idensity] *= rhoref
            else:
                assert isinstance(use_rhoref, bool) or isinstance(use_rhoref, float_types), use_rhoref

            # print('rho2', results[0, :, self.idensity])
            # print('vel2', results[0, :, self.ivelocity])
            # print('mach2', results[0, :, self.imach])

        else:  # pragma: no cover
            raise NotImplementedError(method)

        #-------------------------------------------
        if self.method in ['PK', 'KE']:
            pass
        elif self.method == 'PKNL':
            results = self._set_pknl_results(in_units_dict, results)
        else:  # pragma: no cover
            raise NotImplementedError(self.method)

        if len(self.eigenvector) > 1:
            #print(f'results.shape = {results.shape}')
            nmodes, nvelocity = results.shape[:2]
            #print(f'nmodes={nmodes}; nvelocity={nvelocity}')
            # not sure on the shapes
            #
            # npoint, nvelocity, 3?
            #self.eigr_eigi_velocity = self.eigenvector.reshape(nmodes, nvelocity, 3)
            # nmodes_mpf, npoint, nvelocity?
            eigenvector, eigr_eigi_velocity = reshape_eigenvectors(
                self.eigenvector, self.eigr_eigi_velocity)
            self.eigenvector = eigenvector
            self.eigr_eigi_velocity = eigr_eigi_velocity
            assert len(self.eigenvector) and len(self.eigr_eigi_velocity), (len(self.eigenvector), len(self.eigr_eigi_velocity))
            #= self.eigenvector.reshape(nmodes, nmodes, nvelocity)

            nvelocity = results.shape

        self.results_in = results
        self.results = self.results_in

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
        self.nopoints = False
        self._symbols: list[str] = []
        self._colors: list[str] = []
        self.generate_symbols()
        self.set_symbol_settings()
        self.set_font_settings(font_size=None)

        self._xtick_major_locator_multiple = None
        self._ytick_major_locator_multiple = None

    def set_plot_settings(self, figsize=None,
                          xtick_major_locator_multiple=None,
                          ytick_major_locator_multiple=None):
        """
        Parameters
        ----------
        figsize : tuple[float, float]; default -> (6.4, 4.8)
            the size of the figure
        x/ytick_major_locator_multiple : list[float]; default=None
            the delta spacing for the x/y axis
        """
        if figsize is None:
            figsize = plt.rcParams['figure.figsize']
            # figsize = (6.4, 4.8)

        #print('plt.rcParams = ', plt.rcParams)
        assert xtick_major_locator_multiple is None or isinstance(xtick_major_locator_multiple, (list, tuple)), xtick_major_locator_multiple
        assert ytick_major_locator_multiple is None or isinstance(ytick_major_locator_multiple, (list, tuple)), ytick_major_locator_multiple
        self._xtick_major_locator_multiple = xtick_major_locator_multiple
        self._ytick_major_locator_multiple = ytick_major_locator_multiple
        figure_params = {
            'figsize': figsize,
        }
        matplotlib.rc('figure', **figure_params)

    def set_symbol_settings(self, nopoints: bool=False,
                            show_mode_number: bool=False,
                            point_spacing: int=0,
                            markersize: int=0) -> None:
        if markersize is None:
            markersize = plt.rcParams['lines.markersize']
        #if point_spacing is None:
            #point_spacing = plt.rcParams['lines.markevery']

        #print(list(plt.rcParams.keys()))
        plt.rcParams['lines.markersize'] = markersize
        #plt.rcParams['lines.markevery'] = point_spacing
        #out = {'lines': markersize}
        #matplotlib.rc('lines', **out)

        self.nopoints = nopoints
        self.show_mode_number = show_mode_number
        self.markevery = point_spacing

    def set_font_settings(self, font_size: Optional[int]=13) -> None:
        # TODO: split label, title, tic_marks, legend
        if font_size is None:
            font_size = plt.rcParams['font.size']

        #print(f'font_size = {font_size}')
        font = {
            #'family': 'normal',
            #'weight': 'bold',
            'size': font_size,
        }
        matplotlib.rc('font', **font)
        self.font_size = font_size

    def set_out_units(self, out_units: str | dict[str, str]) -> None:
        self.convert_units(out_units)

    def convert_units(self, out_units: Optional[str | dict[str, str]]) -> None:
        results = self.results_in.copy()
        out_units_dict = get_flutter_units(out_units)
        kvel = _get_unit_factor(self.in_units, out_units_dict, 'velocity')[0]
        keas = _get_unit_factor(self.in_units, out_units_dict, 'eas')[0]
        kdensity = _get_unit_factor(self.in_units, out_units_dict, 'density')[0]
        kpressure = _get_unit_factor(self.in_units, out_units_dict, 'dynamic_pressure')[0]
        kalt = _get_unit_factor(self.in_units, out_units_dict, 'altitude')[0]

        # (imode, istep, iresult)
        results[:, :, self.ivelocity] *= kvel
        if self.method in ['PK', 'KE']:
            pass
        elif self.method == 'PKNL':
            results[:, :, self.idensity] *= kdensity
            results[:, :, self.iq] *= kpressure
            results[:, :, self.ieas] *= keas
            results[:, :, self.ialt] *= kalt
        else:  # pragma: no cover
            raise NotImplementedError(self.method)
        self.out_units = out_units_dict
        self.results = results

    def _set_pknl_results(self,
                          in_units: dict[str, str],
                          results: np.ndarray) -> np.ndarray:
        assert results.shape[2] in [9, 11, 12], results.shape
        density_units_in = in_units['density']

        # in/s
        vel = results[:, :, self.ivelocity]  #.ravel()

        # slinch/in^3 - in_units
        rho = results[:, :, self.idensity]  #.ravel()

        # good
        rho_ref = atm_density(0., R=1716., alt_units='ft',
                              density_units=density_units_in)

        vel_units_in = in_units['velocity']
        q_units_in = in_units['dynamic_pressure']
        if _is_q_units_consistent(density_units_in, vel_units_in, q_units_in):
            q = 0.5 * rho * vel**2
        else:  # pragma: no cover
            raise NotImplementedError((density_units_in, vel_units_in, q_units_in))

        #eas  = (2 * q / rho_ref)**0.5
        # eas = V * sqrt(rho / rhoSL)
        eas = vel * np.sqrt(rho / rho_ref)
        #print('vel = ', vel)
        #print('rho = ', rho)
        #print('rho_ref = ', rho_ref)

        altitude_units = in_units['altitude']

        #print('density_units_in=%r density_units2=%r' % (density_units_in, density_units2))
        kdensityi = convert_density(1., density_units_in, 'slug/ft^3')

        resultsi = results[:, :, :9]
        assert resultsi.shape[2] == 9, resultsi.shape

        rho_in_slug_ft3 = rho * kdensityi
        ft_to_alt_unit = convert_altitude(1., 'ft', altitude_units)
        if self.make_alt:
            alt_ft = []
            for idensity, densityi in enumerate(rho_in_slug_ft3.ravel()):
                try:
                    alt_fti = get_alt_for_density(densityi, density_units='slug/ft^3',
                                                  alt_units='ft', nmax=20)
                except Exception:
                    raise RuntimeError(f'failed to find altitude for density[{idensity}]='
                                       f'{rho.ravel()[idensity]:g}')
                alt_ft.append(alt_fti)
            alt = np.array(alt_ft, dtype=rho.dtype).reshape(vel.shape) * ft_to_alt_unit
        else:
            alt = np.full(vel.shape, np.nan, dtype=vel.dtype)

        results2 = np.dstack([resultsi, eas, q, alt])
        #results2[:, :, self.idensity] = rho
        return results2

    def generate_symbols(self, colors=None, symbols=None,
                         imethod: int=0):
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
            #colors = ['r', 'g', 'b', 'k', 'm']  # 5
            colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']  # 10

        if symbols is None:
            symbols = ['o', '*', 'x', 'v', '>', '<', '^']  # 7
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

    def plot_complex_modes(self):
        fig2 = plt.figure(2)
        fig2.clear()
        axes = fig2.subplots(nrows=2, ncols=2)
        assert axes.shape == (2, 2), axes.shape
        ax11 = axes[0, 0]
        ax12 = axes[1, 0]

        ax21 = axes[1, 0]
        ax22 = axes[1, 1]
        #print(self.eigenvector)
        nmodes1, nmodes2 = self.eigenvector.shape
        nmodes = max(nmodes1, nmodes2)

        for imode1 in range(nmodes1):
            ax11.plot(self.eigenvector[imode1, :].real, label=f'iMode1={imode1+1}')
            ax12.plot(self.eigenvector[imode1, :].imag, label=f'iMode1={imode1+1}')

        for imode2 in range(nmodes2):
            ax21.plot(self.eigenvector[:, imode2].real, label=f'iMode2={imode2+1}')
            ax22.plot(self.eigenvector[:, imode2].imag, label=f'iMode2={imode2+1}')
        plt.legend()

        ax11.set_ylabel('Real')  #, fontsize=self.font_size)
        ax12.set_ylabel('Imaginary')  #, fontsize=self.font_size)
        for ax in axes.ravel():
            ax.grid(True)
            ax.set_xlim([0, nmodes])
        # ax2.set_ylim([0, nmodes])
        plt.show()

    def plot_vg(self, fig=None, modes=None,
                plot_type: str='tas',
                xlim=None, ylim_damping=None,
                v_lines: Optional[list[LineData]]=None,
                ncol: int=0,
                clear: bool=False, legend: bool=True,
                freq_tol: float=-1.0,
                png_filename=None, show: bool=True, **legend_kwargs):
        """
        Make a V-g plot

        See ``plot_root_locus`` for arguments
        """
        ix, xlabel, xunit = self._plot_type_to_ix_xlabel(plot_type)
        ylabel = 'Structural Damping'
        iy = self.idamping
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, xlim=xlim, ylim=ylim_damping,
                       ncol=ncol,
                       show=show, clear=clear, legend=legend,
                       freq_tol=freq_tol,
                       v_lines=v_lines, plot_type=plot_type, xunit=xunit,
                       png_filename=png_filename,
                       **legend_kwargs)

    #@property
    def flutter_speed(self, modes=None,
                      dfreq: float=-1.0,
                      #ddamp: float=-1.0,
                      damping_range: Limit=None,
                      velocity_range: Limit=None):
        """gets the flutter speed"""
        if damping_range is None:
            damping_range = [None, None]
        if velocity_range is None:
            velocity_range = [None, None]

        is_damping_range = damping_range[0] is not None or damping_range[1] is not None
        is_velocity_range = velocity_range[0] is not None or velocity_range[1] is not None

    def get_flutter_crossings(self,
                              damping_required: Optional[list[tuple[float, float]]]=None,
                              modes=None,
                              eas_range: Optional[tuple[float, float]]=None,
                              freq_round: int=2,
                              eas_round: int=3,
                              ) -> dict[int, list[Crossing]]:
        """
        Gets the flutter crossings

        Parameters
        ----------
        damping_required: list[(damping_target, damping_required)]
            A point will be created at the damping point if it meets the criterion
            for damping_required. This is useful when you have slight spikes over the target damping
            target_damping: float
                the target damping for flutter
            damping_required: float
                the required damping for flutter
        modes : (nmode,) int ndarray; default=None -> all modes
            the modes to analyze
        eas_range: [eas_min, eas_max]
            the equivalent airspeed range of the plot to consider
        freq_round: int; default=2
            number of digits for the frequency
        eas_round
            number of digits for the equivalent airspeed

        Returns
        -------
        xcrossing_dict: dict[mode, [(percent_damping, freq, velocity), ...]
            the flutter crossings
        """
        if damping_required is None:
            damping_required = [
                (0.00, 0.01),
                (0.03, 0.03),
            ]
        if eas_range is None:
            eas_range = [None, None]

        eas_min0, eas_max0 = eas_range
        modes, imodes = _get_modes_imodes(self.modes, modes)
        min_damping = _get_min_damping(damping_required)

        # if dfreq > 0. and is_damping_range:
        #     raise NotImplementedError('dfreq > 0. and ddamp > 0.')
        # elif is_damping_range:
        xcrossing_dict = {}
        for imode in imodes:
            mode = imode + 1
            dampi = self.results[imode, :, self.idamping].flatten()
            if dampi.max() < min_damping:
                continue
            easi = self.results[imode, :, self.ieas].flatten()
            freqi = self.results[imode, :, self.ifreq].flatten()

            crossings = []
            for damping_targeti, damping_requiredi in damping_required:
                if dampi.max() < damping_requiredi:
                    continue
                eas0, freq0 = get_zero_crossings(easi, freqi, dampi - damping_targeti)
                freq0 = round(freq0, freq_round)
                eas0 = round(eas0, eas_round)

                eas0, freq0 = check_range(eas_min0, eas_max0, freq0, eas0)
                if np.isnan(eas0):
                    continue
                crossings.append((damping_targeti, float(freq0), float(eas0)))

            if len(crossings) == 0:
                continue
            xcrossing_dict[mode] = crossings
        return xcrossing_dict

    def plot_root_locus(self, modes=None,
                        fig=None, axes=None,
                        eigr_lim=None, eigi_lim=None,
                        ivelocity=None,
                        ncol: int=0,
                        show: bool=True, clear: bool=False,
                        close: bool=False, legend: bool=True,
                        #noline: bool=False, nopoints: bool=False,
                        freq_tol: float=-1.0,
                        png_filename=None,
                        **legend_kwargs):
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
        eigr_lim : list[float/None, float/None]
            the x plot limits
        eigi_lim : list[float/None, float/None]
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
        assert isinstance(freq_tol, float_types), freq_tol
        xlabel = r'Eigenvalue (Real); $\omega \zeta$'
        ylabel = r'Eigenvalue (Imaginary); $\omega$'
        ix = self.ieigr
        iy = self.ieigi
        scatter = True
        self._plot_x_y(ix, iy, xlabel, ylabel, scatter,
                       modes=modes, fig=fig, axes=axes,
                       ivelocity=ivelocity,
                       xlim=eigr_lim, ylim=eigi_lim,
                       ncol=ncol,
                       show=show, clear=clear, close=close, legend=legend,
                       freq_tol=freq_tol,
                       png_filename=png_filename,
                       **legend_kwargs)

    def plot_modal_participation(self, ivel: int, mode: int,
                                 modes=None,
                                 mag_tol: float=-1.0,
                                 fig=None, axes=None,
                                 #eigr_lim=None, eigi_lim=None,
                                 ncol: int=0,
                                 show: bool=True, clear: bool=False,
                                 close: bool=False, legend: bool=True,
                                 #noline: bool=False, nopoints: bool=False,
                                 freq_tol: float=-1.0,
                                 png_filename=None,
                                 **legend_kwargs):
        """
        Plots the modal participations factors

        Parameters
        ----------
        modes : list[int] / int ndarray; (default=None -> all)
            the modes; typically 1 to N
        fig : plt.Figure
            the figure object
        axes : plt.Axes
            the axes object
        mag_tol: float; default=-1.0
            filters eigenvectors that have a magnitude less than a target
            <= 0.0 is ineffective
        eigr_lim : list[float/None, float/None]
            the x plot limits
        eigi_lim : list[float/None, float/None]
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
        assert isinstance(ivel, int), ivel
        assert isinstance(mode, int), mode
        # if not isinstance(mag_tol, float_types):
        #     warnings.warn(f'mag_tol={mag_tol!r} is not a float; default=-1.0')
        #     mag_tol = -1.0
        assert isinstance(mag_tol, float_types), mag_tol
        modes, imodes = _get_modes_imodes(self.modes, modes)
        legend_kwargs = get_legend_kwargs(self.font_size, legend_kwargs)
        if fig is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)

        xlabel = r'Eigenvalue (Real); $\omega \zeta$'
        ylabel = r'Eigenvalue (Imaginary); $\omega$'
        # ix = self.ieigr
        # iy = self.ieigi
        # scatter = True

        #print(f'eigr_eigi_velocity.shape: {self.eigr_eigi_velocity.shape}')
        #print(f'eigenvector.shape: {self.eigenvector.shape}')
        #print(f'eigr_eigi_velocity:\n{self.eigr_eigi_velocity}')

        # MSC
        #eigr_eigi_velocity:
        #[[-9.88553e-02  1.71977e+01  1.52383e+02]
        # [-1.71903e-01  6.60547e+01  1.52383e+02]]
        nvel, nmode = self.eigr_eigi_velocity.shape[:2]
        assert ivel < nvel, f'ivel={ivel} nvel={nvel}'
        assert mode <= nmode, f'mode={mode} nmode={nmode}'

        imode = mode - 1
        assert isinstance(self.eigr_eigi_velocity, np.ndarray), type(self.eigr_eigi_velocity)
        #print('self.eigr_eigi_velocity = ', self.eigr_eigi_velocity)
        try:
            eigr_eigi_velocity = self.eigr_eigi_velocity[ivel, imode, :]
        except IndexError:
            raise RuntimeError(f'eigr_eigi_velocity.shape=(ivel, imode, :)={self.eigr_eigi_velocity.shape}; ivel={ivel} nvel={nvel} imode={imode}')
        eigri, eigii, velocityi = eigr_eigi_velocity

        omega_damping = eigri
        omega = eigii

        #xlabel = r'Eigenvalue (Real); $\omega \gamma$'
        #xlabel = r'Eigenvalue (Imaginary); $\omega$'
        freq = omega / (2 * np.pi)

        # filter division by 0
        damping_g = np.zeros(omega.shape, dtype=omega.dtype)
        iomega = (omega != 0.0)
        damping_g[iomega] = 2 * omega_damping[iomega] / omega[iomega]

        title = f'Subcas {self.subcase}; Modal Participation Factors of Mode {mode}\n'
        title += rf'$\omega$={omega:.2f}; freq={freq:.2f} Hz; g={damping_g:.3g}'
        if np.isfinite(velocityi):
            title += f' V={velocityi:.1f}'
        #print(title)
        axes.set_title(title)  #, fontsize=self.font_size)
        axes.grid(True)
        axes.set_xlabel(xlabel)  #, fontsize=self.font_size)
        axes.set_ylabel(ylabel)  #, fontsize=self.font_size)

        _set_ticks(self, axes, 0)

        #print(f'eigr_eigi_velocity:\n{self.eigr_eigi_velocity}')
        #print(f'eigenvector:\n{self.eigenvector}')

        eig = self.eigenvector[ivel, imode, :]
        #print(f'imodes = {imodes}')
        eigr = eig.real[imodes]
        eigi = eig.imag[imodes]
        #abs_eigr = np.linalg.norm(eigr)
        #abs_eigi = np.linalg.norm(eigi)
        #if abs_eigr == 0.0:
        #    abs_eigr = 1.0
        #if abs_eigi == 0.0:
        #    abs_eigi = 1.0
        mag = np.sqrt(eigr**2 + eigi**2)
        #print(f'eigr = {eigr}')
        #print(f'mag = {mag}')

        # TODO: is this the right axis?
        mag_max = mag.max()  # was imodes
        #print(f'mag_max = {mag_max}')

        # normalize the magnitude, so it's a percentage
        # TODO: is this the right axis?
        mag /= mag_max # [:, np.newaxis]
        #print(f'mag2 = {mag}')
        #print(f'mag_tol = {mag_tol!r}')

        # filter out row (or col) for points that are
        # too close to 0
        ifilter0 = (mag > mag_tol)
        #print(f'ifilter0 = {ifilter0}')
        # TODO: is this the right axis?
        #ifilter = np.any(ifilter0, axis=1)
        #print(f'ifilter = {ifilter}')

        # don't normalize
        abs_eigr = 1.0
        abs_eigi = 1.0

        # reals = eigr / abs_eigr
        # imags = eigi / abs_eigi

        reals = eigr / abs_eigr
        imags = eigi / abs_eigi
        #print(f'reals = {reals}')

        symbols, colors = self._get_symbols_colors_from_modes(modes)
        jcolor = 0
        for i, imodei, mode in zip(count(), imodes, modes):
            symbol = symbols[jcolor]
            color = colors[jcolor]

            reali = reals[i]
            imagi = imags[i]
            magi = mag[i]
            if magi < mag_tol:
                continue
            text = str(mode)
            eig_str = f'{reali:.2g}+{imagi:.2g}j; A={magi:.2g}'.replace('+-', '-')
            label = f'Mode {mode}; {eig_str}'
            #print(label)
            axes.scatter(reali, imagi, label=label, alpha=0.7)
            #print(f'{i}: {reali}, {imagi}, {text!r}')
            axes.text(reali, imagi, text, ha='center', va='center')
                      #fontsize=self.font_size)

        if legend:
            # bbox_to_anchor=(1.125, 1.), ncol=ncol,
            axes.legend(**legend_kwargs)

        _show_save_clear_close(
            fig, show, png_filename, clear, close)
    # for imode1 in range(nmodes1):
        #     ax11.plot(self.eigenvector[imode1, :].real, label=f'iMode1={imode1+1}')
        #     ax12.plot(self.eigenvector[imode1, :].imag, label=f'iMode1={imode1+1}')

    def _plot_x_y(self, ix: int, iy: int,
                  xlabel: str, ylabel: str,
                  scatter: bool,
                  modes=None,
                  fig=None, axes=None,
                  v_lines=None, plot_type: str='', xunit: str='',
                  xlim=None, ylim=None,
                  ivelocity=None,
                  ncol: int=0,
                  show: bool=True, clear: bool=False,
                  close: bool=False, legend: bool=True,
                  freq_tol: float=-1.0,
                  png_filename=None,
                  **legend_kwargs):
        """
        builds the plot for:
         - Vg
         - root-locus
        """
        legend_elements = []
        self.fix()
        #print('plot_xy')
        legend_kwargs = get_legend_kwargs(self.font_size, legend_kwargs)

        modes, imodes = _get_modes_imodes(self.modes, modes)
        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)
        #print(f'plot_xy: modes  = {modes}')
        #print(f'plot_xy: imodes = {imodes}')
        #print(f'plot_xy: ncol   = {ncol}')

        if fig is None:
            fig = plt.figure()
            axes = fig.add_subplot(111)

        symbols, colors = self._get_symbols_colors_from_modes(modes)
        linestyle = 'None' if self.noline else '-'

        jcolor = 0
        #print('starting plot part; jcolor=0')
        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[jcolor]
            color = colors[jcolor]
            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            ys = self.results[imode, :, iy].ravel()
            #print('freq, xs, ys')
            jcolor, color2, linestyle2, symbol2, texti = _increment_jcolor(
                mode, jcolor, color, linestyle, symbol,
                freq, freq_tol=freq_tol,
                show_mode_number=self.show_mode_number)
            #print(f'plot_xy: jcolor={jcolor}; color={color2}; linstyle={linestyle2}; symbol={symbol2}')

            #print(f'freq={freq}')
            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)
            label = _get_mode_freq_label(mode, freq[0])
            #print(f'iplot={iplot}')
            assert len(xs[iplot]) > 0
            line = axes.plot(xs[iplot], ys[iplot],
                             color=color2, marker=symbol2, label=label,
                             linestyle=linestyle2, markersize=0)
            legend_elements.append(line)

            if ivelocity and symbol2 and ivelocity < len(xs):
                # single point selected for modes plotting corresponding to ivelocity
                markersize = 10
                plot_kwargs = {
                    'color': 'k', 'marker': 'o',
                    's': markersize**2, 'alpha': 0.8}
                line = axes.scatter(xs[ivelocity], ys[ivelocity], **plot_kwargs)
                legend_elements.append(line)

            if scatter:
                # used for the root locus
                # velocity increases driving an incrase in the point size
                scatteri = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                #axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=symbol[0], marker=symbol[1])
                line = axes.scatter(xs[iplot], ys[iplot], s=scatteri, color=color, marker=symbol2)
                legend_elements.append(line)

        legend_elementsi = _add_vertical_lines(
            [axes], v_lines, plot_type, xunit)
        legend_elements.extend(legend_elementsi)

        axes.grid(True)
        #axes.set_xlabel(xlabel  + '; _plot_x_y', fontsize=self.font_size)
        axes.set_xlabel(xlabel)  #, fontsize=self.font_size)
        axes.set_ylabel(ylabel)  #, fontsize=self.font_size)
        _set_ticks(self, axes, 0)
        set_xlim(axes, xlim)
        set_ylim(axes, ylim)

        title = f'Subcase {self.subcase:d}'
        if png_filename:
            title += '\n%s' % png_filename
        #print(f'title={title!r}')
        fig.suptitle(title)
        if legend:
            # bbox_to_anchor=(1.125, 1.), ncol=ncol,
            axes.legend(**legend_kwargs)

        _show_save_clear_close(
            fig, show, png_filename, clear, close)
        return axes

    def _plot_x_y2(self, ix: int, iy1: int, iy2: int,
                   xlabel: str, ylabel1: str, ylabel2: str,
                   scatter: bool,
                   modes=None,
                   fig=None, axes1=None, axes2=None,
                   xlim=None, ylim1=None, ylim2=None,
                   ncol: int=0,
                   show: bool=True, clear: bool=False,
                   close: bool=False, legend: bool=True,
                   freq_tol: float=-1.0,
                   png_filename=None,
                   **legend_kwargs) -> tuple[plt.Figure,
                                             tuple[plt.Axes, plt.Axes]]:
        """
        Builds the plot with 2 subplots for:
         - Vg, Vf
         - Vg, Vk
         - kg, kf

        Parameters
        ----------
        scatter : bool
            draw the points with growth along the x-axis

        """
        self.fix()
        legend_kwargs = get_legend_kwargs(self.font_size,
                                          legend_kwargs)

        modes, imodes = _get_modes_imodes(self.modes, modes)
        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)

        if fig is None:
            fig = plt.figure()
            gridspeci = gridspec.GridSpec(2, 4)
            axes1 = fig.add_subplot(gridspeci[0, :3])
            axes2 = fig.add_subplot(gridspeci[1, :3], sharex=axes1)
        #else:
            #$self.log.info('got a fig')

        set_xlim(axes1, xlim)
        set_xlim(axes2, xlim)
        set_ylim(axes1, ylim1)
        set_ylim(axes2, ylim2)

        symbols, colors = self._get_symbols_colors_from_modes(modes)
        #self.log.debug(f'symbols={symbols}')
        #self.log.debug(f'colors={colors}')
        linestyle = 'None' if self.noline else '-'
        if self.nopoints:  # and noline is False:
            scatter = False

        markersize = None
        if self.noline:
            markersize = 0

        legend_elements = []
        jcolor = 0
        for i, imode, mode in zip(count(), imodes, modes):
            symbol = symbols[jcolor]
            color = colors[jcolor]

            freq = self.results[imode, :, self.ifreq].ravel()
            xs = self.results[imode, :, ix].ravel()
            y1s = self.results[imode, :, iy1].ravel()
            y2s = self.results[imode, :, iy2].ravel()
            jcolor, color2, linestyle2, symbol2, texti = _increment_jcolor(
                mode, jcolor, color, linestyle, symbol,
                freq, freq_tol=freq_tol,
                show_mode_number=self.show_mode_number)

            iplot = np.where(freq != np.nan)
            #iplot = np.where(freq > 0.0)

            # plot the line
            label = _get_mode_freq_label(mode, freq[0])
            legend_element = Line2D([0], [0], color=color2,
                                    marker=symbol2, label=label, linestyle=linestyle2)
            if self.nopoints:
                symbol2 = 'None'

            # self.log.info(f'scatter={scatter}; color={color2}; linestyle={linestyle2!r} symbol={symbol2!r}; markersize={markersize}')
            # self.log.info(f'xs={xs[iplot]}')
            # self.log.info(f'y1s={y1s[iplot]}')
            # self.log.info(f'y2s={y2s[iplot]}')
            if scatter:
                # when velocity increases, point size increases
                scatteri = np.linspace(.75, 50., len(xs))
                #assert symbol[2] == '-', symbol
                axes1.scatter(xs[iplot], y1s[iplot],
                              s=scatteri, color=color2, marker=symbol2)
                axes2.scatter(xs[iplot], y2s[iplot],
                              s=scatteri, color=color2, marker=symbol2)

                # Draw the line
                axes1.plot(xs[iplot], y1s[iplot], marker=symbol2, label=label,
                           color=color2, markersize=markersize, linestyle=linestyle2)
                axes2.plot(xs[iplot], y2s[iplot], marker=symbol2,
                           color=color2, markersize=markersize, linestyle=linestyle2)
            else:
                legend_elementsi = _plot_two_axes(
                    axes1, axes2,
                    xs[iplot], y1s[iplot], y2s[iplot],
                    color, symbol2, linestyle, label, texti,
                    self.markevery, markersize=markersize,
                )
                # axes1.plot(xs[iplot], y1s[iplot], marker=symbol2, label=label,
                #            color=color2, markersize=markersize, linestyle=linestyle2)
                # axes2.plot(xs[iplot], y2s[iplot], marker=symbol2,
                #            color=color2, markersize=markersize, linestyle=linestyle2)
            legend_elements.append(legend_element)

        axes1.grid(True)
        axes1.set_xlabel(xlabel)  #, fontsize=self.font_size)
        #axes1.set_xlabel(xlabel + '; _plot_x_y2', fontsize=self.font_size)
        axes1.set_ylabel(ylabel1)  #, fontsize=self.font_size)

        axes2.grid(True)
        axes2.set_xlabel(xlabel)  #, fontsize=self.font_size)
        axes2.set_ylabel(ylabel2)  #, fontsize=self.font_size)

        _set_ticks(self, axes1, 0)
        _set_ticks(self, axes2, 1)

        title = f'Subcase {self.subcase:d}'
        if png_filename:
            title += '\n%s' % png_filename
        fig.suptitle(title)
        if legend:
            #legend_kwargs = {}
            #for key, value in legend_kwargs.items():
                #if key in {'ylim'}:
                    #continue
                #legend_kwargs[key] = value

            axes1.legend(handles=legend_elements, **legend_kwargs)              # TODO: switch to figure...
            #axes1.legend(handles=legend_elements, ncol=ncol, **legend_kwargs)  # TODO: switch to figure...
            #fig.legend(handles=legend_elements, ncol=ncol, **legend_kwargs)

        _show_save_clear_close(
            fig, show, png_filename, clear, close)
        return fig, (axes1, axes2)

    def plot_kfreq_damping(self, modes=None,
                           plot_type: str='tas',
                           fig=None, damp_axes=None, freq_axes=None,
                           xlim=None,
                           ylim_damping=None,
                           ylim_kfreq=None,
                           show: bool=True, clear: bool=False,
                           close: bool=False, legend: bool=True,
                           freq_tol: float=-1.0,
                           png_filename=None,
                           damping_limit=None,
                           v_lines: list[LineData]=None,
                           **kwargs) -> tuple[plt.Figure,
                                             tuple[plt.Axes, plt.Axes]]:
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments
        """
        assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit

        ylabel1 = r'Structural Damping; $g = 2 \gamma $'
        ylabel2 = r'KFreq [rad]; $ \omega c / (2 V)$'

        ix, xlabel, unused_xunit = self._plot_type_to_ix_xlabel(plot_type)
        iy1 = self.idamping
        iy2 = self.ikfreq
        scatter = True
        #print(f"plot_kfreq_damping; plot_type={plot_type}")
        fig, axes = self._plot_x_y2(
            ix, iy1, iy2, xlabel, ylabel1, ylabel2, scatter,
            modes=modes, fig=fig, axes1=damp_axes, axes2=freq_axes,
            xlim=xlim, ylim1=ylim_damping, ylim2=ylim_kfreq,
            show=show, clear=clear, close=close,
            legend=legend,
            freq_tol=freq_tol,
            png_filename=png_filename,
            **kwargs)
        return fig, axes

    def plot_kfreq_damping2(self, modes=None,
                            fig=None, damp_axes=None, freq_axes=None,
                            xlim=None, ylim_damping=None, ylim_freq=None,
                            show: bool=True, clear: bool=False,
                            close: bool=False, legend: bool=True,
                            freq_tol: float=-1.0,
                            png_filename=None,
                            **kwargs) -> tuple[plt.Figure,
                                               tuple[plt.Axes, plt.Axes]]:
        """
        Plots a kfreq vs. damping curve

        See ``plot_root_locus`` for arguments

        """
        xlabel = r'KFreq [rad]; $ \omega c / (2 V)$'
        ylabel1 = r'Structural Damping; $g = 2 \gamma $'
        ylabel2 = 'Frequency [Hz]'
        ix = self.ikfreq
        iy1 = self.idamping
        iy2 = self.ifreq
        scatter = True
        fig, axes = self._plot_x_y2(
            ix, iy1, iy2, xlabel, ylabel1, ylabel2, scatter,
            modes=modes,
            fig=fig, axes1=damp_axes, axes2=freq_axes,
            xlim=xlim, ylim1=ylim_damping, ylim2=ylim_freq,
            show=show, clear=clear, close=close,
            legend=legend,
            freq_tol=freq_tol,
            png_filename=png_filename,
            **kwargs)
        return fig, axes

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

    def _get_symbols_colors_from_modes(self, modes: np.ndarray,
                                       nopoints: Optional[bool]=None) -> tuple[list[str], list[str]]:
        """
        We need to make sure we have a symbol and color for each mode,
        even if we repeat them.

        For the colors, calculate how many more we need N = ceil(nmodes/ncolors)
        and just duplicate colors N times.
        """
        if nopoints is None:
            nopoints = self.nopoints
        nmodes = len(modes)
        symbols, colors = _symbols_colors_from_nlines(self._colors, self._symbols, nmodes)
        if nopoints:
            symbols = ['None'] * len(symbols)
        return symbols, colors

    def plot_vg_vf(self, fig=None, damp_axes=None, freq_axes=None, modes=None,
                   plot_type: str='tas',
                   clear: bool=False, close: bool=False, legend: bool=True,
                   xlim: Optional[Limit]=None,
                   ylim_damping: Optional[Limit]=None,
                   ylim_freq: Optional[Limit]=None,
                   ivelocity: Optional[int]=None,
                   v_lines: list[LineData]=None,
                   damping_limit=None,
                   ncol: int=0,
                   freq_tol: float=-1.0,
                   filter_freq: bool=False,
                   damping_required: list[tuple[float, float]]=None,
                   filter_damping: bool=False,
                   eas_range: Optional[tuple[float, float]]=None,
                   png_filename=None, show: bool=False,
                   ) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes]]:
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
        filter_freq : bool; default=False
            filter modes entirely outside the plot range
            useful for cleaning up the legend; False for interactive mode
        damping_required: list[(damping_target, damping_required)]
            A point will be created at the damping point if it meets the criterion
            for damping_required. This is useful when you have slight spikes over the target damping
            target_damping: float
                the target damping for flutter
            damping_required: float
                the required damping for flutter
        filter_damping: bool; default=False
            filter crossings entirely outside the plot range
            useful for cleaning up the legend
        """
        #assert vl_limit is None or isinstance(vl_limit, float_types), vl_limit
        assert damping_limit is None or isinstance(damping_limit, float_types), damping_limit
        #self.fix()
        if fig is None:
            fig = plt.figure()  # figsize=(12,9), self.subcase
            gridspeci = gridspec.GridSpec(2, 4)
            damp_axes = fig.add_subplot(gridspeci[0, :3])
            freq_axes = fig.add_subplot(gridspeci[1, :3], sharex=damp_axes)

        #self._set_xy_limits(xlim, ylim)
        modes, imodes = _get_modes_imodes(self.modes, modes)
        symbols, colors = self._get_symbols_colors_from_modes(modes)
        symbols_show, colors_show = self._get_symbols_colors_from_modes(modes, nopoints=False)
        linestyle = 'None' if self.noline else '-'

        #plot_type = ['tas', 'eas', 'alt', 'kfreq', '1/kfreq', 'freq', 'damp', 'eigr', 'eigi', 'q', 'mach']
        ix, xlabel, xunit = self._plot_type_to_ix_xlabel(plot_type)

        jcolor = 0
        imodes_crossing = []
        xcrossing_dict = {}

        if hasattr(self, 'ieas') and plot_type == 'eas':
            xcrossing_dict = self.get_flutter_crossings(
                damping_required=damping_required, modes=modes,
                eas_range=eas_range)

        if ylim_freq is None or not isinstance(ylim_freq[0], float_types) or not isinstance(ylim_freq[0], float_types):
            filter_freq = False

        if filter_freq:
            yfreq_min, yfreq_max = ylim_freq

        legend_elements = []
        for i, imode, mode in zip(count(), imodes, modes):
            color = colors[jcolor]
            symbol = symbols[jcolor]

            vel = self.results[imode, :, ix].ravel()
            damping = self.results[imode, :, self.idamping].ravel()
            freq = self.results[imode, :, self.ifreq].ravel()

            jcolor, color, linestyle2, symbol2, texti = _increment_jcolor(
                mode, jcolor, color, linestyle, symbol,
                freq, freq_tol=freq_tol,
                show_mode_number=self.show_mode_number)
            if color != 'gray':
                imodes_crossing.append(imode)

            #iplot = np.where(freq > 0.0)
            #damp_axes.plot(vel, damping, symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel, freq, symbols[i])

            #iplot = np.where(freq != np.nan)
            #damp_axes.plot(vel[iplot], damping[iplot], symbols[i], label='Mode %i' % mode)
            #freq_axes.plot(vel[iplot], freq[iplot], symbols[i])
            #print(color, symbol, linestyle)
            label = _get_mode_freq_label(mode, freq[0])
            if filter_freq and freq.min() > ylim_freq[1]:
                # if we're entirely greater than the max, skip line
                continue
            if filter_freq and freq.max() < ylim_freq[0]:
                # if we're entirely below than the min, skip line
                continue

            # _plot_axes(damp_axes,
            #            vel, damping,
            #            color, symbol2, linestyle2,
            #            label, texti)

            #if texti == '':
                #msgi = str([texti, symbol2, linestyle2]) + 'was unexpected...'
                #warnings.warn(msgi)
            legend_elementsi = _plot_two_axes(
                damp_axes, freq_axes,
                vel, damping, freq,
                color, symbol2, linestyle2,
                label, texti,
                self.markevery,
                markersize=None)
            legend_elements.extend(legend_elementsi)
            if ivelocity and symbol2 and ivelocity < len(vel):
                markersize = 10
                plot_kwargs = {
                    'color': 'k', 'marker': 'o',
                    's': markersize**2, 'alpha': 0.8}
                damp_axes.scatter(vel[ivelocity], damping[ivelocity], **plot_kwargs)
                freq_axes.scatter(vel[ivelocity], freq[ivelocity], **plot_kwargs)

        # add horizontal line
        legend_elementsi = _add_damping_limit(
            plot_type, damp_axes, damping_limit)
        legend_elements.extend(legend_elementsi)

        legend_elementsi = _add_vertical_lines(
            [damp_axes, freq_axes], v_lines, plot_type, xunit)
        legend_elements.extend(legend_elementsi)

        # crossings go on top (aka at the end)
        eas_max = None if (xlim is None or xlim[1] is None) else xlim[1]
        legend_elementsi = self._plot_crossings(
            damp_axes,  # damping_required,
            imodes, modes,
            imodes_crossing, xcrossing_dict,
            colors_show, symbols_show,
            eas_max=eas_max,
            filter_damping=filter_damping)
        legend_elements.extend(legend_elementsi)

        damp_axes.set_xlabel(xlabel)  #, size=self.font_size)
        freq_axes.set_xlabel(xlabel)  #, size=self.font_size)
        damp_axes.set_ylabel(r'Structural Damping; $g = 2 \gamma $')  #, fontsize=self.font_size)
        freq_axes.set_ybound(lower=0.)

        damp_axes.grid(True)
        set_xlim(damp_axes, xlim)
        set_ylim(damp_axes, ylim_damping)
        freq_axes.set_ylabel('Frequency [Hz]')  #, fontsize=self.font_size)
        freq_axes.grid(True)

        set_xlim(freq_axes, xlim)
        set_ylim(freq_axes, ylim_freq)

        _set_ticks(self, damp_axes, 0)
        _set_ticks(self, freq_axes, 1)

        title = f'Subcase {self.subcase}'
        if png_filename:
            title += '\n%s' % png_filename

        damp_axes.set_title(title)  #, fontsize=self.font_size)
        #plt.suptitle(title)

        nmodes = len(modes)
        ncol = _update_ncol(nmodes, ncol)
        if legend:
            damp_axes.legend(
                handles=legend_elements,
                bbox_to_anchor=(1.125, 1.), loc=2, ncol=ncol,
                #fontsize=10,
            )
            #damp_axes.legend(fontsize=10, bbox_to_anchor=(1.125, 1.), loc=2, ncol=ncol)
            #fig.subplots_adjust(hspace=0.25)
            #fig.subplots_adjust(hspace=.5)
            #plt.legend()
            #damp_axes.legend(legend_items, bbox_to_anchor=anchor, ncol=2)
            #fig.subplots_adjust(hspace=0.25)
            #fig.subplots_adjust(hspace=.5)

        _show_save_clear_close(
            fig, show, png_filename, clear, close)
        return fig, (damp_axes, freq_axes)

    def _plot_crossings(self,
                        damp_axes: plt.Axes,
                        imodes: np.ndarray,
                        modes: np.ndarray,
                        imodes_crossing: np.ndarray,
                        xcrossing_dict: dict[int, list[Crossing]],
                        colors: list[str],
                        symbols: list[str],
                        eas_max: Optional[float]=None,
                        filter_damping: bool=False) -> list[Line2D]:
        """
        Parameters
        ----------
        filter_damping: bool; default=False
            filter crossings entirely outside the plot range
            useful for cleaning up the legend
        eas_max : float; default=None
            requires filter_damping
            if velocity is greater than the allowable range, don't plot it

        """
        plot_type = 'eas'
        legend_elements = []
        assert isinstance(xcrossing_dict, dict), xcrossing_dict
        if len(xcrossing_dict) == 0:
            return legend_elements
        jcolor = 0
        xunit = self.out_units[plot_type]
        # TODO: fix the colors...
        for i, imode, mode in zip(count(), imodes, modes):
            if imode not in imodes_crossing:
                continue
            if mode not in xcrossing_dict:
                continue
            color = colors[jcolor]
            symbol2 = symbols[jcolor]
            #(freq0, eas0), (freq3, eas3) = xcrossing_dict[mode]
            # freq = self.results[imode, :, self.ifreq].ravel()
            # jcolor, color, linestyle2, symbol2 = _increment_jcolor(
            #     jcolor, color, linestyle, symbol,
            #     freq, freq_tol)
            for case in xcrossing_dict[mode]:
                damping0, freq0, eas0 = case
                if np.isnan(eas0):
                    continue
                if filter_damping and eas0 > eas_max:
                    continue
                damping_str = f'{damping0*100:.0f}%'
                label = f'Mode {mode}: g={damping_str}; {eas0:.0f} {xunit}; {freq0:.1f} Hz'
                legend_element = Line2D([0], [0],
                                        marker=symbol2, color=color, label=label, linestyle='')
                damp_axes.plot(eas0, damping0, color=color,
                               marker=symbol2, linestyle='', label=label)
                legend_elements.append(legend_element)
            jcolor += 1
        return legend_elements

    def export_to_csv(self, csv_filename: PathLike,
                      modes: Optional[list[int]]=None) -> None:
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
        imodes = _imodes(self.results.shape, modes)
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

    def export_to_veas(self, veas_filename: PathLike,
                       modes: Optional[list[int]]=None) -> None:
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
        modes, nmodes = _modes_nmodes(self.results.shape, modes)

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

    # self.page_num = result.write_f06(
    #     f06, header, page_stamp, page_num=self.page_num,
    #     is_mag_phase=is_mag_phase, is_sort1=is_sort1)
    def write_f06(self, f06_file: TextIO,
                  header: list[str],
                  page_stamp: Optional[str]=None,
                  page_num: int=1, **kwargs) -> int:
        page_num = self.export_to_f06_file(
            f06_file, #modes=modes,
            page_stamp=page_stamp, page_num=page_num)
        # return self.export_to_f06(f06_filename,
        #                           page_stamp=page_stamp, page_num=page_num)
        return page_num

    def export_to_f06(self, f06_filename: PathLike,
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
        imodes = _imodes(self.results.shape, modes)
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

    def export_to_zona(self, zona_filename: PathLike,
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
        ix, unused_xlabel, xunit = self._plot_type_to_ix_xlabel(plot_type)

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

    def _plot_type_to_ix_xlabel(self, plot_type: str) -> tuple[int, str, str]:
        """helper method for ``plot_vg_vf``"""
        plot_type = plot_type.lower()
        #print(f'plot_type={plot_type!r} out_units={self.out_units!r}')
        assert isinstance(self.out_units, dict), f'out_units={self.out_units!r}'
        if plot_type == 'tas':
            ix = self.ivelocity
            velocity_units = self.out_units['velocity']
            xlabel = f'Velocity [{velocity_units}]'
            xunit = velocity_units
        elif plot_type == 'eas':
            ix = self.ieas
            velocity_units = self.out_units['eas']
            xlabel = f'Equivalent Airspeed [{velocity_units}]'
            xunit = velocity_units
        elif plot_type == 'alt':
            ix = self.ialt
            alt_units = self.out_units['altitude']
            xlabel = f'Altitude [{alt_units}]'
            xunit = alt_units
        elif plot_type == 'kfreq':
            ix = self.ikfreq
            xlabel = r'Reduced Frequency [rad]; $\omega c / (2V) $'
            xunit = 'rad'
        elif plot_type == 'rho':
            ix = self.idensity
            density_units = self.out_units['density']
            xlabel = f'Density [{density_units}]'
            xunit = density_units
        elif plot_type == 'q':
            ix = self.iq
            pressure_unit = self.out_units['dynamic_pressure']
            xlabel = f'Dynamic Pressure [{pressure_unit}]'
            xunit = pressure_unit
        elif plot_type == 'mach':
            ix = self.imach
            xlabel = 'Mach'
            xunit = ''
        elif plot_type == 'freq':
            ix = self.ifreq
            xlabel = 'Frequency [Hz]'
            xunit = 'Hz'
        elif plot_type in ['1/kfreq', 'ikfreq', 'inv_kfreq', 'kfreq_inv']:
            ix = self.ikfreq_inv
            xlabel = r'1/KFreq [1/rad]; $2V / (\omega c) $'
            xunit = '1/rad'
        elif plot_type == 'eigr':
            ix = self.ieigr
            xlabel = r'Eigenvalue (Real); $\omega \gamma$'
            xunit = 'rad'
        elif plot_type == 'eigi':
            ix = self.ieigi
            xlabel = r'Eigenvalue (Imaginary); $\omega$'
            xunit = 'rad'
        elif plot_type in ['damp', 'damping']:
            ix = self.idamping
            xlabel = r'Structural Damping; $g = 2 \gamma $'
            xunit = 'g'
        else:  # pramga: no cover
            raise NotImplementedError(f"plot_type={plot_type!r} not in ['tas', 'eas', 'alt', 'kfreq', "
                                      "'1/kfreq', 'freq', 'damp', 'eigr', 'eigi', 'q', 'mach', 'alt']")
        return ix, xlabel, xunit

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


def _imodes(results_shape: tuple[int, int],
            modes: Optional[np.ndarray, slice[int] |
                                        tuple[int] | list[int]]) -> np.ndarray:
    """gets the imodes from the modes"""
    if modes is None:
        nmodes = results_shape[0]
        imodes = range(nmodes)
    elif isinstance(modes, slice):
        nmodes = results_shape[0]
        all_modes = np.arange(0, nmodes, dtype='int32')
        imodes = all_modes[modes]
    elif isinstance(modes, (list, tuple)):
        modes = np.array(modes, dtype='int32')
        imodes = modes - 1
    else:
        imodes = modes - 1
    return imodes


def _modes_nmodes(results_shape: tuple[int, int],
                  modes: Optional[Iterable[int]]) -> tuple[Iterable[int], int]:
    """gets the modes and nmodes"""
    if modes is None:
        nmodes = results_shape[0]
        modes = range(1, nmodes + 1)
    else:
        nmodes = max(modes)
    return modes, nmodes


def _get_modes_imodes(all_modes: np.ndaray,
                      modes: Optional[slice | np.ndarray]):
    """gets the index of the modes to plot"""
    if modes is None:
        modes = all_modes
    elif isinstance(modes, slice):
        start = modes.start
        #if modes.stop is None:
            #stop = len(all_modes) + 1
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


def get_zero_crossings(x: np.ndarray,
                       freq: np.ndarray,
                       y: np.ndarray) -> np.ndarray:
    """https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python

    Haven't tested this yet...
    """
    # zero_crossings = np.where(np.diff(np.sign(y)))[0]
    zero_crossings = np.where(np.diff(np.sign(y) >= 0))[0]
    # zero_crossings = np.where(np.diff(np.signbit(y)))[0]
    #assert len(zero_crossings) == 1, zero_crossings
    if len(zero_crossings) == 0:
        return np.nan, np.nan

    iy = zero_crossings[0]
    yii = y[iy]
    if yii < 0:
        x0 = x[iy]
        x1 = x[iy+1]
        y0 = y[iy]
        y1 = y[iy + 1]
        #y2 = (y1-y0) / (x1-x0) * (x-x0) + y0
        x2 = (x1-x0) / (y1-y0) * (0.-y0) + x0

        f0 = freq[iy]
        f1 = freq[iy+1]
        f2 = (f1-f0) / (y1-y0) * (0.-y0) + f0
    else:
        return np.nan, np.nan
    return x2, f2


def _asarray(results, allow_fix_kfreq: bool=True):
    """casts the results array"""
    allow_fix_kfreq = True
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
                       linewidth: int=2) -> list[Line2D]:
    if damping_limit is None:
        return []
    #damp_label = f'Damping={damping_limit*100:.1f}'
    #plt.axhline(y=1.0, color="black", linestyle="--")
    line1 = damp_axes.axhline(
        y=0., color='k', linestyle='--', linewidth=linewidth,
        label=f'Structural Damping=0%')
    line2 = damp_axes.axhline(
        y=damping_limit, color='k', linestyle='-', linewidth=linewidth,
        label=f'Limit Structural Damping={damping_limit*100:.0f}%')

    legend_elements = [line1, line2]
    return legend_elements


def get_flutter_units(units: Optional[str | dict[str, str]]) -> dict[str, str]:
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
            units_dict = {
                'velocity': 'm/s', 'density': 'kg/m^3',
                'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'm/s'}
        elif units == 'si_mm':
            units_dict = {
                'velocity': 'mm/s', 'density': 'Mg/mm^3',
                'altitude': 'm', 'dynamic_pressure': 'MPa', 'eas': 'mm/s'}
        #elif units == 'si_cmgs':
            #units = {'velocity': 'cm/s', 'density': 'g/cm^3',
                     #'altitude': 'm', 'dynamic_pressure': 'Pa', 'eas': 'cm/s'}
        elif units == 'english_in':
            units_dict = {
                'velocity': 'in/s', 'density': 'slinch/in^3',
                'altitude': 'ft', 'dynamic_pressure': 'psi', 'eas': 'in/s'}
        elif units == 'english_ft':
            units_dict = {
                'velocity': 'ft/s', 'density': 'slug/ft^3',
                'altitude': 'ft', 'dynamic_pressure': 'psf', 'eas': 'ft/s'}
        elif units == 'english_kt':
            units_dict = {
                'velocity': 'knots', 'density': 'slug/ft^3',
                'altitude': 'ft', 'dynamic_pressure': 'psf', 'eas': 'knots'}
        else:  # pragma: no cover
            raise NotImplementedError(f'units={units!r} must be in [si, si_mm, '
                                      'english_in, english_ft, english_kt]')
    else:  # pragma: no cover
        assert isinstance(units, dict), f'units={units!r}'
        required_keys = ['altitude', 'velocity', 'eas', 'density', 'dynamic_pressure']
        for key in required_keys:
            assert key in units, 'key=%r not in units=%s' % (key, units)
        units_dict = units
    return units_dict

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


def _get_mode_freq_label(mode: int, freq: float) -> str:
    # write tiny numbers
    if abs(freq) > 1.0:
        # don't write big numbers in scientific
        freq_num = f'{freq:.1f}'
    else:
        freq_num = f'{freq:.3g}'
        # strip silly scientific notation
        freq_num = freq_num.replace('-0', '-').replace('-0', '-').replace('+0', '+')
    label = f'Mode {mode:d}; freq={freq_num} Hz'
    return label


def _increment_jcolor(mode: int,
                      jcolor: int, color: str,
                      linestyle: str, symbol: str,
                      freq: np.ndarray,
                      freq_tol: float=-1.0,
                      show_mode_number: bool=False,
                      # jcolor, color, linestyle2, symbol2, text
                      ) -> tuple[int, Color, str, str, str]:
    """
    Filters a line if it doesn't change by more than freq_tol.
    Changes the line color and removes the symbol.

    Parameters
    ----------
    linestyle: str
        '-', '--', 'None'
    freq: np.ndarray
        the frequency data
    freq_tol: float; default=-1.0
        -1.0: no filtering (default)
        >0.0: filter is active

    Returns
    -------
    linestyle2: str
        the updated style

    """
    #print(f'freq = {freq}')
    #print(f'freq_tol = {freq_tol!r}')
    assert isinstance(freq, np.ndarray), freq
    assert isinstance(freq_tol, float_types), freq_tol
    is_filtered = False
    if freq.max() - freq.min() <= freq_tol:
        color = 'gray'
        is_filtered = True
        jcolor -= 1
    linestyle2 = '--' if is_filtered else linestyle
    symbol2 = '' if is_filtered else symbol

    text = ''
    if show_mode_number and symbol2:
        symbol2 = ''
        text = str(mode)
    jcolor += 1
    return jcolor, color, linestyle2, symbol2, text

def set_xlim(axes: plt.Axes, xlim: Limit) -> None:
    if xlim == [None, None] or xlim == (None, None):
        xlim = None
    if xlim is not None:
        axes.set_xlim(xlim)

def set_ylim(axes: plt.Axes, ylim: Limit) -> None:
    if ylim == [None, None] or ylim == (None, None):
        ylim = None
    if ylim is not None:
        axes.set_ylim(ylim)


def _get_unit_factor(in_units: dict[str, str],
                     out_units: dict[str, str],
                     name: str) -> tuple[float, str]:
    if not in_units or not out_units:
        msg = f'name={name!r} in_units={in_units} out_units={out_units}'
        raise RuntimeError(msg)
    unit_f06 = in_units[name]
    unit_out = out_units[name]

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

    if out_units is not None:
        units = out_units[name]
    else:
        units = 'units'
    return factor, units

def get_legend_kwargs(font_size: int,
                      legend_kwargs: Optional[dict[str, Any]],
                      ) -> dict[str, Any]:
    if legend_kwargs is None:
        legend_kwargs = {}
    assert isinstance(legend_kwargs, dict), legend_kwargs

    legend_kwargs_check = ['loc', 'fancybox', 'framealpha']
    for key in legend_kwargs:
        assert key in legend_kwargs_check, key

    #if 'prop' not in legend_kwargs:
        #legend_kwargs['prop'] = {'size': font_size}
    #if 'fontsize' not in legend_kwargs:
        #legend_kwargs['fontsize'] = font_size

    return legend_kwargs

def _is_q_units_consistent(rho_units: str, vel_units: str,
                           q_units: str) -> bool:
    units = [
        ('kg/m^3', 'm/s', 'Pa'),
        ('Mg/mm^3', 'mm/s', 'MPa'),
        ('slug/ft^3', 'ft/s', 'psf'),
        ('slinch/in^3', 'in/s', 'psi'),
    ]
    is_consistent = (rho_units, vel_units, q_units) in units
    return is_consistent


def check_range(eas_min0: float, eas_max0: float,
                freq: float, eas: float) -> tuple[float, float]:
    """
    Flutter is only valid if it's between eas_min0 and eas_max0.  In
    other words, why does it matter if flutter is outside the flight
    envelope?
    """
    if np.isnan(eas):
        return eas, freq
    if eas_min0 is not None and eas < eas_min0:
        eas = np.nan
        freq = np.nan
        return eas, freq

    if eas_max0 is not None and eas > eas_max0:
        eas = np.nan
        freq = np.nan
    return eas, freq


def _get_min_damping(damping_required: list[tuple[float, float]]) -> float:
    """
    Get min required damping for 0% and 3%, such that flutter exists.
    The goal of this is to filter all points that are below 4% to
    limit identifying an excessive number of flutter crossings.

    damping_required = [
        (0.00, 0.01),
        (0.03, 0.04),
    ]
    damping = _get_min_damping(damping_required)
    >>> 0.04
    """
    min_damping = damping_required[0][1]
    for dampingi, requiredi in damping_required[1:]:
        min_damping = min(min_damping, requiredi)
    return min_damping


def _plot_two_axes(damp_axes: plt.Axes, freq_axes: plt.Axes,
                   vel: np.ndarray, damping: np.ndarray, freq: np.ndarray,
                   color: str, symbol: str, linestyle: str,
                   label: str, text: str,
                   markevery: Optional[int]=None,
                   markersize=None) -> list[Line2D]:
    legend_element = Line2D([0], [0], color=color,
                            marker=symbol, label=label, linestyle=linestyle)
    legend_elements = [legend_element]

    # setup for plotting every Nth point
    markevery2 = None if markevery is None else markevery + 1
    vel2 = vel[::markevery2]
    damping2 = damping[::markevery2]
    freq2 = freq[::markevery2]

    if markevery2 is None:
        # plot all points and lines (default)
        line = damp_axes.plot(vel, damping, color=color, marker=symbol, markersize=markersize, linestyle=linestyle, label=label)
        freq_axes.plot(vel, freq, color=color, marker=symbol, markersize=markersize, linestyle=linestyle)
    elif symbol or text or linestyle:
        # draw lines with all points
        line = damp_axes.plot(vel, damping, color=color, linestyle=linestyle, label=label)
        freq_axes.plot(vel, freq, color=color, linestyle=linestyle)

        if symbol:
            # plot every other point (for reduced clutter)
            damp_axes.scatter(vel2, damping2, color=color, marker=symbol, s=markersize)
            freq_axes.scatter(vel2, freq2, color=color, marker=symbol, s=markersize)
    if text:
        # annotate the mode number
        for xi, y1i, y2i in zip(vel2, damping2, freq2):
            damp_axes.text(xi, y1i, text, color=color, clip_on=True)
            freq_axes.text(xi, y2i, text, color=color, clip_on=True)
    return legend_elements


def _show_save_clear_close(fig: plt.Figure,
                           show: bool,
                           png_filename: Optional[str],
                           clear: bool,
                           close: bool) -> None:
    #print('_show_save_clear_close')
    if show:
        plt.show()
    if png_filename:
        plt.savefig(png_filename)
    if clear:
        fig.clear()
    if close:
        plt.close()


def reshape_eigenvectors(eigenvectors: np.array,
                         eigr_eigi_vel: np.array,
                         incorrect_shape: bool=False) -> np.ndarray:
    """
    Parameters
    ----------
    incorrect_shape: bool
        helper for testing
    """
    nmodes1, nmodes_nvel = eigenvectors.shape
    nmodes = nmodes1 - 1 if incorrect_shape else nmodes1
    nvel = nmodes_nvel // nmodes
    if 0:  # pragma: no cover
        print(nmodes1, nmodes)
        print(nmodes_nvel, nmodes, nvel)
        print('eigenvectors:')
        print(eigenvectors)
        print('eigr_eigi_vel:')
        print(eigr_eigi_vel)
        print(f'nmodes={nmodes}; nvel={nvel}')
    assert nvel > 0, nvel

    i = 0
    # eigenvectors2 = np.zeros((nmodes, nmodes, nvel), dtype=eigenvectors.dtype)
    # for ivel in range(nvel):
    #     for imode in range(nmodes):
    #         eigenvectors2[:, imode, ivel] = eigenvectors[:, i]
    #         i += 1

    # was 1,2
    eigenvectors3 = eigenvectors.reshape(nmodes1, nvel, nmodes).swapaxes(0, 1).swapaxes(1, 2)
    eigr_eigi_vel3 = eigr_eigi_vel.reshape((nvel, nmodes, 3))
    # assert eigenvectors2.shape == eigenvectors3.shape, (eigenvectors2.shape, eigenvectors3.shape)
    #assert np.allclose(eigenvectors2, eigenvectors3)
    # print(data3[:, :, 0])
    # print(data3[:, :, 1])
    # print(data3[:, :, 2])
    # print(data3[:, :, 3])
    #for ivel in range(nvel):
        #print(f'ivel={ivel}')
        #print(eigenvectors2[:, :, ivel])
        #assert np.allclose(eigenvectors2[:, :, ivel], eigenvectors3[:, :, ivel])

    # we want the rows
    #ivel = 0
    #imode = 1
    #print(eigenvectors.shape)
    #mpf = eigenvectors[:, imode, ivel]
    #eig.scale(mpf)
    #asdf
    return eigenvectors3, eigr_eigi_vel3

def _add_vertical_lines(axes_list: list[Axes],
                        v_lines: Optional[list[LineData]],
                        plot_type: str, xunit: str,
                        linewidth: int=2) -> list[Line2D]:
    """the first plot gets the label"""
    legend_elements = []
    if v_lines is None or plot_type not in {'tas', 'eas'}:
        return legend_elements

    for v_line in v_lines:
        name, velocity, vcolor, linestyle = v_line
        assert linestyle in LINESTYLES, (name, linestyle)
        for iaxis, axes in enumerate(axes_list):

            # put the label on the first plot
            if iaxis == 0:
                if velocity == int(velocity):
                    label = f'{name}={velocity:.0f} [{xunit}]'
                else:
                    label = f'{name}={velocity:.1f} [{xunit}]'
                legend_element = axes.axvline(
                    x=velocity, color=vcolor, linestyle=linestyle,
                    linewidth=linewidth, label=label)
                legend_elements.append(legend_element)
            else:
                axes.axvline(
                    x=velocity, color=vcolor, linestyle=linestyle,
                    linewidth=linewidth)
            legend_elements.append(legend_element)
    return legend_elements


def _set_ticks(self, axes: plt.Axes, iaxis: int) -> None:
    axes.tick_params(axis='both', which='major')  # , labelsize=self.font_size)
    if _is_tick(self._xtick_major_locator_multiple, iaxis):
        axes.xaxis.set_major_locator(MultipleLocator(self._xtick_major_locator_multiple[iaxis]))
    if _is_tick(self._ytick_major_locator_multiple, iaxis):
        axes.yaxis.set_major_locator(MultipleLocator(self._ytick_major_locator_multiple[iaxis]))

def _is_tick(values: Optional[tuple[float, ...]], index: int):
    out = values is not None and values[index] is not None
    return out
