from __future__ import annotations
import os
from pathlib import Path
import shutil
import traceback
from typing import Optional, TYPE_CHECKING
from matplotlib import pyplot as plt
from pyNastran.f06.parse_flutter import make_flutter_response, get_flutter_units

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.f06.flutter_response import FlutterResponse
    from pyNastran.op2.op2 import OP2


X_PLOT_TYPES = ['eas', 'tas', 'rho', 'q', 'mach', 'alt', 'kfreq', 'ikfreq', 'index']
PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus', 'modal-participation',
              'zimmerman']
UNITS_IN = ['english_in', 'english_kt', 'english_ft',
            'si', 'si_mm']
MODE_SWITCH_METHODS = ['None', 'Frequency', 'Damping']

def load_f06_op2(f06_filename: str, log: SimpleLogger,
                 in_units: str,
                 out_units: str,
                 use_rhoref: bool,
                 make_alt: bool=False) -> tuple[OP2, dict[int, FlutterResponse]]:
    """
    load a Vg-Vf plot from:
     - OP2 / F06
     - dict[subcase: int, FlutterResponse]

    From an OP2, load:
     - eigenvectors, vg_vf_response
    From an F06, load:
     - vg_vf_response
    """
    model = None
    responses = {}
    if not os.path.exists(f06_filename):
        log.error(f'Cant find {f06_filename}')
        return model, responses

    in_units_dict = get_flutter_units(in_units)
    out_units_dict = get_flutter_units(out_units)
    ext = os.path.splitext(f06_filename)[1].lower()
    #print(f'use_rhoref={use_rhoref}')
    if ext == '.f06':
        try:
            responses, mass = make_flutter_response(
                f06_filename,
                f06_units=in_units_dict,
                out_units=out_units_dict,
                use_rhoref=use_rhoref,
                make_alt=make_alt,
                log=log)
        except Exception as e:
            log.error(str(e))
            #raise
            return model, responses
    elif ext == '.out':
        from pyNastran.f06.dev.flutter.read_zona_out import read_zona_out
        try:
            responses, mass = read_zona_out(f06_filename)
        except Exception as e:
            log.error(str(e))
            #raise
            return model, responses
    elif ext == '.op2':
        try:
            from pyNastran.op2.op2 import OP2
        except ImportError as e:
            log.error(str(e))
            return model, responses

        assert isinstance(in_units_dict, dict), in_units_dict
        model = OP2(log=log)
        model.in_units = in_units_dict
        results_to_include = ['eigenvectors', 'vg_vf_response']
        model.set_results(results_to_include)
        try:
            model.read_op2(f06_filename, build_dataframe=False)
            responses = model.op2_results.vg_vf_response
            if len(responses) == 0:
                log.error('Could not find OVG table in op2')
        except Exception as e:
            log.error(str(e))
            return model, responses
    else:
        log.error(f'Invalid file type; {f06_filename}')
        return model, responses

    for response in responses.values():
        response.convert_units(out_units_dict)
    return model, responses



def get_png_filename(base: str, x_plot_type: str, plot_type: str,
                     export_to_png: bool) -> tuple[str, Optional[str]]:
    assert isinstance(base, str), base
    assert isinstance(x_plot_type, str), x_plot_type
    assert isinstance(plot_type, str), plot_type
    assert isinstance(export_to_png, bool), export_to_png
    if 'x-' in plot_type:
        # png_filename0 = base + f'_{x_plot_type}-damp-kfreq.png'
        # png_filename0 = base + f'_{x_plot_type}-damp-freq.png'
        plot_type2 = plot_type.replace('x-', x_plot_type + '-')
        png_filename0 = base + f'_{plot_type2}.png'
    else:
        #png_filename0 = base + '_root-locus.png'
        #png_filename0 = base + '_modal-participation.png'
        png_filename0 = base + f'_{plot_type}.png'
    #print(f'png_filename0 = {png_filename0}')
    png_filename = png_filename0 if export_to_png else None
    #print(f'png_filename = {png_filename}')
    return png_filename0, png_filename


def update_ylog_style(fig: plt.Figure,
                      log_scale_x: bool,
                      log_scale_y1: bool,
                      log_scale_y2: bool) -> None:
    ax_list = fig.axes
    xscale = 'log' if log_scale_x else 'linear'
    yscale1 = 'log' if log_scale_y1 else 'linear'
    yscale2 = 'log' if log_scale_y2 else 'linear'
    ax_list[0].set_xscale(xscale)
    ax_list[1].set_xscale(xscale)

    ax_list[0].set_yscale(yscale1)
    ax_list[1].set_yscale(yscale2)

def get_plot_file() -> str:
    home_dirname = Path(os.path.expanduser('~'))
    old_filename = home_dirname / 'plot_flutter.json'
    new_filename = home_dirname / 'plot_145' / 'settings.json'
    #move_filename(old_filename, new_filename)
    return str(old_filename)

def move_filename(old_filename: Path,
                  new_filename: Path) -> None:
    """
    Creates new file from old file is new file doesn't exist.
    Cleans up the mess.
    """
    dirname = os.path.dirname(new_filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    if old_filename.exists() and not new_filename.exists():
        shutil.copyfile(old_filename, new_filename)
    if old_filename.exists():
        os.remove(old_filename)
    return


def get_point_removal_str(point_removal: list[tuple[float, float]]):
    """
    >>> point_removal = [[400.0, 410.0], [450.0, 500.0]]
    point_removal_str = get_point_removal_str(point_removal)
    >>> point_removal_str
    '400:410,450:500'

    >>> point_removal = [[450.0, -1.0]]
    point_removal_str = get_point_removal_str(point_removal)
    >>> point_removal_str
    '450:'

    >>> point_removal = [[-1.0, 500.0]]
    point_removal_str = get_point_removal_str(point_removal)
    >>> point_removal_str
    ':500'
    """
    out = []
    for mini, maxi in point_removal:
        if mini > 0 and maxi > 0:
            outi = f'{mini:g}:{maxi:g}'
        elif mini > 0:
            outi = f'{mini:g}:'
        elif maxi > 0:
            outi = f':{maxi:g}'
        else:
            continue
        out.append(outi)
    point_removal_str = ','.join(out)
    return point_removal_str


def point_removal_str_to_point_removal(point_removal_str: str,
                                       log: SimpleLogger) -> list[tuple[float, float]]:
    point_removal = []
    point_removal_list = point_removal_str.split(',')

    if point_removal_list == ['']:
        return point_removal_list

    try:
        for ipoint, point in enumerate(point_removal_list):
            sline = point.split(':')
            assert len(sline) == 2, f'point_removal[{ipoint}]={sline}; point_removal={str(point_removal_list)}'
            a_str = sline[0].strip()
            b_str = sline[1].strip()
            a = float(a_str) if a_str != '' else -1.0
            b = float(b_str) if b_str != '' else -1.0
            point_float = (a, b)
            point_removal.append(point_float)
    except Exception as e:
        log.error(str(e))
        # print(traceback.print_tb(e))
        print(traceback.print_exception(e))
    return point_removal


def _to_str(value: Optional[int | float]) -> str:
    if value is None:
        str_value = ''
    else:
        str_value = str(value)
    return str_value


def _float_passed_to_default(value: float, is_passed: bool,
                             default: float=-1.0) -> float:
    if is_passed and value is None:
        value = default
    return value


def get_plot_flags(plot_type: str,
                   x_plot_type: str) -> dict[str, bool]:
    show_index_lim = False
    show_eas_lim = False
    show_tas_lim = False
    show_mach_lim = False
    show_alt_lim = False
    show_q_lim = False
    show_rho_lim = False

    show_xlim = False
    show_freq = False
    show_damp = False
    show_root_locus = False
    show_zimmerman = False
    show_modal_participation = False

    # PLOT_TYPES = ['x-damp-freq', 'x-damp-kfreq', 'root-locus']
    assert plot_type in PLOT_TYPES, plot_type

    if x_plot_type == 'kfreq':
        show_kfreq = True
    else:
        show_kfreq = False

    if x_plot_type == 'ikfreq':
        show_ikfreq = True
    else:
        show_ikfreq = False

    if plot_type == 'x-damp-freq':
        show_xlim = True
        show_damp = True
        show_freq = True
    elif plot_type == 'x-damp-kfreq':
        # kfreq-damp-kfreq not handled
        show_xlim = True
        show_damp = True
        show_kfreq = True
    elif plot_type == 'zimmerman':
        show_zimmerman = True
    elif plot_type == 'root-locus':
        show_root_locus = True
        # show_kfreq = False
    elif plot_type == 'modal-participation':
        show_modal_participation = True
        # show_kfreq = False
    else:  # pragma: no cover
        raise RuntimeError(f'plot_type={plot_type!r}')

    if show_xlim:
        if 'index' == x_plot_type:
            show_index_lim = True
        elif 'eas' == x_plot_type:
            show_eas_lim = True
        elif 'tas' == x_plot_type:
            show_tas_lim = True
        elif 'mach' == x_plot_type:
            show_mach_lim = True
        elif 'alt' == x_plot_type:
            show_alt_lim = True
        elif 'q' == x_plot_type:
            show_q_lim = True
        elif 'rho' == x_plot_type:
            show_rho_lim = True
        elif 'kfreq' == x_plot_type:
            show_kfreq_lim = True
        elif 'ikfreq' == x_plot_type:
            show_ikfreq_lim = True
    flags = {
        'show_index_lim': show_index_lim,
        'show_eas_lim': show_eas_lim,
        'show_tas_lim': show_tas_lim,
        'show_mach_lim': show_mach_lim,
        'show_alt_lim': show_alt_lim,
        'show_q_lim': show_q_lim,
        'show_rho_lim': show_rho_lim,

        'show_xlim': show_xlim,
        'show_freq': show_freq,
        'show_damp': show_damp,
        'show_kfreq': show_kfreq,
        'show_ikfreq': show_ikfreq,
        'show_root_locus': show_root_locus,
        'show_modal_participation': show_modal_participation,
        'show_zimmerman': show_zimmerman,
    }
    return flags
