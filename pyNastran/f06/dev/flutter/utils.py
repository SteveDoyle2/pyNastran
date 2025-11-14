from __future__ import annotations
import os
from pathlib import Path
import shutil
from typing import Optional, TYPE_CHECKING
from matplotlib import pyplot as plt
from pyNastran.f06.parse_flutter import make_flutter_response, get_flutter_units

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.f06.flutter_response import FlutterResponse
    from pyNastran.op2.op2 import OP2


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
