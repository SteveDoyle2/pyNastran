from __future__ import annotations
import sys
from typing import Any, TYPE_CHECKING
from docopt import docopt

import pyNastran
from .utils import (
    get_bdf_filename_punch_log, filter_no_args,
    get_bdf_outfilename)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


SHELL_QUALITY = (
    '[--skew SKEW] [--max_theta MAX_THETA] [--min_theta MIN_THETA] '
    '[--max_ar MAX_AR] [--max_taper MAX_TAPER] [--max_warp MAX_WARP]'
)
def cmd_line_delete_bad_shells(argv=None, quiet: bool=False) -> None:
    """command line interface to ``delete_bad_shells``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        f'  bdf delete_bad_shells IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] {SHELL_QUALITY}\n'
        '  bdf delete_bad_shells -h | --help\n'
        '  bdf delete_bad_shells -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME        path to input BDF/DAT/NAS file\n"
        #"  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '  --punch                            flag to identify a *.pch/*.inc file\n'
        "  --skew SKEW            The maximum skew angle (default=70.0)\n"
        "  --max_theta MAX_THETA  The maximum interior angle (default=175.0)\n"
        "  --min_theta MIN_THETA  The minimum interior angle (default=0.1)\n"
        "  --max_ar MAX_AR        The maximum aspect ratio (default=100.0)\n"
        "  --max_taper MAX_TAPER  The maximum taper ratio (default=4.0)\n"
        "  --max_warp MAX_WARP    The maximum warp angle (default=90.0)\n\n"

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)

    bdf_filename_out = get_bdf_outfilename(
        bdf_filename, bdf_filename_out=None,
        tag='fixedquality')

    #TOLERANCE LIMITS ARE:
    #   SA = 30.00
    #   IA(MIN) = 30.00
    #   IA(MAX) = 150.00
    #   WF = 0.05
    #   TR = 0.50
    #   AR = 100.00
    #
    # Skew angle for the quadrilateral element is defined to be the angle between the lines that join
    # midpoints of the opposite sides of the quadrilateral. Skew angle for the triangular element is
    # defined to be the smallest angle at any of the three vertices.
    #
    # Taper ratio for the quadrilateral element is defined to be the absolute value of the ratio of the area
    # of the triangle formed at each corner grid point to one half the area of the quadrilateral minus
    # 1.0. The largest of the four ratios is compared against the tolerance value. Note that as the ratio
    # approaches 0.0, the shape approaches a rectangle.
    # taper = |atri / (0.5 * aquad) - 1 |
    #
    # Surface warping factor for a quadrilateral is defined to be the distance of the corner points of the
    # element to the mean plane of the grid points divided by the average of the element diagonal
    # lengths. For flat elements (such that all the grid points lie in a plane), this factor is zero.

    defaults = {
        '--skew': 70.,
        '--max_theta': 175.,
        '--min_theta': 0.1,
        '--max_ar': 100.,
        '--max_taper': 4.,
        '--max_warp': 90.,

        #'SKEW': 70.,
        #'MAX_THETA': 175.,
        #'MIN_THETA': 0.1,
        #'MAX_AR': 100.,
        #'MAX_TAPER': 4.,
        #'MAX_WARP': 90.,
    }
    _apply_float_values_to_dict(data, defaults)
    try:
        skew = float(data['--skew'])
        max_theta = float(data['--max_theta'])
        min_theta = float(data['--min_theta'])
        max_aspect_ratio = float(data['--max_ar'])
        max_taper_ratio = float(data['--max_taper'])
        max_warping = float(data['--max_warp'])
        #skew = float(data['SKEW'])
        #max_theta = float(data['MAX_THETA'])
        #min_theta = float(data['MIN_THETA'])
        #max_aspect_ratio = float(data['MAX_AR'])
        #max_taper_ratio = float(data['MAX_TAPER'])
        #max_warping = float(data['MAX_WARP'])
    except:
        if not quiet:  # pragma: no cover
            print(data)
        raise

    if not quiet:  # pragma: no cover
        print(data)
    #print('max_aspect_ratio =', max_aspect_ratio)
    #sss
    #if bdf_filename_out is None:
        #bdf_filename_out = 'merged.bdf'
    size = 8
    from pyNastran.bdf.bdf import read_bdf, BDF
    from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells

    model: BDF = read_bdf(bdf_filename, validate=True, xref=True, punch=punch,
                          log=log, debug=True, mode='msc')  # encoding=None,
    delete_bad_shells(model,
                      min_theta=min_theta, max_theta=max_theta, max_skew=skew,
                      max_aspect_ratio=max_aspect_ratio, max_taper_ratio=max_taper_ratio,
                      max_warping=max_warping)
    model.write_bdf(bdf_filename_out, size=size,
                    nodes_size=16, elements_size=16, loads_size=8)


def _apply_float_values_to_dict(data: dict[str, Any],
                                defaults: dict[str, float]) -> None:
    for name, default_value in defaults.items():
        if data[name] is None:
            #print(f'applying {name}')
            data[name] = default_value
