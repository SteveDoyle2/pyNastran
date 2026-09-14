from __future__ import annotations
import sys
import argparse
from typing import TYPE_CHECKING
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args, get_bdf_outfilename
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

    filter_no_args("bdf delete_bad_shells: use 'bdf delete_bad_shells -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf delete_bad_shells',
        description='Delete bad shell elements based on quality metrics',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--skew', type=float, default=70.,
                        help='the maximum skew angle (default=70.0)')
    parser.add_argument('--max_theta', type=float, default=175.,
                        help='the maximum interior angle (default=175.0)')
    parser.add_argument('--min_theta', type=float, default=0.1,
                        help='the minimum interior angle (default=0.1)')
    parser.add_argument('--max_ar', type=float, default=100.,
                        help='the maximum aspect ratio (default=100.0)')
    parser.add_argument('--max_taper', type=float, default=4.,
                        help='the maximum taper ratio (default=4.0)')
    parser.add_argument('--max_warp', type=float, default=90.,
                        help='the maximum warp angle (default=90.0)')

    args = parser.parse_args(argv[2:])

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    bdf_filename_out = get_bdf_outfilename(
        bdf_filename, bdf_filename_out=None,
        tag='fixedquality')

    skew = args.skew
    max_theta = args.max_theta
    min_theta = args.min_theta
    max_aspect_ratio = args.max_ar
    max_taper_ratio = args.max_taper
    max_warping = args.max_warp

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

    if not quiet:  # pragma: no cover
        print(vars(args))

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


