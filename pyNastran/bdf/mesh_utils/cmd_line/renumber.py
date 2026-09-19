from __future__ import annotations
import os
import sys
import argparse
from cpylog import SimpleLogger

import pyNastran


def cmd_line_renumber(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_renumber"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    if len(argv) == 1:
        sys.exit("bdf renumber: use 'bdf renumber -h' for help")

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf renumber',
        description='Renumber a BDF model',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('OUT_BDF_FILENAME', nargs='?', default=None,
                        help='path to output BDF/DAT/NAS file')
    parser.add_argument('--nid', type=int, default=None, metavar='NID',
                        help='starting node id')
    parser.add_argument('--eid', type=int, default=None, metavar='EID',
                        help='starting element id')
    parser.add_argument('--pid', type=int, default=None, metavar='PID',
                        help='starting property id')
    parser.add_argument('--mid', type=int, default=None, metavar='MID',
                        help='starting material id')
    parser.add_argument('--superelement', action='store_true',
                        help='calls superelement_renumber')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('-x', '--xref', action='store_true',
                        help='flag to disable cross-referencing')
    parser.add_argument('--size', type=int, default=16,
                        help='set the field size (default=16)')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    bdf_filename_out = args.OUT_BDF_FILENAME
    if bdf_filename_out is None:
        base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f'{base}.renumber{ext}'

    size = args.size
    xref = not args.xref
    assert size in [8, 16], f'size={size} args={argv}'
    cards_to_skip = []

    starting_id_dict = {}
    for arg in {'nid', 'eid', 'pid', 'mid'}:
        value = getattr(args, arg)
        if value is not None:
            starting_id_dict[arg] = value
    if len(starting_id_dict) == 0:
        starting_id_dict = None
    else:
        log.debug(f'starting_id_dict = {starting_id_dict}')

    from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber, superelement_renumber, _get_bdf_model

    model = _get_bdf_model(
        bdf_filename, punch=punch, xref=xref,
        cards_to_skip=cards_to_skip, log=log, debug=True)
    if args.superelement:
        superelement_renumber(model, bdf_filename_out, size=size, is_double=False,
                              starting_id_dict=starting_id_dict,  #round_ids=False,
                              log=log)
    else:
        bdf_renumber(model, bdf_filename_out, size=size, is_double=False,
                     starting_id_dict=starting_id_dict, round_ids=False,
                     log=log)
