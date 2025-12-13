from __future__ import annotations
import os
import sys
from docopt import docopt

import pyNastran
from .utils import get_bdf_filename_punch_log


def cmd_line_renumber(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_renumber"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    options = '[--nid NID] [--eid EID] [--pid PID] [--mid MID] [--punch]'
    # TODO: add punch?
    msg = (
        "Usage:\n"
        f'  bdf renumber IN_BDF_FILENAME OUT_BDF_FILENAME [--superelement] [--size SIZE] {options}\n'
        f'  bdf renumber IN_BDF_FILENAME                  [--superelement] [--size SIZE] {options}\n'
        '  bdf renumber -h | --help\n'
        '  bdf renumber -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '--nid NID       starting node id\n'
        '--eid EID       starting element id\n'
        '--pid PID       starting property id\n'
        '--mid MID       starting material id\n'
        '--superelement  calls superelement_renumber\n'
        '--punch         flag to identify a *.pch/*.inc file\n'
        '--size SIZE     set the field size (default=16)\n\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)
    bdf_filename_out = data['OUT_BDF_FILENAME']
    if bdf_filename_out is None:
        base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f'{base}.renumber{ext}'

    size = 16
    if data['--size']:
        if 'SIZE' in data:
            size_str = data['SIZE']
        else:
            size_str = data['--size']
        size = int(size_str)

    assert size in [8, 16], f'size={size} args={argv}'
    #punch = data['--punch']
    # cards_to_skip = [
    #     'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
    #     'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []

    starting_id_dict = {}
    for arg in {'nid', 'eid', 'pid', 'mid'}:
        dash_arg = f'--{arg}'
        if dash_arg in data and data[dash_arg] is not None:
            starting_id_dict[arg] = int(data[dash_arg])
    if len(starting_id_dict) == 0:
        starting_id_dict = None
    else:
        log.debug(f'starting_id_dict = {starting_id_dict}')

    from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber, superelement_renumber
    if data['--superelement']:
        superelement_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                              starting_id_dict=starting_id_dict,  #round_ids=False,
                              cards_to_skip=cards_to_skip, log=log)
    else:
        bdf_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False, punch=punch,
                     starting_id_dict=starting_id_dict, round_ids=False,
                     cards_to_skip=cards_to_skip, log=log)
