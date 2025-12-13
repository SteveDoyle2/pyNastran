from __future__ import annotations
import os
import sys
from typing import TYPE_CHECKING
from cpylog import SimpleLogger
import numpy as np

import pyNastran
from .utils import add_argparse_arguments
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_solid_dof(argv=None, quiet: bool=False,
                       ) -> tuple[BDF, np.ndarray]:
    """command line interface to solid_dof"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    from pyNastran.utils.arg_handling import argparse_to_dict, update_message
    parent_parser = argparse.ArgumentParser()
    # positional arguments
    parent_parser.add_argument('IN_BDF_FILENAME', help='path to input BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('-o', '--output', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('--spc', default=100, help='SPC ID (default=100)', type=int)
    add_argparse_arguments(parent_parser, ['--punch', '--lax', '--allow_dup'])
    parent_parser.add_argument('-v', '--version', action='version', version=pyNastran.__version__)
    args = parent_parser.parse_args(args=argv[2:])
    data = argparse_to_dict(args)

    if not quiet:  # pragma: no cover
        print(data)
    #size = 16

    bdf_filename = data['IN_BDF_FILENAME']
    out_filename = data['output']
    spc_id = data['spc']
    punch = args.punch
    is_strict_card_parser = not args.lax
    if out_filename is None:
        base, ext = os.path.splitext(bdf_filename)
        out_filename = base + '.solid_dof_constraint.blk'

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    is_strict_card_parser = not args.lax
    from .utils_bdf import read_lax_bdf
    model = read_lax_bdf(
        bdf_filename, punch=punch, xref=False,
        is_strict_card_parser=is_strict_card_parser,
        log=log)

    from pyNastran.bdf.mesh_utils.solid_dof import solid_dof
    model, out_nids = solid_dof(model, nid_filename=out_filename, spc_id=spc_id)
    model.log.info('done')
    return model, out_nids
