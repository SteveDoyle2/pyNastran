"""interface to ``bdf ...``"""
from __future__ import annotations
import sys
from typing import TYPE_CHECKING

import pyNastran

from pyNastran.bdf.mesh_utils.run_jobs import cmd_line_run_jobs
from pyNastran.bdf.mesh_utils.host_jobs import cmd_line_host_jobs

from .inclzip import cmd_line_inclzip
from .convert import cmd_line_convert
from .diff import cmd_line_diff
from .delete import cmd_line_delete
from .export_mcids import cmd_line_export_mcids
from .equivalence import cmd_line_equivalence
from .merge import cmd_line_merge
from .renumber import cmd_line_renumber
from .scale import cmd_line_scale
from .remove_unused import cmd_line_remove_unused
from .free_faces import cmd_line_free_faces
from .split_cbar_by_pin_flag import cmd_line_split_cbars_by_pin_flag
from .export_caero_mesh import cmd_line_export_caero_mesh
from .create_flutter import cmd_line_create_flutter
from .list_conm2 import cmd_line_list_conm2
from .mirror import cmd_line_mirror
from .solid_dof import cmd_line_solid_dof
from .delete_bad_shells import cmd_line_delete_bad_shells, SHELL_QUALITY
from .remove_comments import cmd_line_remove_comments
from .stats import cmd_line_stats
from .rbe3_to_rbe2 import cmd_line_rbe3_to_rbe2
from .rbe2_to_rbe3 import cmd_line_rbe2_to_rbe3
from .merge_rbe2 import cmd_line_merge_rbe2
from .collapse_quads import cmd_line_collapse_quads
from .super import cmd_line_super
from .flip_shell_normals import cmd_line_flip_shell_normals
from .nsm_split import cmd_line_nsm_split
from .dev import cmd_line_bin, cmd_line_filter, cmd_line_transform
from .utils import filter_no_args


CMD_MAPS = {
    'inclzip': cmd_line_inclzip,
    'diff': cmd_line_diff,
    'merge': cmd_line_merge,
    'equivalence': cmd_line_equivalence,
    'renumber': cmd_line_renumber,
    'remove_comments': cmd_line_remove_comments,
    'mirror': cmd_line_mirror,
    'convert': cmd_line_convert,
    'delete_bad_shells': cmd_line_delete_bad_shells,
    'collapse_quads': cmd_line_collapse_quads,
    'scale': cmd_line_scale,
    'list_conm2': cmd_line_list_conm2,
    'export_mcids': cmd_line_export_mcids,
    'solid_dof': cmd_line_solid_dof,
    'remove_unused': cmd_line_remove_unused,
    'split_cbars_by_pin_flags': cmd_line_split_cbars_by_pin_flag,
    'run_jobs': cmd_line_run_jobs,
    'host_jobs': cmd_line_host_jobs,
    'super': cmd_line_super,
    'rbe3_to_rbe2': cmd_line_rbe3_to_rbe2,
    'rbe2_to_rbe3': cmd_line_rbe2_to_rbe3,
    'merge_rbe2': cmd_line_merge_rbe2,
    'nsm_split': cmd_line_nsm_split,

    'delete': cmd_line_delete,
    #'summary': cmd_line_summary,

    'export_caero_mesh': cmd_line_export_caero_mesh,
    'transform': cmd_line_transform,
    'filter': cmd_line_filter,
    'free_faces': cmd_line_free_faces,
    'flip_shell_normals': cmd_line_flip_shell_normals,
    'flutter': cmd_line_create_flutter,
    'stats': cmd_line_stats,
}

dev = True
if dev:
    CMD_MAPS.update({
        'bin': cmd_line_bin,
    })

SCALES = (
    '[--length LENGTH_SF] [--mass MASS_SF] [--force FORCE_SF] '
    '[--pressure PRESSURE_SF] [--time TIME_SF] [--velocity VEL_SF]')


def cmd_line(argv=None, log=None, quiet: bool=False) -> None:
    """command line interface to multiple other command line scripts"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    usage_bin_msg = '  bdf bin                         IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n' if dev else ''
    help_bin_msg = '  bdf bin                         -h | --help\n' if dev else ''
    msg = (
        'Usage:\n'
        '  bdf diff                        IN_BDF_FILENAME1 IN_BDF_FILENAME2 [--punch]\n'
        '  bdf merge                       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
        '  bdf equivalence                 IN_BDF_FILENAME EQ_TOL [--punch]\n'
        '  bdf inclzip                     IN_BDF_FILENAME EQ_TOL [-o OUT_BDF_FILENAME] [--punch]\n'
        '  bdf renumber                    IN_BDF_FILENAME [OUT_BDF_FILENAME] [--superelement] [--size SIZE]\n'
        '  bdf remove_unused               IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch]\n'
        '  bdf filter                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--x YSIGN X] [--y YSIGN Y] [--z YSIGN Z]\n'
       f'  bdf delete_bad_shells           IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] {SHELL_QUALITY}\n'
        '  bdf collapse_quads              IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--size SIZE]\n'
        '  bdf mirror                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--plane PLANE] [--tol TOL]\n'
        '  bdf convert                     IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
       f'  bdf scale                       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] {SCALES}\n'
        '  bdf export_mcids                IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x | --no_y]\n'
        '  bdf free_faces                  BDF_FILENAME SKIN_FILENAME [-d | -l] [-f] [--encoding ENCODE]\n'
        '  bdf nsm_split                   BDF_FILENAME NSM_ID NSM_VALUE\n'
        '  bdf flutter                     UNITS eas EAS1 EAS2 SWEEP_UNIT N CONST_TYPE CONST_VAL [-o OUT_BDF_FILENAME] [--size SIZE | --clean]'
        '  bdf flip_shell_normals          IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n'
        '  bdf transform                   IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--shift XYZ]\n'
        '  bdf export_caero_mesh           IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--aerobox] [--pid PID]\n'
        '  bdf split_cbars_by_pin_flags    IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [-p PIN_FLAGS_CSV_FILENAME]\n'
        '  bdf solid_dof                   IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--lax] [--allow_dup]\n'
        '  bdf stats                       IN_BDF_FILENAME [--punch]\n'
        '  bdf rbe3_to_rbe2                IN_BDF_FILENAME [--infile INFILE] [--punch] [--lax]\n'
        '  bdf merge_rbe2                  IN_BDF_FILENAME [--infile INFILE] [--punch] [--lax]\n'
        '  bdf run_jobs                    BDF_FILENAME_DIRNAME [FILE...] [--exe NASTRAN_PATH] [--infile INFILE] [--outfile OUTFILE] [--skip SKIP] [--all] [--test] [--debug] [--cleanup]\n'
       f'{usage_bin_msg}'
    )

    msg += (
        # '\n'
        'Mesh Tools:\n'
        '  bdf merge              -h | --help\n'
        '  bdf equivalence        -h | --help\n'
        '  bdf renumber           -h | --help\n'
        '  bdf remove_unused      -h | --help\n'
        '  bdf delete_bad_shells  -h | --help\n'
        '  bdf collapse_quads     -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf mirror             -h | --help\n'
        '  bdf convert            -h | --help\n'
        '  bdf scale              -h | --help\n'
        '  bdf free_faces         -h | --help\n'
        '  bdf flip_shell_normals -h | --help\n'
        '  bdf transform          -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf rbe3_to_rbe2       -h | --help\n'
        '  bdf merge_rbe2         -h | --help\n'
        'Export Tools:\n'
        '  bdf inclzip            -h | --help\n'
        '  bdf export_mcids       -h | --help\n'
        '  bdf export_caero_mesh  -h | --help\n'
        '  bdf split_cbars_by_pin_flags    -h | --help\n'
        '  bdf solid_dof                   -h | --help\n'
        'Comparision Tools:\n'
        '  bdf diff               -h | --help\n'
        '  bdf stats              -h | --help\n'
        f'{usage_bin_msg}'
        'Analysis Tools:\n'
        '  bdf flutter            -h | --help\n'
        '  bdf run_jobs           -h | --help\n'
    )
    msg += '  bdf -v | --version\n'
    msg += '\n'

    filter_no_args(msg + 'Not enough arguments.\n', argv, quiet=quiet)

    # assert sys.argv[0] != 'bdf', msg
    method = argv[1]
    if method in ['-v', '--version']:
        print(pyNastran.__version__)
    else:
        try:
            func = CMD_MAPS[method]
        except KeyError:
            print(argv)
            print(f'method={method!r} not found')
            sys.exit(msg)
        #print('end of cmd_line')
        return func(argv, quiet=quiet)
