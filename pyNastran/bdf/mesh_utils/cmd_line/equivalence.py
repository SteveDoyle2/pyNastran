import os
import sys
import argparse
from cpylog import SimpleLogger

# testing these imports are up to date
# if something is imported and tested, it should be removed from here
import pyNastran
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris


def cmd_line_equivalence(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_equivalence_nodes"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    if len(argv) == 1:
        sys.exit("bdf equivalence: use 'bdf equivalence -h' for help")

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf equivalence',
        description='Equivalence nodes in a BDF model',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('EQ_TOL', type=float,
                        help='the spherical equivalence tolerance')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    bdf_filename_out = args.output
    if bdf_filename_out is None:
        dirname = os.path.dirname(bdf_filename)
        bdf_filename_out = os.path.join(dirname, 'merged.bdf')
    else:
        dirname = os.path.dirname(bdf_filename_out)

    tol = args.EQ_TOL
    punch = args.punch
    size = 16
    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    model = read_bdf(bdf_filename, xref=True, punch=punch, log=log, debug=True)
    bdf_equivalence_nodes(model, bdf_filename_out, tol,
                          renumber_nodes=False,
                          neq_max=10, xref=True,
                          node_set=None, size=size,
                          is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False,
                          log=log, debug=True)

    bdf_filename_out2 = os.path.join(dirname, 'merged_collapsed.bdf')
    model = read_bdf(bdf_filename_out, xref=False, validate=False, log=log)
    convert_bad_quads_to_tris(model, eids_to_check=None, xyz_cid0=None, min_edge_length=0.0)
    model.write_bdf(bdf_filename_out2, size=size)
