"""
defines:
    bdf merge        (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    bdf equivalence  IN_BDF_FILENAME EQ_TOL\n'
    bdf renumber     IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf mirror       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
    bdf export_mcids IN_BDF_FILENAME [-o OUT_GEOM_FILENAME]\n'
    bdf split_cbars_by_pin_flags IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf flutter UNITS [-o OUT_BDF_FILENAME]

"""
import os
import sys
from cpylog import SimpleLogger

# testing these imports are up to date
# if something is imported and tested, it should be removed from here
import pyNastran
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris

def cmd_line_equivalence(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_equivalence_nodes"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf equivalence IN_BDF_FILENAME EQ_TOL [-o OUT_BDF_FILENAME]\n'
        '  bdf equivalence -h | --help\n'
        '  bdf equivalence -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME   path to input BDF/DAT/NAS file\n"
        "  EQ_TOL            the spherical equivalence tolerance\n"
        #"  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

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
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        dirname = os.path.dirname(bdf_filename)
        bdf_filename_out = os.path.join(dirname, 'merged.bdf')
    else:
        dirname = os.path.dirname(bdf_filename_out)

    tol = float(data['EQ_TOL'])
    size = 16
    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
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
