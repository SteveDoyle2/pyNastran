import sys
from io import StringIO
from cpylog import SimpleLogger
from docopt import docopt

import pyNastran
from .utils import (
    get_bdf_filename_punch_log, filter_no_args)


def cmd_line_mirror(argv=None, quiet: bool=False) -> None:
    """command line interface to write_bdf_symmetric"""
    if argv is None:  # pragma: no cover
        argv = sys.argv
    msg = (
        'Usage:\n'
        '  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--plane PLANE] [--tol TOL]\n'
        '  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--plane PLANE] [--noeq]\n'
        '  bdf mirror -h | --help\n'
        '  bdf mirror -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        #"  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '  --punch                             flag to identify a *.pch/*.inc file\n'
        "  --plane PLANE                       the symmetry plane (xz, yz, xy); default=xz\n"
        '  --tol   TOL                         the spherical equivalence tolerance; default=1e-6\n'
        '  --noeq                              disable equivalencing\n'
        "\n"  # (default=0.000001)

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
    if data['--tol'] is False:
        data['TOL'] = 0.000001

    if isinstance(data['TOL'], str):
        data['TOL'] = float(data['TOL'])
    tol = data['TOL']

    assert data['--noeq'] in [True, False]
    if data['--noeq']:
        tol = -1.

    plane = 'xz'
    if data['--plane'] is not None:  # None or str
        plane = data['--plane']

    if not quiet:  # pragma: no cover
        print(data)

    size = 16
    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)
    log.debug(f'plane = {plane!r}')
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'mirrored.bdf'

    #from io import StringIO
    from pyNastran.bdf.bdf import read_bdf, BDF
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
    from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric

    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    model = read_bdf(bdf_filename, punch=punch, log=log)
    # model.read_bdf(bdf_filename, validate=True, xref=False, punch=punch,
    #                read_includes=True, save_file_structure=False, encoding=None)

    # grids = {}
    # for set_id, seti in model.sets.items():
    #     for i in seti.ids:
    #         if i not in grids:
    #             #x = set_id + float(i)
    #             y = float(i)
    #             grids[i] = f'GRID,{i:d},0,0.,{y},1.'
    # for i, grid in sorted(grids.items()):
    #     print(grid)
    # model.cross_reference(
    #     xref=True, xref_nodes=True, xref_elements=True,
    #     xref_nodes_with_elements=False, xref_properties=True,
    #     xref_masses=True, xref_materials=True, xref_loads=True,
    #     xref_constraints=True, xref_aero=True, xref_sets=False,
    #     xref_optimization=True, word='')
    bdf_filename_stringio = StringIO()
    unused_model, unused_nid_offset, eid_offset = write_bdf_symmetric(
        model, bdf_filename_stringio, encoding=None, size=size,
        is_double=False,
        enddata=None, close=False,
        plane=plane, log=log)
    bdf_filename_stringio.seek(0)

    if eid_offset > 0 and tol >= 0.0:
        bdf_equivalence_nodes(bdf_filename_stringio, bdf_filename_out, tol,
                              renumber_nodes=False,
                              neq_max=10, xref=True,
                              node_set=None, size=size,
                              is_double=False,
                              remove_collapsed_elements=False,
                              avoid_collapsed_elements=False,
                              crash_on_collapse=False,
                              debug=True, log=log)
    else:
        log: SimpleLogger = model.log
        if eid_offset == 0:
            log.info(f'writing mirrored model {bdf_filename_out} without equivalencing '
                     'because there are no elements')
        else:
            log.info(f'writing mirrored model {bdf_filename_out} without equivalencing')
        with open(bdf_filename_out, 'w') as bdf_file:
            bdf_file.write(bdf_filename_stringio.getvalue())
