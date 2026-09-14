import sys
import argparse
from io import StringIO
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args


def cmd_line_mirror(argv=None, quiet: bool=False) -> None:
    """command line interface to write_bdf_symmetric"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf mirror: use 'bdf mirror -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf mirror',
        description='Mirror a BDF model across a symmetry plane',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='mirrored.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file (default=mirrored.bdf)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--plane', default='xz',
                        help='the symmetry plane: xz, yz, xy (default=xz)')
    parser.add_argument('--tol', type=float, default=1e-6,
                        help='the spherical equivalence tolerance (default=1e-6)')
    parser.add_argument('--noeq', action='store_true',
                        help='disable equivalencing')

    args = parser.parse_args(argv[2:])

    tol = args.tol
    if args.noeq:
        tol = -1.

    plane = args.plane

    if not quiet:  # pragma: no cover
        print(vars(args))

    size = 16
    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    log.debug(f'plane = {plane!r}')
    bdf_filename_out = args.output

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
