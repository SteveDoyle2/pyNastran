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
from io import StringIO
from typing import Any
from cpylog import SimpleLogger
import pyNastran
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber, superelement_renumber
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.pierce_shells import pierce_shell_model
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

# testing these imports are up to date
# if something is imported and tested, it should be removed from here
from pyNastran.bdf.mesh_utils.shift import update_nodes
from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric
#from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.dev.create_vectorized_numbered import create_vectorized_numbered
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.free_faces import write_skin_solid_faces
from pyNastran.bdf.mesh_utils.get_oml import get_oml_eids

from .cmd_line.bdf_merge import cmd_line_merge
from .cmd_line.bdf_equivalence import cmd_line_equivalence
from .cmd_line.export_caero_mesh import cmd_line_export_caero_mesh
from .cmd_line.create_flutter import cmd_line_create_flutter
from .cmd_line.utils import filter_no_args


def cmd_line_create_vectorized_numbered(argv=None, quiet=False):  # pragma: no cover
    if argv is None:
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf create_vectorized_numbered IN_BDF_FILENAME [OUT_BDF_FILENAME]\n'
        '  bdf create_vectorized_numbered -h | --help\n'
        '  bdf create_vectorized_numbered -v | --version\n'
        '\n'
        'Positional Arguments:\n'
        '  IN_BDF_FILENAME   the model to convert\n'
        "  OUT_BDF_FILENAME  the converted model name (default=IN_BDF_FILENAME + '_convert.bdf')"
        '\n'
        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(argv) == 1:
        sys.exit(msg)

    from docopt import docopt
    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    bdf_filename_in = data['IN_BDF_FILENAME']
    if data['OUT_BDF_FILENAME']:
        bdf_filename_out = data['OUT_BDF_FILENAME']
    else:
        base, ext = os.path.splitext(bdf_filename_in)
        bdf_filename_out = base + '_convert' + ext
    create_vectorized_numbered(bdf_filename_in, bdf_filename_out)

def cmd_line_delete_bad_shells(argv=None, quiet: bool=False) -> None:
    """command line interface to ``delete_bad_shells``"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf delete_bad_shells IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--skew SKEW] [--max_theta MAX_THETA] [--min_theta MIN_THETA] [--max_ar MAX_AR] [--max_taper MAX_TAPER] [--max_warp MAX_WARP]\n'

        '  bdf delete_bad_shells -h | --help\n'
        '  bdf delete_bad_shells -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME        path to input BDF/DAT/NAS file\n"
        #"  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
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
    try:
        bdf_filename = data['IN_BDF_FILENAME']
    except:
        if not quiet:  # pragma: no cover
            print(data)
        raise

    bdf_out_filename = 'fixed_quality.bdf'

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
    # lengths. For flat elements (such that all of the grid points lie in a plane), this factor is zero.

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
    #from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
    from pyNastran.bdf.bdf import read_bdf

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, validate=True, xref=True, punch=False,
                     encoding=None, log=log, debug=True, mode='msc')
    delete_bad_shells(model,
                      min_theta=min_theta, max_theta=max_theta, max_skew=skew,
                      max_aspect_ratio=max_aspect_ratio, max_taper_ratio=max_taper_ratio,
                      max_warping=max_warping)
    model.write_bdf(bdf_out_filename, size=size,
                    nodes_size=16, elements_size=16, loads_size=8)

def _apply_float_values_to_dict(data: dict[str, Any], defaults: dict[str, float]) -> None:
    for name, default_value in defaults.items():
        if data[name] is None:
            #print(f'applying {name}')
            data[name] = default_value

def cmd_line_bin(argv=None, quiet=False):  # pragma: no cover
    """bins the model into nbins"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        "Usage:\n"
        #"  bdf bin IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n"
        "  bdf bin IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--nbins NBINS]\n"
        '  bdf bin -h | --help\n'
        '  bdf bin -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        "  AXIS1              axis to loop over\n"
        "  AXIS2              axis to bin\n"
        '\n'

        'Options:\n'
        "  --cid CID   the coordinate system to bin (default:0)\n"
        "  --step SIZE   the step size for binning\n\n"
        "  --nbins NBINS the number of bins\n\n"

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n\n"

        'Plot z (2) as a function of y (1) in y-stepsizes of 0.1:\n'
        '  bdf bin fem.bdf 1 2 --cid 0 --step 0.1\n\n'

        'Plot z (2) as a function of y (1) with 50 bins:\n'
        '  bdf bin fem.bdf 1 2 --cid 0 --nbins 50'
    )

    if len(argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    bdf_filename = data['IN_BDF_FILENAME']
    axis1 = int(data['AXIS1'])
    axis2 = int(data['AXIS2'])
    cid = 0
    if data['--cid']:
        cid = int(data['--cid'])

    #stepsize = 0.1
    #if data['--step']:
        #stepsize = float(data['--step'])

    nbins = 10
    if data['--nbins']:
        nbins = int(data['--nbins'])
    assert nbins >= 2, nbins
    if not quiet:  # pragma: no cover
        print(data)

    import numpy as np
    import matplotlib.pyplot as plt
    from pyNastran.bdf.bdf import read_bdf
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)

    model = read_bdf(bdf_filename, log=log)
    xyz_cid = model.get_xyz_in_coord(cid=cid, fdtype='float64')
    y = xyz_cid[:, axis1]
    z = xyz_cid[:, axis2]

    plt.figure(1)
    #n, bins, patches = plt.hist( [x0,x1,x2], 10, weights=[w0, w1, w2], histtype='bar')
    ys = []
    #zs = []
    zs_min = []
    zs_max = []
    y0 = y.min()
    y1 = y.max()
    dy = (y1 - y0) / nbins
    y0i = y0
    y1i = y0 + dy
    for unused_i in range(nbins):
        j = np.where((y0i <= y) & (y <= y1i))[0]
        if not len(j):
            continue
        ys.append(y[j].mean())
        zs_min.append(z[j].min())
        zs_max.append(z[j].max())
        y0i += dy
        y1i += dy
    zs_max = np.array(zs_max)
    zs_min = np.array(zs_min)
    if not quiet:  # pragma: no cover
        print('ys = %s' % ys)
        print('zs_max = %s' % zs_max)
        print('zs_min = %s' % zs_min)
    plt.plot(ys, zs_max, 'r-o', label='max')
    plt.plot(ys, zs_min, 'b-o', label='min')
    plt.plot(ys, zs_max - zs_min, 'g-o', label='delta')
    #plt.xlim([y0, y1])
    plt.xlabel('Axis %s' % axis1)
    plt.ylabel('Axis %s' % axis2)
    plt.grid(True)
    plt.legend()
    plt.show()



def cmd_line_renumber(argv=None, quiet=False):
    """command line interface to bdf_renumber"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        "Usage:\n"
        '  bdf renumber IN_BDF_FILENAME OUT_BDF_FILENAME [--superelement] [--size SIZE]\n'
        '  bdf renumber IN_BDF_FILENAME                  [--superelement] [--size SIZE]\n'
        '  bdf renumber -h | --help\n'
        '  bdf renumber -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '--superelement  calls superelement_renumber\n'
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
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['OUT_BDF_FILENAME']
    if bdf_filename_out is None:
        bdf_filename_out = 'renumber.bdf'

    size = 16
    if data['--size']:
        if 'SIZE' in data:
            size_str = data['SIZE']
        else:
            size_str = data['--size']
        size = int(size_str)

    assert size in [8, 16], f'size={size} args={argv}'
    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    if data['--superelement']:
        superelement_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                              starting_id_dict=None, #round_ids=False,
                              cards_to_skip=cards_to_skip, log=log)
    else:
        bdf_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                     starting_id_dict=None, round_ids=False,
                     cards_to_skip=cards_to_skip, log=log)


def cmd_line_mirror(argv=None, quiet: bool=False):
    """command line interface to write_bdf_symmetric"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        "  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--plane PLANE] [--tol TOL]\n"
        "  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--plane PLANE] [--noeq]\n"
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
        "\n" #  (default=0.000001)

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
    punch = data['--punch']
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'mirrored.bdf'

    #from io import StringIO
    from pyNastran.bdf.bdf import read_bdf, BDF
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    model = read_bdf(bdf_filename, punch=punch, log=log)
    #model.read_bdf(bdf_filename, validate=True, xref=False, punch=punch, read_includes=True, save_file_structure=False, encoding=None)

    #grids = {}
    #for set_id, seti in model.sets.items():
        #for i in seti.ids:
            #if i not in grids:
                ##x = set_id + float(i)
                #y = float(i)
                #grids[i] = f'GRID,{i:d},0,0.,{y},1.'
    #for i, grid in sorted(grids.items()):
        #print(grid)
    #model.cross_reference(xref=True, xref_nodes=True, xref_elements=True, xref_nodes_with_elements=False, xref_properties=True,
                          #xref_masses=True, xref_materials=True, xref_loads=True, xref_constraints=True, xref_aero=True, xref_sets=False, xref_optimization=True, word='')
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
        if eid_offset == 0:
            model.log.info('writing mirrored model %s without equivalencing because there are no elements' % bdf_filename_out)
        else:
            model.log.info('writing mirrored model %s without equivalencing' % bdf_filename_out)
        with open(bdf_filename_out, 'w') as bdf_file:
            bdf_file.write(bdf_filename_stringio.getvalue())


def cmd_line_flip_shell_normals(argv=None, quiet: bool=False):
    """command line interface to flip_shell_normals"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        "  bdf flip_shell_normals IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n"
        "  bdf flip_shell_normals IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n"
        '  bdf flip_shell_normals -h | --help\n'
        '  bdf flip_shell_normals -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        #"  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '  --punch                             flag to identify a *.pch/*.inc file\n'
        "\n" #  (default=0.000001)

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

    if not quiet:  # pragma: no cover
        print(data)
    size = 16
    punch = data['--punch']
    zero_zoffset = data['--zero_zoffset']
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'flipped_shell_normals.bdf'

    #from io import StringIO
    from pyNastran.bdf.bdf import read_bdf, BDF

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    model = read_bdf(bdf_filename, punch=punch, log=log, xref=False)

    skip_elements = {
        # nothing to convert (verified)
        'CCONEAX',
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
        'CVISC', 'CBUSH', 'CBUSH1D', 'CBUSH2D',
        'CROD', 'CTUBE',
        'CBAR', 'CBEAM',
        'CSHEAR', 'CQUADX', 'CTRIAX', 'CTRIAX6',
        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        'CRAC2D', 'CRAC3D',
        'CHBDYG', 'CHBDYE', 'CHBDYP',

        # TODO: NX-verify
        'CTRAX3', 'CTRAX6',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
        'CQUADX4', 'CQUADX8',

        # acoustic
        'CHACAB',

    }

    shells = {
        'CTRIA3', 'CTRIA6', 'CTRIAR',
        'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
    }
    for eid, elem in model.elements.items():
        elem_type = elem.type
        if elem_type in shells:
            elem.flip_normal()
            if zero_zoffset:
                elem.zoffset = 0.
        elif elem_type in skip_elements:
            pass
        else:
            log.warning(f'cannot flip {elem_type}')

    model.write_bdf(bdf_filename_out, encoding=None,
                    size=size, nodes_size=16, elements_size=8, loads_size=8,
                    is_double=False, interspersed=False, enddata=None, write_header=True, close=True)


def cmd_line_convert(argv=None, quiet=False):
    """command line interface to bdf_merge"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        "Usage:\n"
        '  bdf convert IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
        '  bdf convert -h | --help\n'
        '  bdf convert -v | --version\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n'
        '  --in_units  IN_UNITS                length,mass\n'
        '  --out_units  OUT_UNITS              length,mass\n\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n\n"

        'Example:\n'
        '  bdf convert model.bdf --in_units m,kg  --out_units in,lbm\n'
        '  bdf convert model.bdf --in_units m,kg  --out_units in,slinch\n'
        '  bdf convert model.bdf --in_units m,kg  --out_units ft,slug\n'
        '  bdf convert model.bdf --in_units m,kg  --out_units ft,lbm\n'
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
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        #bdf_filename_out = 'merged.bdf'
        bdf_filename_out = bdf_filename + '.convert.bdf'

    in_units = data['IN_UNITS']
    if in_units is None:
        in_units = 'm,kg'

    out_units = data['OUT_UNITS']
    if out_units is None:
        out_units = 'm,kg'

    length_in, mass_in = in_units.split(',')
    length_out, mass_out = out_units.split(',')
    units_to = [length_out, mass_out, 's']
    units = [length_in, mass_in, 's']

    length_to_mass = {
        'in': {'lbm', 'slinch'},
        'ft': {'lbm', 'slug'},
        'm': {'g', 'kg', 'Mg'},
        'cm': {'g', 'kg', 'Mg'},
        'mm': {'g', 'kg', 'Mg'},
    }
    length_allowed = {'in', 'ft', 'm', 'cm', 'mm'}
    assert length_in in length_allowed, f'mass_out={mass_out!r} allowed={length_allowed}'
    assert length_out in {'in', 'ft', 'm', 'cm', 'mm'}, f'mass_out={mass_out!r} allowed={length_allowed}'
    assert mass_in in length_to_mass[length_in], f'mass_out={mass_out!r} allowed={length_to_mass[length_in]}'
    assert mass_out in length_to_mass[length_out], f'mass_out={mass_out!r} allowed={length_to_mass[length_out]}'

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.convert import convert

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, validate=True, xref=True,
                     punch=False, save_file_structure=False,
                     skip_cards=None, read_cards=None,
                     encoding=None, log=log, debug=True, mode='msc')
    convert(model, units_to, units=units)
    for prop in model.properties.values():
        prop.comment = ''
    model.write_bdf(bdf_filename_out)


def cmd_line_scale(argv=None, quiet=False):
    if argv is None:
        argv = sys.argv

    import argparse
    #import textwrap
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # positional arguments
    parent_parser.add_argument('scale', type=str)
    parent_parser.add_argument('INPUT', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)

    #'  --l  LENGTH_SF                    length scale factor\n'
    #'  --m  MASS_SF                      mass scale factor\n'
    #'  --f  FORCE_SF                     force scale factor\n'
    #'  --p  PRESSURE_SF                  pressure scale factor\n'
    #'  --t  TIME_SF                      time scale factor\n'
    #'  --v  VEL_SF                       velocity scale factor\n'

    parent_parser.add_argument('-l', '--length', help='length scale factor')
    parent_parser.add_argument('-m', '--mass', help='mass scale factor')
    parent_parser.add_argument('-f', '--force', help='force scale factor')
    parent_parser.add_argument('-p', '--pressure', help='pressure scale factor')
    parent_parser.add_argument('-t', '--time', help='time scale factor')
    parent_parser.add_argument('-V', '--velocity', help='velocity scale factor')
    #parent_parser.add_argument('--user_geom', type=str, help='log msg')


    #parent_parser.add_argument('-q', '--quiet', help='prints debug messages (default=True)', action='store_true')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    scales = []
    terms = []
    bdf_filename = args.INPUT
    bdf_filename_out = args.OUTPUT
    if bdf_filename_out is None:
        bdf_filename_base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = '%s.scaled%s' % (bdf_filename_base, ext)

    #assert bdf_filename_out is not None
    if args.mass:
        scale = float(args.mass)
        scales.append(scale)
        terms.append('M')
    if args.length:
        scale = float(args.length)
        scales.append(scale)
        terms.append('L')
    if args.time:
        scale = float(args.time)
        scales.append(scale)
        terms.append('T')
    if args.force:
        scale = float(args.force)
        scales.append(scale)
        terms.append('F')
    if args.pressure:
        scale = float(args.pressure)
        scales.append(scale)
        terms.append('P')
    if args.velocity:
        scale = float(args.velocity)
        scales.append(scale)
        terms.append('V')

    from pyNastran.bdf.mesh_utils.convert import scale_by_terms

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    scale_by_terms(bdf_filename, terms, scales, bdf_filename_out=bdf_filename_out, log=log)


def cmd_line_export_mcids(argv=None, quiet=False):
    """command line interface to export_mcids"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf export_mcids IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--iplies PLIES] [--no_x | --no_y]\n'
        '  bdf export_mcids -h | --help\n'
        '  bdf export_mcids -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_CSV_FILENAME  path to output CSV file\n'
        '  --iplies PLIES                      the plies indices to export; comma separated (default=0)\n'
        '\n'

        'Data Suppression:\n'
        "  --no_x,  don't write the x axis\n"
        "  --no_y,  don't write the y axis\n"
        '\n'

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
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    csv_filename_in = data['--output']
    if csv_filename_in is None:
        csv_filename_in = 'mcids.csv'

    export_xaxis = True
    export_yaxis = True
    if data['--no_x']:
        export_xaxis = False
    if data['--no_y']:
        export_yaxis = False
    csv_filename_base = os.path.splitext(csv_filename_in)[0]
    iplies = [0]
    if data['--iplies']:
        iplies = data['--iplies'].split(',')
        iplies = [int(iply) for iply in iplies]
        if not quiet:  # pragma: no cover
            print('iplies = %s' % iplies)

    from pyNastran.bdf.bdf import read_bdf

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log, xref=False)
    model.safe_cross_reference()

    for iply in iplies:
        csv_filename = csv_filename_base + '_ply=%i.csv' % iply
        export_mcids(model, csv_filename,
                     export_xaxis=export_xaxis, export_yaxis=export_yaxis, iply=iply)
        model.log.info('wrote %s' % csv_filename)

def cmd_line_remove_unused(argv=None, quiet=False):
    """command line interface to remove_unused"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf remove_unused IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf remove_unused -h | --help\n'
        '  bdf remove_unused -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF file\n'
        '\n'

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
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    out_bdf_filename = data['--output']
    if out_bdf_filename is None:
        abs_name = os.path.abspath(bdf_filename)
        dirname = os.path.dirname(abs_name)
        basename = os.path.basename(abs_name)
        out_bdf_filename = os.path.join(dirname, 'clean_' + basename)

    from pyNastran.bdf.bdf import read_bdf

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log, xref=False)
    #model.cross_reference()
    remove_unused(model,
                  remove_nids=True, remove_cids=True,
                  remove_pids=True, remove_mids=True,
                  remove_spcs=True, remove_mpcs=True,
                  remove_optimization=True,
                  reset_type_to_id_map=False)
    model.write_bdf(out_bdf_filename,
                    nodes_size=None,
                    is_double=False, interspersed=False)
    #for iply in iplies:
        #csv_filename = csv_filename_base + '_ply=%i.csv' % iply
        #export_mcids(model, csv_filename,
                     #export_xaxis=export_xaxis, export_yaxis=export_yaxis, iply=iply)
        #model.log.info('wrote %s' % csv_filename)

def cmd_line_free_faces(argv=None, quiet=False):
    """command line interface to bdf free_faces"""
    if argv is None:
        argv = sys.argv

    encoding = sys.getdefaultencoding()
    usage = (
        'Usage:\n'
        '  bdf free_faces BDF_FILENAME SKIN_FILENAME [-d] [-l] [-f] [--encoding ENCODE]\n'
        '  bdf free_faces -h | --help\n'
        '  bdf free_faces -v | --version\n'
        '\n'
    )
    arg_msg = (
        "Positional Arguments:\n"
        "  BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        "  SKIN_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        '  -l, --large        writes the BDF in large field, single precision format (default=False)\n'
        '  -d, --double       writes the BDF in large field, double precision format (default=False)\n'
        f'  --encoding ENCODE  the encoding method (default=None -> {encoding!r})\n'
        '\n'
        'Developer:\n'
        '  -f, --profile    Profiles the code (default=False)\n'
        '\n'
        "Info:\n"
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n"
    )
    filter_no_args(arg_msg, argv, quiet=quiet)

    arg_msg += '\n'

    examples = (
        'Examples\n'
        '--------\n'
        '  bdf free_faces solid.bdf skin.bdf\n'
        '  bdf free_faces solid.bdf skin.bdf --large\n'
    )
    import argparse
    parent_parser = argparse.ArgumentParser()
    # positional arguments
    parent_parser.add_argument('BDF_FILENAME', help='path to input BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('SKIN_FILENAME', help='path to output BDF/DAT/NAS file', type=str)

    size_group = parent_parser.add_mutually_exclusive_group()
    size_group.add_argument('-d', '--double', help='writes the BDF in large field, single precision format', action='store_true')
    size_group.add_argument('-l', '--large', help='writes the BDF in large field, double precision format', action='store_true')
    size_group.add_argument('--encoding', help='the encoding method (default=None -> {repr(encoding)})', type=str)
    parent_parser.add_argument('--profile', help='Profiles the code', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version', version=pyNastran.__version__)

    from pyNastran.utils.arg_handling import argparse_to_dict, update_message

    update_message(parent_parser, usage, arg_msg, examples)
    if not quiet:
        print(argv)
    args = parent_parser.parse_args(args=argv[2:])
    data = argparse_to_dict(args)

    if not quiet:  # pragma: no cover
        for key, value in sorted(data.items()):
            print("%-12s = %r" % (key.strip('--'), value))

    import time
    time0 = time.time()

    is_double = False
    if data['double']:
        size = 16
        is_double = True
    elif data['large']:
        size = 16
    else:
        size = 8
    bdf_filename = data['BDF_FILENAME']
    skin_filename = data['SKIN_FILENAME']


    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    tol = 1e-005
    bdf_filename_merged = 'merged.bdf'
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    bdf_equivalence_nodes(bdf_filename, bdf_filename_merged, tol,
                          renumber_nodes=False, neq_max=10, xref=True,
                          node_set=None,
                          size=8, is_double=is_double,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False, log=log, debug=True)
    if not quiet:  # pragma: no cover
        print('done with equivalencing')
    write_skin_solid_faces(
        bdf_filename_merged, skin_filename,
        write_solids=False, write_shells=True,
        size=size, is_double=is_double, encoding=None, log=log,
    )
    if not quiet:  # pragma: no cover
        print('total time:  %.2f sec' % (time.time() - time0))


def cmd_line_split_cbars_by_pin_flag(argv=None, quiet=False):
    """command line interface to split_cbars_by_pin_flag"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
        '  bdf split_cbars_by_pin_flags -h | --help\n'
        '  bdf split_cbars_by_pin_flags -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file\n'
        ' -p PIN, --pin     PIN_FLAGS_CSV_FILENAME  path to pin_flags_csv file\n'
        '\n'

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
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    bdf_filename_in = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'model_new.bdf'

    pin_flags_filename = data['--pin']
    if pin_flags_filename is None:
        pin_flags_filename = 'pin_flags.csv'

    split_cbars_by_pin_flag(bdf_filename_in, pin_flags_filename=pin_flags_filename,
                            bdf_filename_out=bdf_filename_out)

def cmd_line_transform(argv=None, quiet=False):
    """command line interface to export_caero_mesh"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf transform IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--shift XYZ]\n'
        '  bdf transform -h | --help\n'
        '  bdf transform -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file\n'
        '\n'

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
    if not quiet:  # pragma: no cover
        print(data)

    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'transform.bdf'

    dxyz = None
    import numpy as np
    if data['--shift']:
        dxyz = np.array(data['XYZ'].split(','), dtype='float64')
        assert len(dxyz) == 3, dxyz

    from pyNastran.bdf.bdf import read_bdf

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log)

    nid_cp_cd, xyz_cid0, unused_xyz_cp, unused_icd_transform, unused_icp_transform = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int32')

    update_nodes_flag = False
    # we pretend to change the SPOINT location
    if dxyz is not None:
        xyz_cid0 += dxyz
        update_nodes_flag = True

    if update_nodes_flag:
        update_nodes(model, nid_cp_cd, xyz_cid0)
        model.write_bdf(bdf_filename_out)

def cmd_line_filter(argv=None, quiet=False):  # pragma: no cover
    """command line interface to bdf filter"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    msg = (
        'Usage:\n'
        '  bdf filter IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf filter IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--x YSIGN_X] [--y YSIGN_Y] [--z YSIGN_Z]\n'
        '  bdf filter -h | --help\n'
        '  bdf filter -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME    path to output BDF file (default=filter.bdf)\n'
        " --x YSIGN_X                           a string (e.g., '< 0.')\n"
        " --y YSIGN_Y                           a string (e.g., '< 0.')\n"
        " --z YSIGN_Z                           a string (e.g., '< 0.')\n"
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
        '\n'
        'Examples\n'
        '1. remove unused cards:\n'
        '   >>> bdf filter fem.bdf'
        '2. remove GRID points and associated cards with y value < 0:\n'
        "   >>> bdf filter fem.bdf --y '< 0.'"
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'filter.bdf'

    import numpy as np
    func_map = {
        '<' : np.less,
        '>' : np.greater,
        '<=' : np.less_equal,
        '>=' : np.greater_equal,
    }
    xsign = None
    ysign = None
    zsign = None
    if data['--x']:
        xsign, xval = data['--x'].split(' ')
        xval = float(xval)
        assert xsign in ['<', '>', '<=', '>='], xsign
    if data['--y']: # --y < 0
        ysign, yval = data['--y'].split(' ')
        yval = float(yval)
        assert ysign in ['<', '>', '<=', '>='], ysign
    if data['--z']:
        zsign, zval = data['--z'].split(' ')
        zval = float(zval)
        assert zsign in ['<', '>', '<=', '>='], zsign

    from pyNastran.bdf.bdf import read_bdf

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log)

    #nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array(
        #cid=0, fdtype='float64', idtype='int32')

    eids = []
    xyz_cid0 = []
    for eid, elem in sorted(model.elements.items()):
        xyz = elem.Centroid()
        xyz_cid0.append(xyz)
        eids.append(eid)
    xyz_cid0 = np.array(xyz_cid0)
    eids = np.array(eids)

    # we pretend to change the SPOINT location
    update_nodesi = False
    # we pretend to change the SPOINT location
    iunion = None
    if xsign:
        xvals = xyz_cid0[:, 0]
        xfunc = func_map[xsign]
        ix = xfunc(xvals, xval)
        iunion = _union(xval, ix, iunion)
        update_nodesi = True
    if ysign:
        yvals = xyz_cid0[:, 1]
        yfunc = func_map[ysign]
        iy = yfunc(yvals, yval)
        iunion = _union(yval, iy, iunion)
        update_nodesi = True
    if zsign:
        zvals = xyz_cid0[:, 2]
        zfunc = func_map[zsign]
        iz = zfunc(zvals, zval)
        iunion = _union(zval, iz, iunion)
        update_nodesi = True

    if update_nodesi:
        eids_to_remove = eids[iunion]
        for eid in eids_to_remove:
            etype = model.elements[eid].type
            model._type_to_id_map[etype].remove(eid)
            del model.elements[eid]

    #update_nodes(model, nid_cp_cd, xyz_cid0)
    # unxref'd model
    remove_unused(model, remove_nids=True, remove_cids=True,
                  remove_pids=True, remove_mids=True)
    model.write_bdf(bdf_filename_out)


def _union(xval, iunion, ix):
    """helper method for ``filter``"""
    import numpy as np
    if xval:
        if iunion:
            iunion = np.union1d(iunion, ix)
        else:
            pass
    return iunion

def cmd_line(argv=None, quiet: bool=False):
    """command line interface to multiple other command line scripts"""
    if argv is None:
        argv = sys.argv

    dev = True
    msg = (
        'Usage:\n'
        '  bdf merge                       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
        '  bdf equivalence                 IN_BDF_FILENAME EQ_TOL\n'
        '  bdf renumber                    IN_BDF_FILENAME [OUT_BDF_FILENAME] [--superelement] [--size SIZE]\n'
        '  bdf remove_unused               IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf filter                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--x YSIGN X] [--y YSIGN Y] [--z YSIGN Z]\n'
        '  bdf delete_bad_shells           IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--skew SKEW] [--max_theta MAX_THETA] [--min_theta MIN_THETA] [--max_ar MAX_AR] [--max_taper MAX_TAPER] [--max_warp MAX_WARP]\n'
        '  bdf mirror                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
        '  bdf convert                     IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
        '  bdf scale                       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--lsf LENGTH_SF] [--msf MASS_SF] [--fsf FORCE_SF] [--psf PRESSURE_SF] [--tsf TIME_SF] [--vsf VEL_SF]\n'
        '  bdf export_mcids                IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x | --no_y]\n'
        '  bdf free_faces                  BDF_FILENAME SKIN_FILENAME [-d | -l] [-f] [--encoding ENCODE]\n'
        '  bdf flutter                     UNITS eas EAS1 EAS2 SWEEP_UNIT N CONST_TYPE CONST_VAL [-o OUT_BDF_FILENAME] [--size SIZE | --clean]'
        '  bdf flip_shell_normals          IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n'
        '  bdf transform                   IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--shift XYZ]\n'
        '  bdf export_caero_mesh           IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--subpanels] [--pid PID]\n'
        '  bdf split_cbars_by_pin_flags    IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
    )

    if dev:
        msg += '  bdf create_vectorized_numbered  IN_BDF_FILENAME [OUT_BDF_FILENAME]\n'
        msg += '  bdf bin                         IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n'

    msg += (
        #'\n'
        '  bdf merge              -h | --help\n'
        '  bdf equivalence        -h | --help\n'
        '  bdf renumber           -h | --help\n'
        '  bdf remove_unused      -h | --help\n'
        '  bdf delete_bad_shells  -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf mirror             -h | --help\n'
        '  bdf convert            -h | --help\n'
        '  bdf scale              -h | --help\n'
        '  bdf export_mcids       -h | --help\n'
        '  bdf free_faces         -h | --help\n'
        '  bdf flip_shell_normals -h | --help\n'
        '  bdf transform          -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf flutter            -h | --help\n'
        '  bdf export_caero_mesh  -h | --help\n'
        '  bdf split_cbars_by_pin_flags  -h | --help\n'
    )
    if dev:
        msg += (
            '  bdf create_vectorized_numbered  -h | --help\n'
            '  bdf bin                         -h | --help\n'
        )
    msg += '  bdf -v | --version\n'
    msg += '\n'

    filter_no_args(msg + 'Not enough arguments.\n', argv, quiet=quiet)

    #assert sys.argv[0] != 'bdf', msg

    method = argv[1]

    maps = {
        'merge': cmd_line_merge,
        'equivalence': cmd_line_equivalence,
        'renumber': cmd_line_renumber,
        'mirror': cmd_line_mirror,
        'convert': cmd_line_convert,
        'delete_bad_shells': cmd_line_delete_bad_shells,
        'scale': cmd_line_scale,
        'export_mcids': cmd_line_export_mcids,
        'remove_unused': cmd_line_remove_unused,
        'split_cbars_by_pin_flags': cmd_line_split_cbars_by_pin_flag,

        'export_caero_mesh': cmd_line_export_caero_mesh,
        'transform': cmd_line_transform,
        'filter': cmd_line_filter,
        'free_faces': cmd_line_free_faces,
        'flip_shell_normals': cmd_line_flip_shell_normals,
    }
    if dev:
        maps.update({
            'bin': cmd_line_bin,
            'create_vectorized_numbered': cmd_line_create_vectorized_numbered,
            'flutter': cmd_line_create_flutter,
        })

    if method in ['-v', '--version']:
        print(pyNastran.__version__)
    else:
        try:
            func = maps[method]
        except KeyError:
            print(argv)
            sys.exit(msg)
        func(argv, quiet=quiet)


if __name__ == '__main__':  # pragma: no cover
    # for the exe, we pass all the args, but we hack them to have the bdf prefix
    from copy import deepcopy
    argv = deepcopy(sys.argv)
    argv[0] = 'bdf'
    cmd_line(argv=argv)
