"""
defines:
    bdf merge        (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    bdf equivalence  IN_BDF_FILENAME EQ_TOL\n'
    bdf renumber     IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf mirror       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
    bdf export_mcids IN_BDF_FILENAME [-o OUT_GEOM_FILENAME]\n'
    bdf split_cbars_by_pin_flags IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'

"""
from __future__ import print_function
import os
import sys
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.pierce_shells import pierce_shell_model

# testing these imports are up to date
from pyNastran.bdf.mesh_utils.shift import update_nodes
from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.dev.create_vectorized_numbered import create_vectorized_numbered
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

def cmd_line_create_vectorized_numbered():  # pragma: no cover
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
    if len(sys.argv) == 1:
        sys.exit(msg)

    from docopt import docopt
    import pyNastran
    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    print(data)
    bdf_filename_in = data['IN_BDF_FILENAME']
    if data['OUT_BDF_FILENAME']:
        bdf_filename_out = data['OUT_BDF_FILENAME']
    else:
        base, ext = os.path.splitext(bdf_filename_in)
        bdf_filename_out = base + '_convert' + ext
    create_vectorized_numbered(bdf_filename_in, bdf_filename_out)


def cmd_line_equivalence():  # pragma: no cover
    """command line interface to bdf_equivalence_nodes"""
    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        "  bdf equivalence IN_BDF_FILENAME EQ_TOL  [-o OUT_BDF_FILENAME]\n"

        '  bdf equivalence -h | --help\n'
        '  bdf equivalence -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME   path to input BDF/DAT/NAS file\n"
        "  EQ_TOL            the spherical equivalence tolerance\n\n"
        #"  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'merged.bdf'
    tol = data['EQ_TOL']
    size = 16
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol,
                          renumber_nodes=False,
                          neq_max=10, xref=True,
                          node_set=None, size=size,
                          is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False,
                          debug=True)

def cmd_line_bin():  # pragma: no cover
    """bins the model into nbins"""
    from docopt import docopt
    import pyNastran
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

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
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
    print(data)
    #asdf

    import numpy as np
    import matplotlib.pyplot as plt
    from pyNastran.bdf.bdf import read_bdf

    model = read_bdf(bdf_filename)
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
    for i in range(nbins):
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



def cmd_line_renumber():  # pragma: no cover
    """command line interface to bdf_renumber"""
    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        "  bdf renumber IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n"
        '  bdf renumber -h | --help\n'
        '  bdf renumber -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
       #"  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'renumber.bdf'

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []
    bdf_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                 starting_id_dict=None, round_ids=False,
                 cards_to_skip=cards_to_skip)

def cmd_line_mirror():  # pragma: no cover
    """command line interface to write_bdf_symmetric"""
    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        "  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n"
        "  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--noeq]\n"
        '  bdf mirror -h | --help\n'
        '  bdf mirror -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
       #"  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"
        "  --plane PLANE                       the symmetry plane (xz, yz, xy)\n\n"
        '  --tol   TOL                         the spherical equivalence tolerance\n'
        '  --noeq                              disable equivalencing\n'
        "\n" #  (default=0.000001)

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    if data['--tol'] is None:
        data['TOL'] = 0.000001

    tol = data['TOL']
    if data['--noeq'] is not None:
        tol = -1.

    plane = data['--plane']

    print(data)
    size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'mirrored.bdf'

    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
    model = read_bdf(bdf_filename)
    size = 16
    bdf_filename_temp = '__temp.bdf__'
    write_bdf_symmetric(model, bdf_filename_temp, encoding=None, size=size,
                        is_double=False,
                        enddata=None, close=True,
                        plane=plane)

    bdf_equivalence_nodes(bdf_filename_temp, bdf_filename_out, tol,
                          renumber_nodes=False,
                          neq_max=10, xref=True,
                          node_set=None, size=size,
                          is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False,
                          debug=True)
    os.remove(bdf_filename_temp)

def cmd_line_merge():  # pragma: no cover
    """command line interface to bdf_merge"""
    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        '  bdf merge (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
        '  bdf merge -h | --help\n'
        '  bdf merge -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAMES   path to input BDF/DAT/NAS files\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
    bdf_filenames = data['IN_BDF_FILENAMES']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'merged.bdf'

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []
    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip)


def cmd_line_export_mcids():  # pragma: no cover
    """command line interface to export_mcids"""
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf export_mcids IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--iplies PLIES] [--no_x] [--no_y]\n'
        '  bdf export_mcids -h | --help\n'
        '  bdf export_mcids -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_CSV_FILENAME  path to output CSV file\n'
        '  --iplies PLIES                      the plies to export; comma separated (default=0)\n'
        '\n'

        'Data Suppression:\n'
        "  --no_x,  don't write the x axis\n"
        "  --no_y,  don't write the y axis\n"
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
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
        print('iplies = %s' % iplies)

    from pyNastran.bdf.bdf import read_bdf
    model = read_bdf(bdf_filename, xref=False) #, log=log, debug=debug)
    model.safe_cross_reference()

    for iply in iplies:
        csv_filename = csv_filename_base + '_ply=%i.csv' % iply
        export_mcids(model, csv_filename,
                     export_xaxis=export_xaxis, export_yaxis=export_yaxis, iply=iply)
        model.log.info('wrote %s' % csv_filename)


def cmd_line_split_cbars_by_pin_flag():  # pragma: no cover
    """command line interface to split_cbars_by_pin_flag"""
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
        '  bdf split_cbars_by_pin_flags -h | --help\n'
        '  bdf split_cbars_by_pin_flags -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        " -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file\n"
        " -p PIN, --pin     PIN_FLAGS_CSV_FILENAME  path to pin_flags_csv file\n\n"
        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
    bdf_filename_in = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'model_new.bdf'

    pin_flags_filename = data['--pin']
    if pin_flags_filename is None:
        pin_flags_filename = 'pin_flags.csv'

    split_cbars_by_pin_flag(bdf_filename_in, pin_flags_filename=pin_flags_filename,
                            bdf_filename_out=bdf_filename_out)

def cmd_line_transform():  # pragma: no cover
    """command line interface to export_caero_mesh"""
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf transform IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--shift XYZ]\n'
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
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
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
    model = read_bdf(bdf_filename)

    nid_cp_cd, xyz_cid0, xyz_cp, icd_transform, icp_transform = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int32')

    update_nodes = False
    # we pretend to change the SPOINT location
    if dxyz is not None:
        xyz_cid0 += dxyz
        update_nodes = True

    if update_nodes:
        update_nodes(model, nid_cp_cd, xyz_cid0)
        model.write_bdf(bdf_filename_out)

def cmd_line_filter():  # pragma: no cover
    """command line interface to export_caero_mesh"""
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf filter IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--x YSIGN_X] [--y YSIGN_Y] [--z YSIGN_Z]\n'
        '  bdf filter -h | --help\n'
        '  bdf filter -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file\n'
        " --x YSIGN_X                                a string (e.g., '< 0.')\n"
        " --y YSIGN_Y                                a string (e.g., '< 0.')\n"
        " --z YSIGN_Z                                a string (e.g., '< 0.')\n"
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
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
    model = read_bdf(bdf_filename)

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
    update_nodes = False
    # we pretend to change the SPOINT location
    iunion = None
    if xsign:
        xvals = xyz_cid0[:, 0]
        xfunc = func_map[xsign]
        ix = xfunc(xvals, xval)
        iunion = _union(xval, ix, iunion)
        update_nodes = True
    if ysign:
        yvals = xyz_cid0[:, 1]
        yfunc = func_map[ysign]
        iy = yfunc(yvals, yval)
        iunion = _union(yval, iy, iunion)
        update_nodes = True
    if zsign:
        zvals = xyz_cid0[:, 2]
        zfunc = func_map[zsign]
        iz = xfunc(zvals, zval)
        iunion = _union(zval, iz, iunion)
        update_nodes = True

    if update_nodes:
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

def cmd_line_export_caero_mesh():  # pragma: no cover
    """command line interface to export_caero_mesh"""
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf export_caero_mesh IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf export_caero_mesh -h | --help\n'
        '  bdf export_caero_mesh -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_CAERO_BDF_FILENAME  path to output BDF file\n'
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    caero_bdf_filename = data['--output']
    if caero_bdf_filename is None:
        caero_bdf_filename = 'caero.bdf'

    from pyNastran.bdf.bdf import read_bdf
    model = read_bdf(bdf_filename)
    model.write_caero_model(caero_bdf_filename)

def cmd_line():  # pragma: no cover
    """command line interface to multiple other command line scripts"""
    dev = True
    msg = (
        'Usage:\n'
        '  bdf merge                       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
        '  bdf equivalence                 IN_BDF_FILENAME EQ_TOL\n'
        '  bdf renumber                    IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf mirror                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
        '  bdf export_mcids                IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x] [--no_y]\n'
        '  bdf transform                   IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--shift XYZ]\n'
        '  bdf export_caero_mesh           IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
        '  bdf split_cbars_by_pin_flags    IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
    )

    if dev:
        msg += '  bdf create_vectorized_numbered  IN_BDF_FILENAME [OUT_BDF_FILENAME]\n'
        msg += '  bdf filter                      IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--x YSIGN X] [--y YSIGN Y] [--z YSIGN Z]\n'
        msg += '  bdf bin                         IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n'

    msg += (
        '\n'
        '  bdf merge              -h | --help\n'
        '  bdf equivalence        -h | --help\n'
        '  bdf renumber           -h | --help\n'
        '  bdf mirror             -h | --help\n'
        '  bdf export_mcids       -h | --help\n'
        '  bdf transform          -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf export_caero_mesh  -h | --help\n'
        '  bdf split_cbars_by_pin_flags  -h | --help\n'
    )
    #bdf create_vectorized_numbered -h | --help
    #bdf create_vectorized_numbered -v | --version


    if dev:
        msg += '  bdf create_vectorized_numbered  -h | --help\n'
        msg += '  bdf filter                      -h | --help\n'
        msg += '  bdf bin                         -h | --help\n'
    msg += '  bdf -v | --version\n'
    msg += '\n'

    if len(sys.argv) == 1:
        sys.exit(msg)

    #assert sys.argv[0] != 'bdf', msg

    if sys.argv[1] == 'merge':
        cmd_line_merge()
    elif sys.argv[1] == 'equivalence':
        cmd_line_equivalence()
    elif sys.argv[1] == 'renumber':
        cmd_line_renumber()
    elif sys.argv[1] == 'mirror':
        cmd_line_mirror()
    elif sys.argv[1] == 'export_mcids':
        cmd_line_export_mcids()
    elif sys.argv[1] == 'split_cbars_by_pin_flags':
        cmd_line_split_cbars_by_pin_flag()
    elif sys.argv[1] == 'export_caero_mesh':
        cmd_line_export_caero_mesh()
    elif sys.argv[1] == 'transform':
        cmd_line_transform()
    elif sys.argv[1] == 'filter' and dev:  # TODO: make better name
        cmd_line_filter()
    elif sys.argv[1] == 'bin' and dev:
        cmd_line_bin()
    elif sys.argv[1] == 'create_vectorized_numbered' and dev:
        cmd_line_create_vectorized_numbered()
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':  # pragma: no cover
    sys.argv = sys.argv[1:]
    cmd_line()
