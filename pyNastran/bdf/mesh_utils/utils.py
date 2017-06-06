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

# testing these imports are up to date
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag

def cmd_line_equivalence():  # pragma: no cover
    """command line interface to bdf_equivalence_nodes"""
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf equivalence IN_BDF_FILENAME EQ_TOL  [-o OUT_BDF_FILENAME]\n"

    msg += '  bdf equivalence -h | --help\n'
    msg += '  bdf equivalence -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME   path to input BDF/DAT/NAS file\n"
    msg += "  EQ_TOL            the spherical equivalence tolerance\n\n"
    #msg += "  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  -o OUT, --output OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

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
    msg = "Usage:\n"
    #msg += "  bdf bin IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n"
    msg += "  bdf bin IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--nbins NBINS]\n"
    msg += '  bdf bin -h | --help\n'
    msg += '  bdf bin -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    msg += "  AXIS1              axis to loop over\n"
    msg += "  AXIS2              axis to bin\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  --cid CID   the coordinate system to bin (default:0)\n"
    #msg += "  --step SIZE   the step size for binning\n\n"
    msg += "  --nbins NBINS the number of bins\n\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n\n"

    msg += 'Plot z (2) as a function of y (1) in y-stepsizes of 0.1:\n'
    msg += '  bdf bin fem.bdf 1 2 --cid 0 --step 0.1\n\n'

    msg += 'Plot z (2) as a function of y (1) with 50 bins:\n'
    msg += '  bdf bin fem.bdf 1 2 --cid 0 --nbins 50'

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
    msg = "Usage:\n"
    msg += "  bdf renumber IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n"
    msg += '  bdf renumber -h | --help\n'
    msg += '  bdf renumber -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    #msg += "  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

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

    cards_to_skip = [
        'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    bdf_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                 starting_id_dict=None, round_ids=False,
                 cards_to_skip=cards_to_skip)

def cmd_line_mirror():  # pragma: no cover
    """command line interface to write_bdf_symmetric"""
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf mirror IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n"
    msg += '  bdf mirror -h | --help\n'
    msg += '  bdf mirror -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    #msg += "  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"
    msg += "  --plane PLANE                       the symmetry plane (xz, ???)\n\n"
    msg += "  --tol   TOL                         the spherical equivalence tolerance (default=0.000001)\n\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    if data['--tol'] is None:
        data['--tol'] = 0.000001
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
    model.write_bdf_symmetric(bdf_filename_temp, encoding=None, size=size,
                              is_double=False,
                              enddata=None,
                              close=True,
                              plane='xz')
    tol = 0.000001
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
    msg = "Usage:\n"
    msg += "  bdf merge (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n"
    msg += '  bdf merge -h | --help\n'
    msg += '  bdf merge -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAMES   path to input BDF/DAT/NAS files\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

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

    cards_to_skip = [
        'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip)


def cmd_line_export_mcid():  # pragma: no cover
    """command line interface to export_mcids"""
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf export_mcids IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x] [--no_y]\n"
    msg += '  bdf export_mcids -h | --help\n'
    msg += '  bdf export_mcids -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += "  -o OUT, --output  OUT_CSV_FILENAME  path to output CSV file\n\n"

    msg += 'Data Suppression:\n'
    msg += "  --no_x,  don't write the x axis\n"
    msg += "  --no_y,  don't write the y axis\n"

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

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
    csv_filename = data['--output']
    if csv_filename is None:
        csv_filename = 'mcids.csv'

    export_xaxis = True
    export_yaxis = True
    if data['--no_x']:
        export_xaxis = False
    if data['--no_y']:
        export_yaxis = False
    export_mcids(bdf_filename, csv_filename,
                 export_xaxis=export_xaxis, export_yaxis=export_yaxis)

def cmd_line_split_cbars_by_pin_flag():  # pragma: no cover
    """command line interface to split_cbars_by_pin_flag"""
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += '  bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
    msg += '  bdf split_cbars_by_pin_flags -h | --help\n'
    msg += '  bdf split_cbars_by_pin_flags -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    msg += " -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file\n"
    msg += " -p PIN, --pin     PIN_FLAGS_CSV_FILENAME  path to pin_flags_csv file\n\n"
    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

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

def cmd_line():  # pragma: no cover
    """command line interface to multiple other command line scripts"""
    dev = True
    msg = 'Usage:\n'
    msg += '  bdf merge         (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    msg += '  bdf equivalence   IN_BDF_FILENAME EQ_TOL\n'
    msg += '  bdf renumber      IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    msg += '  bdf mirror        IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
    msg += '  bdf export_mcids  IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x] [--no_y]\n'
    msg += '  bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
    if dev:
        msg += '  bdf bin          IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n'
    msg += '\n'
    msg += '  bdf merge         -h | --help\n'
    msg += '  bdf equivalence   -h | --help\n'
    msg += '  bdf renumber      -h | --help\n'
    msg += '  bdf mirror        -h | --help\n'
    msg += '  bdf export_mcids  -h | --help\n'
    msg += '  bdf split_cbars_by_pin_flags  -h | --help\n'

    if dev:
        msg += '  bdf bin          -h | --help\n'
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
        cmd_line_export_mcid()
    elif sys.argv[1] == 'split_cbars_by_pin_flags':
        cmd_line_split_cbars_by_pin_flag()
    elif sys.argv[1] == 'bin' and dev:
        cmd_line_bin()
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':  # pragma: no cover
    sys.argv = sys.argv[1:]
    cmd_line()

