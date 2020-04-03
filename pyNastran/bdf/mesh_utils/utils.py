"""
defines:
    bdf merge        (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    bdf equivalence  IN_BDF_FILENAME EQ_TOL\n'
    bdf renumber     IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    bdf mirror       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
    bdf export_mcids IN_BDF_FILENAME [-o OUT_GEOM_FILENAME]\n'
    bdf split_cbars_by_pin_flags IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'

"""
import os
import sys
from cpylog import SimpleLogger
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber, superelement_renumber
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.bdf.mesh_utils.export_mcids import export_mcids
from pyNastran.bdf.mesh_utils.pierce_shells import pierce_shell_model

# testing these imports are up to date
# if something is imported and tested, it should be removed from here
from pyNastran.bdf.mesh_utils.shift import update_nodes
from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric
from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells
from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused
from pyNastran.bdf.mesh_utils.free_faces import write_skin_solid_faces


def cmd_line_equivalence(argv=None, quiet=False):
    """command line interface to bdf_equivalence_nodes"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
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
        bdf_filename_out = 'merged.bdf'
    tol = data['EQ_TOL']
    size = 16
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


def cmd_line_bin(argv=None, quiet=False):  # pragma: no cover
    """bins the model into nbins"""
    if argv is None:
        argv = sys.argv

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
    import pyNastran
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
        size = int(data['SIZE'])

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


def cmd_line_mirror(argv=None, quiet=False):
    """command line interface to write_bdf_symmetric"""
    if argv is None:
        argv = sys.argv

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
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        "  --plane PLANE                       the symmetry plane (xz, yz, xy); default=xz\n"
        '  --tol   TOL                         the spherical equivalence tolerance; default=1e-6\n'
        '  --noeq                              disable equivalencing\n'
        "\n" #  (default=0.000001)

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
    if data['--tol'] is None:
        data['TOL'] = 0.000001

    tol = data['TOL']
    if data['--noeq'] is not None:
        tol = -1.

    plane = 'xz'
    if data['--plane'] is not None:
        plane = data['--plane']

    if not quiet:  # pragma: no cover
        print(data)
    size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'mirrored.bdf'

    #from io import StringIO
    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log)
    size = 16
    bdf_filename_temp = '__temp__.bdf'
    write_bdf_symmetric(model, bdf_filename_temp, encoding=None, size=size,
                        is_double=False,
                        enddata=None, close=True,
                        plane=plane, log=log)

    bdf_equivalence_nodes(bdf_filename_temp, bdf_filename_out, tol,
                          renumber_nodes=False,
                          neq_max=10, xref=True,
                          node_set=None, size=size,
                          is_double=False,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False,
                          debug=True, log=log)
    os.remove(bdf_filename_temp)


def cmd_line_merge(argv=None, quiet=False):
    """command line interface to bdf_merge"""
    if argv is None:
        argv = sys.argv

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
    if len(argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
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


def cmd_line_convert(argv=None, quiet=False):
    """command line interface to bdf_merge"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        '  bdf convert IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
        '  bdf convert -h | --help\n'
        '  bdf convert -v | --version\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n'
        '  --in_units  IN_UNITS                length,mass\n\n'
        '  --out_units  OUT_UNITS              length,mass\n\n'

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
    import pyNastran
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
    import pyNastran
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


def cmd_line_free_faces(argv=None, quiet=False):
    """command line interface to bdf free_faces"""
    if argv is None:
        argv = sys.argv

    encoding = sys.getdefaultencoding()
    import pyNastran
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
        '  --encoding ENCODE  the encoding method (default=None -> %r)\n'
        '\n'
        'Developer:\n'
        '  -f, --profile    Profiles the code (default=False)\n'
        '\n'
        "Info:\n"
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n" % encoding
    )
    if len(argv) == 1:
        sys.exit(arg_msg)

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
        print("total time:  %.2f sec" % (time.time() - time0))


def cmd_line_split_cbars_by_pin_flag(argv=None, quiet=False):
    """command line interface to split_cbars_by_pin_flag"""
    if argv is None:
        argv = sys.argv

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
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf filter IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME]\n'
        '  bdf filter IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--x YSIGN_X] [--y YSIGN_Y] [--z YSIGN_Z]\n'
        '  bdf filter -h | --help\n'
        '  bdf filter -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME         path to output BDF file (default=filter.bdf)\n'
        " --x YSIGN_X                                a string (e.g., '< 0.')\n"
        " --y YSIGN_Y                                a string (e.g., '< 0.')\n"
        " --z YSIGN_Z                                a string (e.g., '< 0.')\n"
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


def cmd_line_export_caero_mesh(argv=None, quiet=False):
    """command line interface to export_caero_mesh"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  bdf export_caero_mesh IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--subpanels] [--pid PID]\n'
        '  bdf export_caero_mesh -h | --help\n'
        '  bdf export_caero_mesh -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_CAERO_BDF_FILENAME  path to output BDF file\n'
        '  --subpanels                               write the subpanels (default=False)\n'
        '  --pid PID                                 sets the pid; {aesurf, caero, paero} [default: aesurf]\n'
        '\n'

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
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    caero_bdf_filename = data['--output']
    if caero_bdf_filename is None:
        caero_bdf_filename = 'caero.bdf'
    is_subpanel_model = data['--subpanels']

    pid_method = 'aesurf'
    if data['--pid']:
        pid_method = data['--pid']

    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.export_caero_mesh import export_caero_mesh
    skip_cards = [
        # elements
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CONM2',
        'CROD', 'CTUBE', 'CONROD', 'CBAR', 'CBEAM',
        'CQUAD4', 'CTRIA3',
        'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
        'RBE1', 'RBE2', 'RBE3', 'RBAR',

        # properties
        'PELAS', 'PDAMP', 'PROD', 'PTUBE',
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG', 'PSOLID',
        'MAT1', 'MAT8',

        # loads
        'PLOAD', 'PLOAD2', 'PLOAD4', 'FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1', 'MOMENT2',
        'GRAV', 'ACCEL', 'ACCEL1',
        # constraints
        'SPC', 'SPC1', 'MPC', 'SPCADD', 'MPCADD', 'DEQATN',

        #  optimization
        'DVPREL1', 'DVPREL2', 'DVMREL1', 'DVMREL2', 'DVCREL1', 'DVCREL2', 'DCONADD',
        'DRESP1', 'DRESP2', 'DRESP3', 'DESVAR',
        #  aero: mabye enable later
        'TRIM', 'AESTAT', 'FLUTTER', 'FLFACT',
    ]
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log, skip_cards=skip_cards)
    export_caero_mesh(model, caero_bdf_filename,
                      is_subpanel_model=is_subpanel_model, pid_method=pid_method)

def cmd_line(argv=None, quiet=False):
    """command line interface to multiple other command line scripts"""
    if argv is None:
        argv = sys.argv

    dev = True
    msg = (
        'Usage:\n'
        '  bdf merge                       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
        '  bdf equivalence                 IN_BDF_FILENAME EQ_TOL\n'
        '  bdf renumber                    IN_BDF_FILENAME [OUT_BDF_FILENAME] [--superelement] [--size SIZE]\n'
        '  bdf mirror                      IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]\n'
        '  bdf convert                     IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--in_units IN_UNITS] [--out_units OUT_UNITS]\n'
        '  bdf scale                       IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--lsf LENGTH_SF] [--msf MASS_SF] [--fsf FORCE_SF] [--psf PRESSURE_SF] [--tsf TIME_SF] [--vsf VEL_SF]\n'
        '  bdf export_mcids                IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x | --no_y]\n'
        '  bdf free_faces                  BDF_FILENAME SKIN_FILENAME [-d | -l] [-f] [--encoding ENCODE]\n'
        '  bdf transform                   IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--shift XYZ]\n'
        '  bdf export_caero_mesh           IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--subpanels] [--pid PID]\n'
        '  bdf split_cbars_by_pin_flags    IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]\n'
    )

    if dev:
        msg += '  bdf filter                      IN_BDF_FILENAME [-o OUT_CAERO_BDF_FILENAME] [--x YSIGN X] [--y YSIGN Y] [--z YSIGN Z]\n'
        msg += '  bdf bin                         IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n'

    msg += (
        #'\n'
        '  bdf merge              -h | --help\n'
        '  bdf equivalence        -h | --help\n'
        '  bdf renumber           -h | --help\n'
        '  bdf mirror             -h | --help\n'
        '  bdf convert            -h | --help\n'
        '  bdf scale              -h | --help\n'
        '  bdf export_mcids       -h | --help\n'
        '  bdf free_faces         -h | --help\n'
        '  bdf transform          -h | --help\n'
        '  bdf filter             -h | --help\n'
        '  bdf export_caero_mesh  -h | --help\n'
        '  bdf split_cbars_by_pin_flags  -h | --help\n'
    )

    if dev:
        msg += '  bdf filter                      -h | --help\n'
        msg += '  bdf bin                         -h | --help\n'
    msg += '  bdf -v | --version\n'
    msg += '\n'

    if len(argv) == 1:
        sys.exit(msg + 'Not enough arguments.\n')

    #assert sys.argv[0] != 'bdf', msg

    if argv[1] == 'merge':
        cmd_line_merge(argv, quiet=quiet)
    elif argv[1] == 'equivalence':
        cmd_line_equivalence(argv, quiet=quiet)
    elif argv[1] == 'renumber':
        cmd_line_renumber(argv, quiet=quiet)
    elif argv[1] == 'mirror':
        cmd_line_mirror(argv, quiet=quiet)
    elif argv[1] == 'convert':
        cmd_line_convert(argv, quiet=quiet)
    elif argv[1] == 'scale':
        cmd_line_scale(argv, quiet=quiet)
    elif argv[1] == 'export_mcids':
        cmd_line_export_mcids(argv, quiet=quiet)
    elif argv[1] == 'split_cbars_by_pin_flags':
        cmd_line_split_cbars_by_pin_flag(argv, quiet=quiet)
    elif argv[1] == 'export_caero_mesh':
        cmd_line_export_caero_mesh(argv, quiet=quiet)
    elif argv[1] == 'transform':
        cmd_line_transform(argv, quiet=quiet)
    elif argv[1] == 'filter' and dev:  # TODO: make better name
        cmd_line_filter(argv, quiet=quiet)
    elif argv[1] == 'free_faces':
        cmd_line_free_faces(argv, quiet=quiet)
    elif argv[1] == 'bin' and dev:
        cmd_line_bin(argv, quiet=quiet)
    elif argv[1] in ['-v', '--version']:
        import pyNastran
        print(pyNastran.__version__)
    else:
        print(argv)
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':  # pragma: no cover
    # for the exe, we pass all the args, but we hack them to have the bdf prefix
    from copy import deepcopy
    argv = deepcopy(sys.argv)
    argv[0] = 'bdf'
    cmd_line(argv=argv)
