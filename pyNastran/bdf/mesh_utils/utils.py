from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
from pyNastran.bdf.mesh_utils.delete_bad_elements import delete_bad_shells, get_bad_shells


def cmd_line_equivalence():
    import sys
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf equivalence IN_BDF_FILENAME EQ_TOL\n"

    msg += '  bdf equivalence -h | --help\n'
    msg += '  bdf equivalence -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAMES   path to input BDF/DAT/NAS files\n"
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
    asdf
    size = 16
    bdf_filenames = data['IN_BDF_FILENAMES']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'merged.bdf'

    cards_to_skip = ['AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2', 'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    aaa
    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip)

def cmd_line_bin():
    import sys
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
    xyz_cid = model.get_xyz_in_coord(cid=cid, dtype='float64')
    y = xyz_cid[:, axis1]
    z = xyz_cid[:, axis2]

    plt.figure(1)
    #n, bins, patches = plt.hist( [x0,x1,x2], 10, weights=[w0, w1, w2], histtype='bar')
    ys = []
    zs = []
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



def cmd_line_renumber():
    import sys
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
    raise NotImplementedError(data)
    size = 16
    bdf_filenames = data['IN_BDF_FILENAMES']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'renumber.bdf'

    cards_to_skip = ['AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2', 'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    aaa
    bdf_renumber(bdf_filename, bdf_filename_out, size=size, is_double=False,
                starting_id_dict=None, round_ids=False,
                cards_to_skip=cards_to_skip)

def cmd_line_merge():
    import sys
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf merge (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-d] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-l] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C]      [-p] [-r] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-s] [-f] [--encoding ENCODE] BDF_FILENAME\n"

    #msg += "  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n" #
    msg += '  bdf merge -h | --help\n'
    msg += '  bdf merge -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  IN_BDF_FILENAMES   path to input BDF/DAT/NAS files\n"
    msg += "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
    msg += "  EQ_TOL             spherical equivalence tolerance\n"
    #msg += "  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
    msg += '\n'

    msg += 'Options:\n'
    #msg += '  --crash C,     Crash on specific cards (e.g. CGEN,EGRID)\n'
    #msg += '  -q, --quiet    prints debug messages (default=False)\n'
    #msg += '  -x, --xref     disables cross-referencing and checks of the BDF.\n'
    #msg += '                 (default=True -> on)\n'
    #msg += '  -p, --punch    disables reading the executive and case control decks in the BDF\n'
    #msg += '                 (default=False -> reads entire deck)\n'
    #msg += '  -c, --check    disables BDF checks.  Checks run the methods on \n'
    #msg += '                 every element/property to test them.  May fails if a \n'
    #msg += '                 card is fully not supported (default=False)\n'
    #msg += '  -l, --large    writes the BDF in large field, single precision format (default=False)\n'
    #msg += '  -d, --double   writes the BDF in large field, double precision format (default=False)\n'
    #msg += '  -L, --loads    Disables forces/moments summation for the different subcases (default=True)\n'
    #msg += '  -r, --reject   rejects all cards with the appropriate values applied (default=False)\n'
    #msg += '  -D, --dumplines  Writes the BDF exactly as read with the INCLUDES processed (pyNastran_dump.bdf)\n'
    #msg += '  -i, --dictsort  Writes the BDF with exactly as read with the INCLUDES processed (pyNastran_dict.bdf)\n'
    #msg += '  -f, --profile   Profiles the code (default=False)\n'
    #msg += '  -s, --stop      Stop after first read/write (default=False)\n'
    #msg += '  -e E, --nerrors E  Allow for cross-reference errors (default=100)\n'
    #msg += '  --encoding ENCODE  the encoding method\n'
    #msg += '  -o <VAR_VAL>, --openmdao <VAR_VAL>   rejects all cards with the appropriate values applied;\n'
    #msg += '                 Uses the OpenMDAO %var syntax to replace it with value.\n'
    #msg += '                 So test_bdf -r var1=val1 var2=val2\n'

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

    cards_to_skip = ['AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2', 'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip)

def cmd_line():
    import sys
    dev = True
    msg = 'Usage:\n'
    msg += '  bdf merge       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n'
    msg += '  bdf equivalence IN_BDF_FILENAME EQ_TOL\n'
    msg += '  bdf renumber    IN_BDF_FILENAME [-o OUT_BDF_FILENAME]\n'
    if dev:
        msg += '  bdf bin         IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]\n'
    msg += '\n'
    msg += '  bdf merge       -h | --help\n'
    msg += '  bdf equivalence -h | --help\n'
    msg += '  bdf renumber    -h | --help\n'

    if dev:
        msg += '  bdf bin         -h | --help\n'
    msg += '  bdf -v | --version\n'
    msg += '\n'

    if len(sys.argv) == 1:
        sys.exit(msg)

    #assert sys.argv[0] != 'bdf', msg

    if sys.argv[1] == 'merge':
        cmd_line_merge()
    elif sys.argv[1] == 'equivalence':
        cmd_line_equivalence()
    elif sys.argv[1] == 'equivalence':
        cmd_line_renumber()
    elif sys.argv[1] == 'bin' and dev:
        cmd_line_bin()
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':
    main()

