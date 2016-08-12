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


def cmd_line_renumber():
    import sys
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  bdf renumber IN_BDF_FILENAME\n"
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
        bdf_filename_out = 'merged.bdf'

    cards_to_skip = ['AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2', 'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    aaa
    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip)

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
    msg = "Usage:\n"
    msg += "  bdf merge       (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]\n"
    msg += "  bdf equivalence IN_BDF_FILENAME EQ_TOL\n"
    msg += "  bdf renumber    IN_BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-d] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-c] [-L] [-l] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C]      [-p] [-r] [-f] [--encoding ENCODE] BDF_FILENAME\n"
    #msg += "  test_bdf [-q] [-D] [-i] [-e E] [--crash C] [-x] [-p] [-s] [-f] [--encoding ENCODE] BDF_FILENAME\n"

    #msg += "  test_bdf [-q] [-p] [-o [<VAR=VAL>]...] BDF_FILENAME\n" #
    msg += '\n'
    msg += '  bdf merge       -h | --help\n'
    msg += '  bdf equivalence -h | --help\n'
    msg += '  bdf renumber    -h | --help\n'
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
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':
    main()

