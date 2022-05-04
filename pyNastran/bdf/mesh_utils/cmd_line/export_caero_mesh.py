import sys
from cpylog import SimpleLogger
import pyNastran
from .utils import filter_no_args

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
        #  aero: maybe enable later
        'TRIM', 'AESTAT', 'FLUTTER', 'FLFACT',
    ]
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = read_bdf(bdf_filename, log=log, skip_cards=skip_cards)
    export_caero_mesh(model, caero_bdf_filename,
                      is_subpanel_model=is_subpanel_model, pid_method=pid_method)
