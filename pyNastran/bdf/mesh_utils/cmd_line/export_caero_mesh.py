import os.path
import sys
import argparse
import pyNastran
from .utils import filter_no_args


def cmd_line_export_caero_mesh(argv=None, quiet=False):
    """command line interface to export_caero_mesh"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf export_caero_mesh: use 'bdf export_caero_mesh -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf export_caero_mesh',
        description='Export CAERO mesh to a BDF file',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_CAERO_BDF_FILENAME',
                        help='path to output BDF file')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('-x', '--xref', action='store_true',
                        help='flag to disable xref (default=False)')
    parser.add_argument('--aerobox', action='store_true',
                        help='write the aeroboxes (default=False)')
    parser.add_argument('--pid', default='aesurf', metavar='PID',
                        help="sets the pid; {aesurf, caero, paero} (default=aesurf)")
    parser.add_argument('--skip_zero_check', action='store_true',
                        help='flag to skip W2GJ, WKK, etc. checks (default=False)')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    xref = not args.xref
    caero_bdf_filename = args.output
    base = os.path.splitext(bdf_filename)[0]
    if caero_bdf_filename is None:
        caero_bdf_filename = base + '.caero.bdf'
    is_aerobox_model = args.aerobox
    skip_zero_check = args.skip_zero_check

    pid_method = args.pid

    # from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.aero.export_caero_mesh import export_caero_mesh
    # skip_cards = [
    #     # elements
    #     'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CONM2',
    #     'CROD', 'CTUBE', 'CONROD', 'CBAR', 'CBEAM',
    #     'CQUAD4', 'CTRIA3',
    #     'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM',
    #     'RBE1', 'RBE2', 'RBE3', 'RBAR',
    #
    #     # properties
    #     'PELAS', 'PDAMP', 'PROD', 'PTUBE',
    #     'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
    #     'PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG', 'PSOLID',
    #     'MAT1', 'MAT8',
    #
    #     # loads
    #     'PLOAD', 'PLOAD2', 'PLOAD4',
    #     'FORCE', 'FORCE1', 'FORCE2',
    #     'MOMENT', 'MOMENT1', 'MOMENT2',
    #     'GRAV', 'ACCEL', 'ACCEL1',
    #     # constraints
    #     'SPC', 'SPC1', 'SPCAX', 'SPCADD', 'DEQATN',
    #     'MPC', 'MPCAX', 'MPCADD',
    #     'NSM', 'NSM1', 'NSML', 'NSML1', 'NSMADD',
    #
    #     #  optimization
    #     'DVPREL1', 'DVPREL2', 'DVMREL1', 'DVMREL2', 'DVCREL1', 'DVCREL2', 'DCONADD',
    #     'DRESP1', 'DRESP2', 'DRESP3', 'DESVAR', 'DCONSTR',
    #     #  aero: maybe enable later
    #     'TRIM', 'AESTAT', 'FLUTTER', 'FLFACT',
    # ]
    # level = 'debug' if not quiet else 'warning'
    # log = SimpleLogger(level=level, encoding='utf-8')
    # model = read_bdf(bdf_filename, punch=punch, xref=xref,
    #                  log=log, skip_cards=skip_cards)

    export_caero_mesh(bdf_filename, caero_bdf_filename,
                      is_aerobox_model=is_aerobox_model,
                      pid_method=pid_method,
                      xref=xref,
                      skip_zero_check=skip_zero_check)
