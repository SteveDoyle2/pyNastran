from __future__ import annotations
import sys
import argparse
from cpylog import SimpleLogger


def cmd_line_convert(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_merge"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    if len(argv) == 1:
        sys.exit("bdf convert: use 'bdf convert -h' for help")

    import pyNastran
    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf convert',
        description='Convert BDF units',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  bdf convert model.bdf --in_units m,kg  --out_units in,lbm\n'
            '  bdf convert model.bdf --in_units m,kg  --out_units in,slinch\n'
            '  bdf convert model.bdf --in_units m,kg  --out_units ft,slug\n'
            '  bdf convert model.bdf --in_units m,kg  --out_units ft,lbm\n'
        ),
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file')
    parser.add_argument('--in_units', default='m,kg', metavar='IN_UNITS',
                        help='length,mass (default=m,kg)')
    parser.add_argument('--out_units', default='m,kg', metavar='OUT_UNITS',
                        help='length,mass (default=m,kg)')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    bdf_filename_out = args.output
    if bdf_filename_out is None:
        bdf_filename_out = bdf_filename + '.convert.bdf'

    in_units = args.in_units
    out_units = args.out_units

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

    # cards_to_skip = [
    #     'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
    #     'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.convert import convert

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    model = read_bdf(bdf_filename, validate=True, xref=True,
                     punch=False, save_file_structure=False,
                     skip_cards=None, read_cards=None,
                     encoding=None, log=log, debug=True, mode='msc')
    convert(model, units_to, units=units)
    for prop in model.properties.values():
        prop.comment = ''
    model.write_bdf(bdf_filename_out)
