from __future__ import annotations
import sys
from cpylog import SimpleLogger


def cmd_line_convert(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_merge"""
    if argv is None:  # pragma: no cover
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

    import pyNastran
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
