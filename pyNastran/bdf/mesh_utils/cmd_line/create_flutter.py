#import os
import sys
#from io import StringIO
#from typing import Dict, Any
from cpylog import SimpleLogger
#import pyNastran
from .utils import filter_no_args


def cmd_line_create_flutter(argv=None, quiet: bool=False):
    """command line interface to flip_shell_normals"""
    if argv is None:
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        # SWEEP_UNIT
        '  bdf flutter UNITS eas  EAS1  EAS2  N CONST_TYPE CONST_VAL [-o OUT_BDF_FILENAME] [--size SIZE] [--clean]\n'
        '  bdf flutter UNITS tas  TAS1  TAS2  N CONST_TYPE CONST_VAL [--eas_limit EAS EAS_UNITS] [-o OUT_BDF_FILENAME] [--size SIZE | --clean]\n'
        '  bdf flutter UNITS alt  ALT1  ALT2  N CONST_TYPE CONST_VAL [--eas_limit EAS EAS_UNITS] [-o OUT_BDF_FILENAME] [--size SIZE | --clean]\n'
        '  bdf flutter UNITS mach MACH1 MACH2 N CONST_TYPE CONST_VAL [--eas_limit EAS EAS_UNITS] [-o OUT_BDF_FILENAME] [--size SIZE | --clean]\n'
        '  bdf flutter -h | --help\n'
        '  bdf flutter -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  ALT, ALT1, ALT2     altitude;            sweep_units = [m, ft]\n'
        '  EAS1, EAS2          equivalent airspeed; sweep_units = [m/s, ft/s, ft/s, knots]\n'
        '  TAS1, EAS2          true airspeed;       sweep_units = [m/s, ft/s, ft/s, knots]\n'
        '  MACH1, MACH2        mach number\n'
        '  SWEEP UNIT          the unit for sweeping across'
        "  N                   the number of points in the sweep\n"
        '  CONST_TYPE          the parameter to be held constant when sweeping (alt, mach)'
        '  CONST_VAL           the value corresponding to CONST_TYPE'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file (default=flutter_cards.inc)\n'
        ' --size SIZE                          size of the BDF (8/16; default=16)\n'
        ' --clean                              writes a BDF with at least 1 whitespace in an FLFACT field (for readability)\n'
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
        '\n'
        'Examples:\n'
        '  bdf flutter english_in mach .05 0.5 101 alt 2500\n'
        '  bdf flutter english_in mach .05 0.5 101 alt 2500 --eas_limit 300 knots --out flutter_cards_temp.inc --size 16\n'
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    cmd = ' '.join(argv[1:])
    data = docopt(msg, version=ver, argv=argv[1:])

    if not quiet:  # pragma: no cover
        print(data)

    import numpy as np
    from pyNastran.bdf.bdf import BDF

    size = 16
    if data['--size']:
        size = int(data['--size'])
    units = data['UNITS']
    npoints = int(data['N'])
    clean = data['--clean']

    assert clean in [True, False], clean

    const_type = data['CONST_TYPE'].lower()
    assert const_type in {'alt', 'mach', 'eas', 'tas'}, const_type
    const_value = float(data['CONST_VAL'])
    if const_type == 'alt':
        alt = const_value
    elif const_type == 'mach':
        mach = const_value
    else:
        raise NotImplementedError(const_type)

    eas_units = ''
    eas_limit = 1_000_000.
    if data['--eas_limit']:
        eas_limit = float(data['EAS'])
        eas_units = data['EAS_UNITS']
        assert eas_units not in [None, ''], eas_units
        eas_units = eas_units.lower()
        assert eas_units in {'m/s', 'cm/s', 'in/s', 'ft/s', 'knots'}

    sweep_method = ''
    if data['alt']:
        sweep_method = 'alt'
        alt1 = float(data['ALT1'])
        alt2 = float(data['ALT2'])
        alts = np.linspace(alt1, alt2, num=npoints)
        #sweep_unit = data['SWEEP_UNIT'].lower()
        #assert sweep_unit in ['m', 'ft'], f'sweep_unit={sweep_unit}'

    elif data['mach']:
        sweep_method = 'mach'
        mach1 = float(data['MACH1'])
        mach2 = float(data['MACH2'])
        machs = np.linspace(mach1, mach2, num=npoints)

    elif data['eas']:
        sweep_method = 'eas'
        eas1 = float(data['EAS1'])
        eas2 = float(data['EAS2'])
        eass = np.linspace(eas1, eas2, num=npoints)
        #sweep_unit = data['SWEEP_UNIT'].lower()
        #assert sweep_unit in ['m/s', 'ft/s', 'in/s', 'knots'], f'sweep_unit={sweep_unit}'

    elif data['tas']:
        sweep_method = 'tas'
        tas1 = float(data['TAS1'])
        tas2 = float(data['TAS2'])
        tass = np.linspace(tas1, tas2, num=npoints)
        #sweep_unit = data['SWEEP_UNIT'].lower()
        #assert sweep_unit in ['m/s', 'ft/s', 'in/s', 'knots'], f'sweep_unit={sweep_unit}'
    else:
        raise NotImplementedError(data)

    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'flutter_cards.inc'

    #from io import StringIO
    #from pyNastran.bdf.cards.aero.dynamic_loads import
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8', log_func=None)
    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    flutter_method = 'PKNL'
    sid = 1

    flfact_density = sid + 1
    flfact_mach = sid + 2
    flfact_velocity = sid + 3
    #flfact_eas = sid + 4

    flutter = model.add_flutter(
        sid, flutter_method, flfact_density, flfact_mach, flfact_velocity,
        imethod='L', nvalue=None, omax=None, epsilon=1.0e-3,
        comment='', validate=True)

    units_map = {
        # (alt, velocity, density, eas, pressure)
        'english_in': ('ft', 'in/s', 'slinch/in^3', 'knots', 'psi'),
        'english_ft': ('ft', 'ft/s', 'slug/ft^3', 'knots', 'psf'),
        'si': ('m', 'm/s', 'kg/m^3', 'knots', 'Pa'),
    }
    try:
        unitsi = units_map[units.lower()]
    except KeyError:
        raise NotImplementedError(units)
    alt_units, velocity_units, density_units, eas_units_default, pressure_units = unitsi

    if eas_units in [None, '']:
        eas_units = eas_units_default

    log.info(f'alt_units={alt_units!r}')
    log.info(f'velocity_units={velocity_units!r}')
    log.info(f'density_units={density_units!r}')
    log.info(f'eas_units={eas_units!r}')

    # option 1: overwrite the eas/tas/alt unit...would work for alt
    #           we'll overwrite the
    # option 2: pass another flag in...I don't wanna
    pairs = [
        ('eas', 'alt'),
        ('mach', 'alt'),
        ('alt', 'mach'),
        ('tas', 'alt'),
        ('eas', 'mach'),
    ]
    assert alt_units != '', alt_units
    assert velocity_units != '', velocity_units
    assert density_units != '', density_units
    assert eas_units != '', eas_units
    assert pressure_units != '', pressure_units

    if sweep_method == 'eas' and const_type == 'alt':
        #eas_units = sweep_unit
        flutter.make_flfacts_eas_sweep_constant_alt(
            model, alt, eass,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    elif sweep_method == 'eas' and const_type == 'mach':
        gamma = 1.4
        flutter.make_flfacts_eas_sweep_constant_mach(
            machs, eass,
            gamma=gamma, alt_units=alt_units,
            density_units=density_units,
            pressure_units=pressure_units,
            eas_units=eas_units)

    elif sweep_method == 'mach' and const_type == 'alt':
        flutter.make_flfacts_mach_sweep_constant_alt(
            model, alt, machs,
            eas_limit=eas_limit,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    elif sweep_method == 'alt' and const_type == 'mach':
        #alt_units = sweep_unit
        flutter.make_flfacts_alt_sweep_constant_mach(
            model, mach, alts,
            eas_limit=eas_limit,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    elif sweep_method == 'tas' and const_type == 'alt':
        #velocity_units = sweep_unit
        flutter.make_flfacts_tas_sweep_constant_alt(
            model, alt, tass,
            alt_units=alt_units,
            eas_limit=eas_limit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    else:
        raise NotImplementedError((sweep_method, const_type))

    model.punch = True
    flutter.comment = cmd
    if clean:
        # makes a "clean" deck by writing the data in small field
        # we take advantage of truncation to get a more readable deck
        #
        # the downsides are we have to write twice and we lose extra precision
        model.write_bdf(bdf_filename_out, encoding=None, size=8,
                        nodes_size=None, elements_size=None, loads_size=None,
                        is_double=False, interspersed=False, enddata=None, write_header=True, close=True)

        model2 = BDF(log=log)
        model2.read_bdf(bdf_filename_out, validate=True, xref=False, punch=True, read_includes=True,
                       save_file_structure=False, encoding=None)
        model2.write_bdf(bdf_filename_out, encoding=None, size=16,
                        nodes_size=None, elements_size=None, loads_size=None,
                        is_double=False, interspersed=False, enddata=None, write_header=True, close=True)
    else:
        model.write_bdf(bdf_filename_out, encoding=None, size=size,
                        nodes_size=None, elements_size=None, loads_size=None,
                        is_double=False, interspersed=False, enddata=None, write_header=True, close=True)
    print(cmd)
