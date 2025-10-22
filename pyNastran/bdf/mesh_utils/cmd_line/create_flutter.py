from __future__ import annotations
import sys
from typing import Any, TYPE_CHECKING
import numpy as np

from cpylog import SimpleLogger
from .utils import filter_no_args
from pyNastran.utils.convert import convert_altitude, convert_velocity
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF

UNITS_MAP = {
    # units must be consistent
    #
    # (alt, velocity, density, eas, pressure)
    'english_in': ('ft', 'in/s', 'slinch/in^3', 'knots', 'psi'),
    'english_ft': ('ft', 'ft/s', 'slug/ft^3', 'knots', 'psf'),
    'si': ('m', 'm/s', 'kg/m^3', 'm/s', 'Pa'),
    'si_mm': ('m', 'mm/s', 'Mg/mm^3', 'mm/s', 'MPa'),
}

ALT_UNITS = ['m', 'ft', 'kft']
VELOCITY_UNITS = ['m/s', 'cm/s', 'mm/s', 'in/s', 'ft/s', 'knots']


def cmd_line_create_flutter(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf flutter"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    options = '[-o OUT_BDF_FILENAME] [--size SIZE | --clean] [--sid SID] [--rhoref]'
    msg = (
        'Usage:\n'
        # SWEEP_UNIT
        # CONST_TYPEs = [mach, alt]

        # CONST_TYPE = mach
        #'  bdf flutter gui\n'
        f'  bdf flutter UNITS eas  EAS1  EAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT {options}\n'
        f'  bdf flutter UNITS tas  TAS1  TAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS] {options}\n'
        f'  bdf flutter UNITS alt  ALT1  ALT2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS] {options}\n'

        # CONST_TYPE = alt
        #f'  bdf flutter UNITS eas  EAS1  EAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT {options}\n'
        #f'  bdf flutter UNITS tas  TAS1  TAS2  SWEEP_UNIT N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS] {options}\n'
        f'  bdf flutter UNITS mach MACH1 MACH2            N CONST_TYPE CONST_VAL CONST_UNIT [--eas_limit EAS EAS_UNITS] {options}\n'

        '  bdf flutter -h | --help\n'
        '  bdf flutter -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        #'  gui                 enables the gui\n'
        '  alt, ALT1, ALT2     altitude;            units = [m, ft, kft]\n'
        '  eas, EAS1, EAS2     equivalent airspeed; units = [m/s, cm/s, in/s, ft/s, knots]\n'
        '  tas, TAS1, EAS2     true airspeed;       units = [m/s, cm/s, in/s, ft/s, knots]\n'
        '  mach, MACH1, MACH2  mach number;         units = [none, na]\n'
        '  SWEEP_UNIT          the unit for sweeping across\n'
        "  N                   the number of points in the sweep\n"
        '  alt, mach           the parameter to be held constant when sweeping (alt, mach)\n'
        '  CONST_VAL           the value corresponding to CONST_TYPE\n'
        "  CONST_UNIT          the unit for the altitude that is held constant\n"
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file (default=flutter_cards.inc)\n'
        ' --size SIZE                          size of the BDF (8/16; default=16)\n'
        ' --clean                              writes a BDF with at least 1 whitespace in an FLFACT field (for readability)\n'
        ' --sid SID                            updates the flutter ID\n'
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
        '\n'
        'Examples:\n'
        '  bdf flutter english_in eas  1   800. knots 101 mach 0.8 na\n'
        '  bdf flutter english_in tas  .1  800. ft/s  101 alt 2500 m\n'
        '  bdf flutter english_in mach .05 0.5        101 alt 2500\n'
        '  bdf flutter english_in mach .05 0.5        101 alt 2500 m --eas_limit 300 knots --out flutter_cards_temp.inc --size 16\n'
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    if 'gui' in argv:
        from pyNastran.bdf.mesh_utils.gui_tools.gui_flutter import cmd_line_gui
        data = cmd_line_gui()
    else:
        argv = [str(arg) for arg in argv]
        cmd = 'bdf ' + ' '.join(argv[1:])
        data = docopt(msg, version=ver, argv=argv[1:])

    if not quiet:  # pragma: no cover
        print(data)

    size = 16
    if data['--size']:
        size = _int(data, '--size')

    sid = 1
    if data['--sid']:
        sid = _int(data, '--sid')

    units_out = data['UNITS']
    if units_out.lower() not in UNITS_MAP:  # pragma: no cover
        raise NotImplementedError(units_out)

    rhoref_flag = data['--rhoref']
    npoints = _int(data, 'N')
    clean = data['--clean']
    assert clean in [True, False], clean

    const_type = data['CONST_TYPE'].lower()
    assert const_type in {'alt', 'mach', 'eas', 'tas'}, f'const_type={const_type!r}'
    const_value = _float(data, 'CONST_VAL')
    const_unit = data['CONST_UNIT'].lower()

    eas_units = ''
    eas_limit = 1_000_000.
    if data['--eas_limit']:
        eas_limit = _float(data, 'EAS')
        eas_units = data['EAS_UNITS']
        assert eas_units not in {None, ''}, eas_units
        eas_units = eas_units.lower()
        #assert eas_units in VELOCITY_UNITS, f'eas_unit={eas_unit!r}; allowed={VELOCITY_UNITS}'

    sweep_unit = ''
    if data['alt']:
        sweep_method = 'alt'
        value1 = _float(data, 'ALT1')
        value2 = _float(data, 'ALT2')
        sweep_unit = data['SWEEP_UNIT'].lower()
    elif data['mach']:
        sweep_method = 'mach'
        value1 = _float(data, 'MACH1')
        value2 = _float(data, 'MACH2')
    elif data['eas']:
        sweep_method = 'eas'
        value1 = _float(data, 'EAS1')
        value2 = _float(data, 'EAS2')
        sweep_unit = data['SWEEP_UNIT'].lower()
    elif data['tas']:
        sweep_method = 'tas'
        value1 = _float(data, 'TAS1')
        value2 = _float(data, 'TAS2')
        sweep_unit = data['SWEEP_UNIT'].lower()
    else:  # pragma: no cover
        raise NotImplementedError(data)

    #alts = np.linspace(alt1, alt2, num=npoints)

    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'flutter_cards.inc'

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    create_flutter(log,
                   sweep_method, value1, value2, sweep_unit, npoints,
                   const_type, const_value, const_unit,
                   eas_limit=eas_limit, eas_units=eas_units,
                   units_out=units_out,
                   rhoref_flag=rhoref_flag,
                   sid=sid,
                   size=size, clean=clean,
                   bdf_filename_out=bdf_filename_out,
                   comment=cmd)
    if not quiet:
        print(cmd)


def create_flutter(log: SimpleLogger,
                   sweep_method: str, value1: float, value2: float, sweep_unit: str, npoints: int,
                   const_type: str, const_value: float, const_unit: str,
                   sid: int=1,
                   eas_limit: float=1_000_000.0, eas_units: str='m/s',
                   units_out: str='si',
                   rhoref_flag: bool=False,
                   size: int=8,
                   clean: bool=False,
                   bdf_filename_out: str='flutter_cards.inc',
                   comment: str='') -> tuple[BDF, str, str]:

    unitsi = UNITS_MAP[units_out.lower()]
    alt_units, velocity_units, density_units, eas_units_default, pressure_units = unitsi
    if const_type == 'alt':
        alt = const_value
        alt_unit = const_unit
        assert alt_unit in ALT_UNITS, f'alt_unit={alt_unit!r}; allowed={ALT_UNITS}'
        alt = convert_altitude(alt, alt_unit, alt_units)
    elif const_type == 'mach':
        mach = const_value
        assert const_unit in {'none', 'na'}, f'const_unit={const_unit!r}; allowed=[none]'
    elif const_type == 'eas':
        eas = const_value
        assert const_unit in VELOCITY_UNITS, f'const_unit={const_unit!r}; allowed={VELOCITY_UNITS}'
    elif const_type == 'tas':
        tas = const_value
        assert const_unit in VELOCITY_UNITS, f'const_unit={const_unit!r}; allowed={VELOCITY_UNITS}'
    else:  # pragma: no cover
        raise NotImplementedError(f'const_type={const_type} is not supported')

    sweep_unit = sweep_unit.lower()
    values = np.linspace(value1, value2, num=npoints)
    dvalue1 = values[1] - values[0]
    if sweep_method == 'alt':
        assert sweep_unit in ALT_UNITS, f'sweep_unit={sweep_unit!r}; allowed={ALT_UNITS}'
        alts = convert_altitude(values, sweep_unit, alt_units)
        unit2 = alt_units

    elif sweep_method == 'mach':
        machs = values
        unit2 = 'na'
    elif sweep_method == 'eas':
        eass = values
        eas_units = sweep_unit
        unit2 = eas_units
        #assert sweep_unit in VELOCITY_UNITS, f'sweep_unit={sweep_unit!r}; allowed={VELOCITY_UNITS}'
        #eass = convert_velocity(eass, sweep_unit, velocity_units)

    elif sweep_method == 'tas':
        assert sweep_unit in VELOCITY_UNITS, f'sweep_unit={sweep_unit!r}; allowed={VELOCITY_UNITS}'
        tass = convert_velocity(values, sweep_unit, velocity_units)
        unit2 = velocity_units
    elif sweep_method == 'alt':
        assert sweep_unit in ALT_UNITS, f'sweep_unit={sweep_unit!r}; allowed={ALT_UNITS}'
        alts = convert_altitude(values, sweep_unit, alt_units)
        unit2 = alt_units
    else:  # pragma: no cover
        raise NotImplementedError(sweep_method)
    dvalue2 = values[1] - values[0]

    #------------------------------------------------------------------
    from pyNastran.bdf.bdf import BDF
    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    flutter_method = 'PKNL'

    flfact_density = sid + 1
    flfact_mach = sid + 2
    flfact_velocity = sid + 3
    #flfact_eas = sid + 4

    flutter = model.add_flutter(
        sid, flutter_method, flfact_density, flfact_mach, flfact_velocity,
        imethod='L', nvalue=None, omax=None, epsilon=1.0e-3,
        comment='', validate=True)

    if eas_units in {None, ''}:
        log.debug(f'setting eas_units to default; eas_units={eas_units_default!r}')
        eas_units = eas_units_default

    log.info(f'sweep_method={sweep_method!r}')
    log.debug(f'  d{sweep_unit} = {dvalue1!r} {sweep_unit} = {dvalue2!r} {unit2}')
    log.debug(f'  alt_units={alt_units!r}')
    log.debug(f'  velocity_units={velocity_units!r}')
    log.debug(f'  density_units={density_units!r}')
    log.debug(f'  eas_units={eas_units!r}')
    #del alt_unit

    # option 1: overwrite the eas/tas/alt unit...would work for alt
    #           we'll overwrite the
    # option 2: pass another flag in...I don't wanna
    pairs = [
        # sweep, const
        ('eas', 'alt'),
        ('eas', 'mach'),
        ('mach', 'alt'),
        ('alt', 'mach'),
        ('tas', 'alt'),
        #('tas', 'mach'),  # undefined b/c sos doesn't map 1:1 to alt

        # new
        #('tas', 'eas'),  # kind of dumb, but should work...
        ('alt', 'tas'),
        #('alt', 'eas'), # kind of dumb
    ]

    assert alt_units != '', alt_units
    assert velocity_units != '', velocity_units
    assert density_units != '', density_units
    assert eas_units != '', eas_units
    assert pressure_units != '', pressure_units

    if sweep_method == 'eas' and const_type == 'alt':
        flutter.make_flfacts_eas_sweep_constant_alt(
            model, alt, eass,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    elif sweep_method == 'eas' and const_type == 'mach':
        gamma = 1.4
        flutter.make_flfacts_eas_sweep_constant_mach(  # TODO: need to test this; seems wrong
            model, mach, eass,
            gamma=gamma,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            pressure_units=pressure_units,
            eas_units=eas_units)

    elif sweep_method == 'mach' and const_type == 'alt':
        flutter.make_flfacts_mach_sweep_constant_alt(
            model, alt, machs,
            eas_limit=eas_limit,
            alt_units=alt_unit,
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
    elif sweep_method == 'alt' and const_type == 'tas':
        flutter.make_flfacts_alt_sweep_constant_tas(
            model, tas, alts,
            alt_units=alt_units,
            eas_limit=eas_limit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    elif sweep_method == 'tas' and const_type == 'alt':
        flutter.make_flfacts_tas_sweep_constant_alt(
            model, alt, tass,
            alt_units=alt_units,
            eas_limit=eas_limit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    #elif sweep_method == 'tas' and const_type == 'mach':
        # undefined b/c sos doesn't map 1:1 to alt
        #flutter.make_flfacts_tas_sweep_constant_mach(
            #model, mach, tass,
            #alt_units=alt_units,
            #eas_limit=eas_limit,
            #velocity_units=velocity_units,
            #density_units=density_units,
            #eas_units=eas_units)
    #elif sweep_method == 'tas' and const_type == 'eas':
        #flutter.make_flfacts_tas_sweep_constant_eas(
            #model, eas, tass,
            #alt_units=alt_units,
            #eas_limit=eas_limit,
            #velocity_units=velocity_units,
            #density_units=density_units,
            #eas_units=eas_units)
    elif sweep_method == 'alt' and const_type == 'eas':
        flutter.make_flfacts_alt_sweep_constant_eas(
            model, eas, alts,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units)
    else:  # pragma: no cover
        raise NotImplementedError((sweep_method, const_type))

    if rhoref_flag:
        from pyNastran.utils.atmosphere import atm_density
        rho0 = atm_density(alt=0., density_units=density_units)
        cref = 1.0
        velocity = 0.0
        model.add_aero(velocity, cref, rho0)
        flfact = model.flfacts[flfact_density]
        # print(flfact.get_stats())
        flfact.factors /= rho0

    model.punch = True
    flutter.comment = comment
    clean = True
    if bdf_filename_out:
        if clean:
            # makes a "clean" deck by writing the data in small field
            # we take advantage of truncation to get a more readable deck
            #
            # the downsides are we have to write twice and we lose extra precision
            model.write_bdf(bdf_filename_out, encoding=None, size=8,
                            nodes_size=None, elements_size=None, loads_size=None,
                            is_double=False, interspersed=False, enddata=None, write_header=True, close=True,
                            flfact_size=8)

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
    return model, density_units, velocity_units


def _float(data: dict[str, Any], name: str):
    svalue = data[name]
    try:
        value = float(svalue)
    except:
        raise SyntaxError(f'name={name} value={svalue!r} is not a float')
    return value


def _int(data: dict[str, Any], name: str):
    svalue = data[name]
    try:
        value = int(svalue)
    except:
        raise SyntaxError(f'name={name} value={svalue!r} is not an integer')
    return value
