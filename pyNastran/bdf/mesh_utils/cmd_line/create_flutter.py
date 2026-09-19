from __future__ import annotations
import sys
from typing import Optional, TYPE_CHECKING
import numpy as np

from cpylog import SimpleLogger
from .utils import filter_no_args
from pyNastran.utils.atmosphere import atm_density, atm_speed_of_sound, atm_temperature
from pyNastran.utils.convert import convert_altitude, convert_velocity

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

UNITS_MAP = {
    # units must be consistent
    #
    # (alt, velocity, density, eas, pressure)
    "english_in": ("ft", "in/s", "slinch/in^3", "knots", "psi"),
    "english_ft": ("ft", "ft/s", "slug/ft^3", "knots", "psf"),
    "si": ("m", "m/s", "kg/m^3", "m/s", "Pa"),
    "si_mm": ("m", "mm/s", "Mg/mm^3", "mm/s", "MPa"),
}

ALT_UNITS = ["m", "ft", "kft"]
VELOCITY_UNITS = ["m/s", "cm/s", "mm/s", "in/s", "ft/s", "knots"]


def cmd_line_create_flutter(argv=None, quiet: bool = False) -> None:
    """command line interface to bdf flutter"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    import pyNastran

    filter_no_args("bdf flutter: use 'bdf flutter -h' for help", argv, quiet=quiet)

    if "gui" in argv or '--gui' in argv:
        from pyNastran.bdf.mesh_utils.gui_tools.gui_flutter import cmd_line_gui
        data = cmd_line_gui()
        return

    ver = str(pyNastran.__version__)

    # shared optional arguments — defined on a parent parser so that
    # each subparser inherits them and they are recognised even when
    # they appear *after* the subparser's positional arguments
    option_parser = argparse.ArgumentParser(add_help=False)
    option_parser.add_argument('-o', '--output', default='flutter_cards.inc',
                               metavar='OUT_BDF_FILENAME',
                               help='path to output BDF/DAT/NAS file (default=flutter_cards.inc)')

    size_group = option_parser.add_mutually_exclusive_group()
    size_group.add_argument('--size', type=int, default=None,
                            help='size of the BDF (8/16; default=16)')
    size_group.add_argument('--clean', action='store_true',
                            help='writes a BDF with at least 1 whitespace in an FLFACT field (for readability)')

    option_parser.add_argument('--sid', type=int, default=None,
                               help='updates the flutter ID')
    option_parser.add_argument('--rhoref', action='store_true',
                               help='rhoref flag')
    option_parser.add_argument('--minus_eas', default=None, metavar='MINUS_EAS',
                               help="request flutter mode shapes at the closest point (e.g., '400,500')")
    option_parser.add_argument('--zaero', action='store_true',
                               help='zaero flag')
    option_parser.add_argument('--eas_limit', nargs=2, default=None,
                               metavar=('EAS', 'EAS_UNITS'),
                               help='EAS limit value and units')

    parser = argparse.ArgumentParser(
        prog='bdf flutter',
        description='Create flutter analysis cards',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  bdf flutter english_in eas  1   800. knots 101 mach 0.8 na\n"
            "  bdf flutter english_in eas  1   800. knots 101 mach 0.8 na --minus_eas '100,200'\n"
            "  bdf flutter english_in tas  .1  800. ft/s  101 alt 2500 m\n"
            "  bdf flutter english_in mach .05 0.5        101 alt 2500\n"
            "  bdf flutter english_in mach .05 0.5        101 alt 2500 m "
            "--eas_limit 300 knots --output flutter_cards_temp.inc --size 16\n"
        ),
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('units', metavar='UNITS',
                        choices=list(UNITS_MAP.keys()),
                        help='unit system: english_in, english_ft, si, si_mm')

    subparsers = parser.add_subparsers(dest='sweep_method',
                                       help='sweep method: eas, tas, alt, mach')
    subparsers.required = True

    def _add_sweep_positionals(sub, val1_name, val2_name, has_sweep_unit=True):
        """add the common positional arguments to a sweep subparser"""
        sub.add_argument('value1', type=float, metavar=val1_name)
        sub.add_argument('value2', type=float, metavar=val2_name)
        if has_sweep_unit:
            sub.add_argument('sweep_unit', metavar='SWEEP_UNIT',
                             help='unit for sweep variable')
        sub.add_argument('n', type=int, metavar='N',
                         help='number of points in the sweep')
        sub.add_argument('const_type', metavar='CONST_TYPE',
                         help='parameter held constant (alt, mach, eas, tas)')
        sub.add_argument('const_val', type=float, metavar='CONST_VAL',
                         help='value for constant parameter')
        sub.add_argument('const_unit', metavar='CONST_UNIT',
                         help='unit for constant parameter')

    eas_sub = subparsers.add_parser('eas', parents=[option_parser],
                                    help='sweep equivalent airspeed')
    _add_sweep_positionals(eas_sub, 'EAS1', 'EAS2')

    tas_sub = subparsers.add_parser('tas', parents=[option_parser],
                                    help='sweep true airspeed')
    _add_sweep_positionals(tas_sub, 'TAS1', 'TAS2')

    alt_sub = subparsers.add_parser('alt', parents=[option_parser],
                                    help='sweep altitude')
    _add_sweep_positionals(alt_sub, 'ALT1', 'ALT2')

    mach_sub = subparsers.add_parser('mach', parents=[option_parser],
                                     help='sweep Mach number')
    _add_sweep_positionals(mach_sub, 'MACH1', 'MACH2', has_sweep_unit=False)

    argv = [str(arg) for arg in argv]
    cmd = "bdf " + " ".join(argv[1:])
    args = parser.parse_args(argv[2:])

    if not quiet:  # pragma: no cover
        print(vars(args))

    size = args.size if args.size is not None else 16
    sid = args.sid if args.sid is not None else 1

    minus_eas = []
    if args.minus_eas is not None:
        svalue = args.minus_eas
        if "," in svalue:
            svalues = svalue.strip(",").split(",")
            minus_eas = [float(sv) for sv in svalues]
        else:
            minus_eas = [float(svalue)]

    units_out = args.units
    if units_out.lower() not in UNITS_MAP:  # pragma: no cover
        raise NotImplementedError(units_out)

    is_zaero = args.zaero
    rhoref_flag = args.rhoref
    npoints = args.n
    clean = args.clean

    const_type = args.const_type.lower()
    assert const_type in {"alt", "mach", "eas", "tas"}, f"const_type={const_type!r}"
    const_value = args.const_val
    const_unit = args.const_unit.lower()

    eas_units = ""
    eas_limit = 1_000_000.0
    if args.eas_limit is not None:
        eas_limit = float(args.eas_limit[0])
        eas_units = args.eas_limit[1]
        assert eas_units not in {None, ""}, eas_units
        eas_units = eas_units.lower()

    sweep_method = args.sweep_method
    sweep_unit = getattr(args, 'sweep_unit', '')
    if sweep_unit is None:
        sweep_unit = ''
    value1 = args.value1
    value2 = args.value2

    bdf_filename_out = args.output

    level = "debug" if not quiet else "warning"
    log = SimpleLogger(level=level, encoding="utf-8")
    create_flutter(
        log,
        sweep_method,
        value1,
        value2,
        sweep_unit,
        npoints,
        const_type,
        const_value,
        const_unit,
        eas_limit=eas_limit,
        eas_units=eas_units,
        units_out=units_out,
        rhoref_flag=rhoref_flag,
        sid=sid,
        size=size,
        clean=clean,
        bdf_filename_out=bdf_filename_out,
        minus_eas=minus_eas,
        is_zaero=is_zaero,
        comment=cmd,
    )
    if not quiet:
        print(cmd)


def create_flutter(
    log: SimpleLogger,
    sweep_method: str,
    value1: float,
    value2: float,
    sweep_unit: str,
    npoints: int,
    const_type: str,
    const_value: float,
    const_unit: str,
    sid: int = 1,
    eas_limit: float = 1_000_000.0,
    eas_units: str = "m/s",
    units_out: str = "si",
    rhoref_flag: bool = False,
    size: int = 8,
    clean: bool = False,
    minus_eas: Optional[list[float]] = None,
    is_zaero: bool = False,
    bdf_filename_out: str = "flutter_cards.inc",
    comment: str = "",
) -> tuple[BDF, str, str]:
    if minus_eas is None:
        minus_eas = []

    unitsi = UNITS_MAP[units_out.lower()]
    alt_units, velocity_units, density_units, eas_units_default, pressure_units = unitsi
    if const_type == "alt":
        alt = const_value
        alt_unit = const_unit
        assert alt_unit in ALT_UNITS, f"alt_unit={alt_unit!r}; allowed={ALT_UNITS}"
        alt = convert_altitude(alt, alt_unit, alt_units)
    elif const_type == "mach":
        mach = const_value
        assert const_unit in {"none", "na"}, f"const_unit={const_unit!r}; allowed=[none]"
    elif const_type == "eas":
        eas = const_value
        assert const_unit in VELOCITY_UNITS, f"const_unit={const_unit!r}; allowed={VELOCITY_UNITS}"
    elif const_type == "tas":
        tas = const_value
        assert const_unit in VELOCITY_UNITS, f"const_unit={const_unit!r}; allowed={VELOCITY_UNITS}"
    else:  # pragma: no cover
        raise NotImplementedError(f"const_type={const_type} is not supported")

    sweep_unit = sweep_unit.lower()
    values = np.linspace(value1, value2, num=npoints)
    dvalue1 = values[1] - values[0]
    if sweep_method == "alt":
        assert sweep_unit in ALT_UNITS, f"sweep_unit={sweep_unit!r}; allowed={ALT_UNITS}"
        alts = convert_altitude(values, sweep_unit, alt_units)
        unit2 = alt_units

    elif sweep_method == "mach":
        machs = values
        unit2 = "na"
    elif sweep_method == "eas":
        eass = values
        eas_units = sweep_unit
        unit2 = eas_units
        # assert sweep_unit in VELOCITY_UNITS, f'sweep_unit={sweep_unit!r}; allowed={VELOCITY_UNITS}'
        # eass = convert_velocity(eass, sweep_unit, velocity_units)

    elif sweep_method == "tas":
        assert sweep_unit in VELOCITY_UNITS, f"sweep_unit={sweep_unit!r}; allowed={VELOCITY_UNITS}"
        tass = convert_velocity(values, sweep_unit, velocity_units)
        unit2 = velocity_units
    elif sweep_method == "alt":
        assert sweep_unit in ALT_UNITS, f"sweep_unit={sweep_unit!r}; allowed={ALT_UNITS}"
        alts = convert_altitude(values, sweep_unit, alt_units)
        unit2 = alt_units
    else:  # pragma: no cover
        raise NotImplementedError(sweep_method)
    dvalue2 = values[1] - values[0]

    # ------------------------------------------------------------------
    from pyNastran.bdf.bdf import BDF

    model = BDF(log=log)
    model.set_error_storage(
        nparse_errors=100, stop_on_parsing_error=True, nxref_errors=100, stop_on_xref_error=False
    )
    flutter_method = "PKNL"

    flfact_density = sid + 1
    flfact_mach = sid + 2
    flfact_velocity = sid + 3
    # flfact_eas = sid + 4

    flutter = model.add_flutter(
        sid,
        flutter_method,
        flfact_density,
        flfact_mach,
        flfact_velocity,
        imethod="L",
        nvalue=None,
        omax=None,
        epsilon=1.0e-3,
        comment="",
        validate=True,
    )

    if eas_units in {None, ""}:
        log.debug(f"setting eas_units to default; eas_units={eas_units_default!r}")
        eas_units = eas_units_default

    log.info(f"sweep_method={sweep_method!r}")
    log.debug(f"  d{sweep_unit} = {dvalue1!r} {sweep_unit} = {dvalue2!r} {unit2}")
    log.debug(f"  alt_units={alt_units!r}")
    log.debug(f"  velocity_units={velocity_units!r}")
    log.debug(f"  density_units={density_units!r}")
    log.debug(f"  eas_units={eas_units!r}")
    # del alt_unit

    # option 1: overwrite the eas/tas/alt unit...would work for alt
    #           we'll overwrite the
    # option 2: pass another flag in...I don't wanna
    pairs = [
        # sweep, const
        ("eas", "alt"),
        ("eas", "mach"),
        ("mach", "alt"),
        ("alt", "mach"),
        ("tas", "alt"),
        # ('tas', 'mach'),  # undefined b/c sos doesn't map 1:1 to alt
        # new
        # ('tas', 'eas'),  # kind of dumb, but should work...
        ("alt", "tas"),
        # ('alt', 'eas'), # kind of dumb
    ]

    assert alt_units != "", alt_units
    assert velocity_units != "", velocity_units
    assert density_units != "", density_units
    assert eas_units != "", eas_units
    assert pressure_units != "", pressure_units

    if sweep_method == "eas" and const_type == "alt":
        alts, machs, eass = flutter.make_flfacts_eas_sweep_constant_alt(
            model, alt, eass,
            minus_eas=minus_eas,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    elif sweep_method == "eas" and const_type == "mach":
        gamma = 1.4
        alts, machs, eass = (
            flutter.make_flfacts_eas_sweep_constant_mach(  # TODO: need to test this; seems wrong
                model, mach, eass,
                gamma=gamma,
                minus_eas=minus_eas,
                alt_units=alt_units,
                velocity_units=velocity_units,
                density_units=density_units,
                pressure_units=pressure_units,
                eas_units=eas_units,
            )
        )

    elif sweep_method == "mach" and const_type == "alt":
        alts, machs, eass = flutter.make_flfacts_mach_sweep_constant_alt(
            model, alt, machs,
            minus_eas=minus_eas,
            eas_limit=eas_limit,
            alt_units=alt_unit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    elif sweep_method == "alt" and const_type == "mach":
        # alt_units = sweep_unit
        alts, machs, eass = flutter.make_flfacts_alt_sweep_constant_mach(
            model, mach, alts,
            minus_eas=minus_eas,
            eas_limit=eas_limit,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    elif sweep_method == "alt" and const_type == "tas":
        alts, machs, eass = flutter.make_flfacts_alt_sweep_constant_tas(
            model, tas, alts,
            minus_eas=minus_eas,
            alt_units=alt_units,
            eas_limit=eas_limit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    elif sweep_method == "tas" and const_type == "alt":
        alts, machs, eass = flutter.make_flfacts_tas_sweep_constant_alt(
            model, alt, tass,
            minus_eas=minus_eas,
            alt_units=alt_units,
            eas_limit=eas_limit,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    # elif sweep_method == 'tas' and const_type == 'mach':
    # undefined b/c sos doesn't map 1:1 to alt
    # flutter.make_flfacts_tas_sweep_constant_mach(
    # model, mach, tass,
    # alt_units=alt_units,
    # eas_limit=eas_limit,
    # velocity_units=velocity_units,
    # density_units=density_units,
    # eas_units=eas_units)
    # elif sweep_method == 'tas' and const_type == 'eas':
    # flutter.make_flfacts_tas_sweep_constant_eas(
    # model, eas, tass,
    # alt_units=alt_units,
    # eas_limit=eas_limit,
    # velocity_units=velocity_units,
    # density_units=density_units,
    # eas_units=eas_units)
    elif sweep_method == "alt" and const_type == "eas":
        alts, machs, eass = flutter.make_flfacts_alt_sweep_constant_eas(
            model, eas, alts,
            minus_eas=minus_eas,
            alt_units=alt_units,
            velocity_units=velocity_units,
            density_units=density_units,
            eas_units=eas_units,
        )
    else:  # pragma: no cover
        raise NotImplementedError((sweep_method, const_type))

    if is_zaero:
        model = BDF()
        model.set_as_zaero()
        atm_id = sid
        mass_unit = "slinch"
        length_unit = "in"
        temperature_units = "R"
        velocity_units = f"{length_unit}/s"
        density_units = f"{mass_unit}/{length_unit}^3"
        nalt = len(alts)
        sos = np.array([
            atm_speed_of_sound(alt, alt_units=alt_units, velocity_units=velocity_units)
            for alt in alts])
        density = np.array([
            atm_density(alt, 1716.0, alt_units=alt_units, density_units=density_units)
            for alt in alts])
        temperature = np.array([
            atm_temperature(alt, alt_units=alt_units, temperature_units=temperature_units)
            for alt in alts])

        assert length_unit in ['in', 'm'], length_unit
        alt_scale = 1.
        if length_unit == 'in':
            alt_scale = 12.
        # density0 = atm_density(0.0, 1716.0, alt_units=alt_units, density_units=density_units)
        # round density to 4 sig figs
        density_magnitude = np.floor(np.log10(np.abs(density)))
        decimals = 3 - density_magnitude
        decimals = np.nan_to_num(decimals, nan=0.0, posinf=0.0, neginf=0.0).astype(int)
        density2 = np.array([np.round(x, d) for x, d in zip(density, decimals)])
        alts2 = alts.round(0) * alt_scale
        # print(f'alts = {alts} ({alt_units})')
        sos2 = sos.round(0)
        atmosphere_table = np.column_stack([
                alts2, sos2, density2, temperature.round(2),])
        # low to high
        atmosphere_table_reversed = atmosphere_table[::-1, :]

        atmosphere_list = atmosphere_table_reversed.ravel().tolist()
        msg = (
            f"\n alt: ({alts2.min()}, {alts2.max()}) ({length_unit})"
            f"\n sos: ({sos2.min()}, {sos2.max()}) ({velocity_units})"
            f"\n rho: ({density2.min()}, {density2.max()}) ({density_units})"
        )
        sweep_id = sid + 1
        mkaeroz_id = sid + 2
        fluttf_id = 0
        print_flag = 0
        if sweep_method == 'alt':
            comment = (
                f"\n alt: ({alts.min()}, {alts.max()}) ({alt_units})\n"
                f"\n alt: ({alts2.min()}, {alts2.max()}) ({length_unit})\n")
            alts2.sort()
            model.zaero.add_fixmatm(
                sweep_id, mkaeroz_id,  mass_unit, length_unit, alts2,
                atm_id=0, fluttf_id=fluttf_id, print_flag=print_flag, vref=1.0, comment=comment)
        else:
            model.zaero.add_atmos(
                atm_id,
                mass_unit,
                length_unit,
                temperature_units,
                atmosphere_list,
                comment=comment + msg,
            )
            alts2.sort()
            model.zaero.add_fixmatm(
                sweep_id, mkaeroz_id,  mass_unit, length_unit, alts2,
                atm_id=atm_id,
                fluttf_id=fluttf_id, print_flag=print_flag, vref=1.0, comment="")

        freqs = [0.1, 0.2, 0.5, 1.]
        model.zaero.add_mkaeroz(mkaeroz_id, mach, freqs, comment="")
        # velocity = mach * sos
        # model.zaero.add_fixmach(
        #     sid, mkaeroz_id, mass_unit, length_unit,
        #     fluttf_id, print_flag, velocity, density2,
        #     vref=1.0, comment="")

    if rhoref_flag:
        rho0 = atm_density(alt=0.0, density_units=density_units)
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
            model.write_bdf(
                bdf_filename_out,
                encoding=None,
                size=8,
                nodes_size=None,
                elements_size=None,
                loads_size=None,
                is_double=False,
                interspersed=False,
                enddata=None,
                write_header=True,
                close=True,
                flfact_size=8,
            )

            model2 = BDF(log=log)
            model2.read_bdf(
                bdf_filename_out,
                validate=True,
                xref=False,
                punch=True,
                read_includes=True,
                save_file_structure=False,
                encoding=None,
            )
            model2.write_bdf(
                bdf_filename_out,
                encoding=None,
                size=16,
                nodes_size=None,
                elements_size=None,
                loads_size=None,
                is_double=False,
                interspersed=False,
                enddata=None,
                write_header=True,
                close=True,
            )
        else:
            model.write_bdf(
                bdf_filename_out,
                encoding=None,
                size=size,
                nodes_size=None,
                elements_size=None,
                loads_size=None,
                is_double=False,
                interspersed=False,
                enddata=None,
                write_header=True,
                close=True,
            )
    return model, density_units, velocity_units


