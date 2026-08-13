import os
from io import StringIO
from typing import Optional, cast
import numpy as np
from cpylog import SimpleLogger, get_logger
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.bdf.bdf import BDF
from pyNastran.f06.flutter_response import FlutterResponse


def read_zaero_out(zaero_out_filename: PathLike,
                   log: Optional[SimpleLogger] = None,
                   debug: Optional[str | bool] = True,) -> tuple[dict[int, FlutterResponse], dict]:
    """
    Returns
    -------
    data_dict : dict
        opgwg : dict[str, Any]
            ref_point : int
            MO: (6?,6?) float np.ndarray
            cg : (3,1) float np.ndarray
            mass : float
            mass_error: float
            cg_error: float...why???
            I(S): (6?,6?) float np.ndarray
        matrices
            freq : (nmodes,) float np.ndarray
            MHH : (nmodes,nmodes) float np.ndarray
            BHH : (nmodes,nmodes) float np.ndarray
            KHH : (nmodes,nmodes) float np.ndarray

    """
    log = get_logger(log, debug)
    assert os.path.exists(zaero_out_filename), print_bad_path(zaero_out_filename)
    with open(zaero_out_filename, "r") as zaero_out_file:
        lines = zaero_out_file.readlines()

    (case_dict, ref_dict, exec_control_lines, case_control_lines, bulk_data_lines, modal_data, trim_data) = (
        zaero_lines_to_out(log, lines)
    )

    bdf_model = _build_bdf_from_echo(
        exec_control_lines, case_control_lines, bulk_data_lines, log)

    if trim_data is not None:
        # TRIM discipline — no flutter V-G-F results
        responses = {}
        data_dict = {
            "matrices": modal_data,
            "model": bdf_model,
            "trim": trim_data,
        }
        return responses, data_dict

    for key, data in case_dict.items():
        log.debug(f"{key}:\n{str(data)}")

    # data = out_to_data()
    modes, result, in_units = out_dict_to_results(case_dict, ref_dict, log)

    # FlutterResponse.from_nx()
    assert isinstance(in_units, dict), in_units
    log.debug(f'creating flutter response in_units={in_units}')
    resp = FlutterResponse.from_zaero(modes, result, in_units, zaero_out_filename)
    responses = {
        1: resp,
    }

    data_dict = {
        "matrices": modal_data,
        "model": bdf_model,
    }
    return responses, data_dict


def _build_bdf_from_echo(exec_control_lines: list[str],
                         case_control_lines: list[str],
                         bulk_data_lines: list[str],
                         log: SimpleLogger,) -> BDF:
    """Build a BDF from ZAERO .out executive/case control and bulk data echo."""
    bdf_lines = []
    for line in exec_control_lines:
        if line.strip().upper().startswith("ASSIGN"):
            continue
        bdf_lines.append(line)
    bdf_lines.append("CEND")
    for line in case_control_lines:
        bdf_lines.append(line)
    bdf_lines.append("BEGIN BULK")
    for line in bulk_data_lines:
        bdf_lines.append(line)
    bdf_lines.append("ENDDATA")
    # print('\n'.join(bdf_lines))

    bdf_string = "\n".join(bdf_lines)
    model = BDF(log=log)
    model.set_as_zaero()
    model.read_bdf(StringIO(bdf_string), xref=False, punch=False)
    return model


def out_dict_to_results(
    out_dict: dict, ref_dict: dict[str, tuple[float, str]], log: SimpleLogger
) -> tuple[list[int], np.ndarray, dict[str, str]]:
    header_name = out_dict["header"]
    header_data_strs = out_dict[header_name]  # vvref,density,mach

    # v/vref, density, mach
    header_data = np.array(header_data_strs, dtype="float64")
    nrows = len(header_data)

    if header_name == "v/vref,density,mach":
        vref, velocity_units = ref_dict["VREF"]
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        mach = header_data[:, 2]
        density_units = out_dict["density_units"]
    elif header_name == "v/vref,density,q":
        vref, velocity_units = ref_dict["VREF"]
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        q = header_data[:, 2]
        mach = q
        density_units = out_dict["density_units"]
    elif header_name == "v/vref,density,alt":
        # ref_dict = {'MACH': (0.8, ''), 'ATMOS TABLE': ('STANDARD', ''),
        #             'REFERENCE LENGTH (L)': (1.0, 'in'),
        #             'VREF': (2.0, 'in/s')}
        # print('ref_dict = ', ref_dict)
        vref, velocity_units = ref_dict["VREF"]
        atmos_table = ref_dict["ATMOS TABLE"][0]
        # print(f'atmos_table = {atmos_table!r}')
        if isinstance(atmos_table, str):
            assert atmos_table == "STANDARD", atmos_table
        else:
            assert isinstance(atmos_table, int), atmos_table
        mach = ref_dict["MACH"][0]
        # alt, altitude_units = ref_dict['ALT']
        # print(f'vref={vref}; velocity_units={velocity_units}')
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        q = 0.5 * density * velocity**2
        density_units = out_dict["density_units"]
    elif header_name == "v/vref,v,q":
        # print('ref_dict = ', ref_dict)
        # velocity0, velocity_units = ref_dict['VELOCITY']
        # atmos_table = ref_dict['ATMOS TABLE'][0]
        # assert atmos_table == 'STANDARD', atmos_table
        # mach = ref_dict['MACH'][0]
        # alt, altitude_units = ref_dict['ALT']
        # print(f'vref={vref}; velocity_units={velocity_units}')
        v_vref = header_data[:, 0]
        # print((velocity - velocity0).max())
        q = header_data[:, 1]
        if "DENSITY" in ref_dict and "VREF" in ref_dict and "MACH" in ref_dict:
            vref, velocity_units = ref_dict["VREF"]
            velocity = v_vref * vref
            density, density_units = ref_dict["DENSITY"]
            mach = ref_dict["MACH"][0]
        else:
            raise RuntimeError(f"V,q; {ref_dict}")
        log.debug(f"density = {density}")
    else:  # pragma: no cover
        raise NotImplementedError(header_name)

    if velocity_units in "in/s":
        altitude_units = "ft"
        dynamic_pressure_units = "psi"
    elif velocity_units in "ft/s":
        altitude_units = "ft"
        dynamic_pressure_units = "psf"
    else:
        assert velocity_units == "m/s", velocity_units
        altitude_units = "m"
        dynamic_pressure_units = "Pa"

    # print(f'out_dict = {list(out_dict.keys())}')
    # print(f'out_dict[header] = {out_dict['header']}')
    in_units = {
        "velocity": velocity_units,
        "density": density_units,
        "eas": velocity_units,
        "altitude": altitude_units,
        "dynamic_pressure": dynamic_pressure_units,
    }

    modes = [mode for mode in out_dict if isinstance(mode, int)]
    nmodes = len(modes)
    log.info(f"modes = {modes}")
    results_list = []
    result = np.zeros((nmodes, nrows, 9), dtype="float64")
    log.debug(f"nmodes={nmodes}; nrows={nrows}")
    for imode, mode in enumerate(modes):
        # g, f(Hz), k
        data = out_dict[mode]
        # print(f'data = {data}')
        g_freq_k = np.array(data, dtype="float64")
        damping = g_freq_k[:, 0]
        freq = g_freq_k[:, 1]
        kfreq = g_freq_k[:, 2]

        omega = 2 * np.pi * freq
        eigr = omega * damping / 2
        eigi = omega

        # fdata   [                rho,     mach, velocity, damping, freq]
        # results [kfreq, 1/kfreq, density, mach, velocity, damping, freq, eigr, eigi]
        # log.debug(len(kfreq))
        result[imode, :, 0] = kfreq
        result[imode, :, 1] = 1 / kfreq
        result[imode, :, 2] = density
        result[imode, :, 3] = mach
        result[imode, :, 4] = velocity
        result[imode, :, 5] = damping
        result[imode, :, 6] = freq
        result[imode, :, 7] = eigr
        result[imode, :, 8] = eigi
    # for mode_name in out_dict.keys():
    #     if isinstance()
    # for mode in out_dict['modes']:
    return modes, result, in_units


def _read_ruler_section(lines: list[str], iline: int, nlines: int,
                        out_lines: list[str], stop_on: str) -> int:
    """Read lines between a |...1...| ruler and a stop keyword (e.g. CEND, BEGIN BULK)."""
    found_ruler = False
    while iline < nlines:
        stripped = lines[iline].strip()
        if stripped.startswith("|...1...|"):
            found_ruler = True
            iline += 1
            continue
        if not found_ruler:
            iline += 1
            continue
        if stripped == "":
            iline += 1
            continue
        if stripped.startswith(stop_on):
            iline += 1
            break
        out_lines.append(stripped)
        iline += 1
    return iline


def _parse_modal_table(lines: list[str], iline: int, nlines: int,
                       log: SimpleLogger) -> dict[str, np.ndarray]:
    """Parse the ZAERO modal eigenvalue table into arrays.

    Returns
    -------
    modal_data : dict[str, np.ndarray]
        mode : (nmodes,) int
        extraction_order : (nmodes,) int
        eigenvalue : (nmodes,) float — (rad/s)^2
        omega : (nmodes,) float — rad/s
        freq : (nmodes,) float — Hz
        generalized_mass : (nmodes,) float
        generalized_stiffness : (nmodes,) float
    """
    # find the header
    while iline < nlines:
        if "MODE      EXTRACTION      EIGENVALUE" in lines[iline]:
            iline += 2  # skip header + units line
            break
        if "****************************************" in lines[iline]:
            return {}
        iline += 1

    rows: list[list[float]] = []
    while iline < nlines:
        stripped = lines[iline].strip()
        if stripped == "":
            iline += 1
            break
        if stripped.startswith("MASS UNIT"):
            iline += 1
            continue
        if "========" in stripped:
            break
        values = stripped.split()
        if len(values) >= 7:
            rows.append([float(v) for v in values[:7]])
        iline += 1

    if not rows:
        return {}, iline

    data = np.array(rows, dtype="float64")
    log.debug(f"parsed {len(rows)} modal eigenvalues")
    Mhh = data[:, 5]
    Khh = data[:, 6]
    isort = np.argsort(Khh)

    nmodes = len(rows)
    modal_data = {
        "mode": data[:, 0].astype(int),
        "extraction_order": data[:, 1].astype(int),
        "eigenvalue": data[:, 2],
        "omega": data[:, 3],
        "freq": data[isort, 4],
        "generalized_mass": Mhh,
        "generalized_stiffness": Khh,
        "MHH": np.diag(Mhh[isort]),
        "BHH": np.zeros((nmodes, nmodes), dtype="float64"),
        "KHH": np.diag(Khh[isort]),
    }
    return modal_data, iline


def _parse_trim_section(
    lines: list[str], iline: int, nlines: int, log: SimpleLogger,
) -> dict:
    """Parse the ZAERO trim results section.

    Returns a dict with keys such as ``mach``, ``dynamic_pressure``,
    ``dynamic_pressure_units``, ``stability_derivatives``,
    ``trim_results``, ``trim_variables``, and ``aero_forces``.
    """
    trim: dict = {}

    # --- scan for MACH NUMBER and DYNAMIC PRESSURE header ---
    while iline < nlines:
        line = lines[iline]
        if "MACH NUMBER =" in line:
            # ' MACH NUMBER =  0.9540. STEADY ...'
            parts = line.split("MACH NUMBER =")[1]
            mach_str = parts.split(".")[0] + "." + parts.split(".")[1]
            # mach_str is e.g. '  0.9540'
            trim["mach"] = float(mach_str.strip().rstrip("."))
            iline += 1
            break
        iline += 1

    while iline < nlines:
        line = lines[iline]
        if "DYNAMIC PRESSURE=" in line:
            # ' DYNAMIC PRESSURE= 0.12000E+04  LBF/IN**2'
            after = line.split("DYNAMIC PRESSURE=")[1].strip()
            tokens = after.split()
            trim["dynamic_pressure"] = float(tokens[0])
            trim["dynamic_pressure_units"] = tokens[1] if len(tokens) > 1 else ""
            iline += 1
            break
        iline += 1

    # --- stability derivatives table ---
    # The table looks like:
    #  -----...-----
    #  | IDVAR  | LABEL  | DRAG COEFFICIENT| ...
    #  |        |        |  RIGID |FLEXIBLE| ...
    #  -----...-----
    #  |     100|ALPHA   | 0.00000| 0.00000| ...
    #  |  UNITS= 1/DEG   | E/R= ...
    #  -----...-----
    stability_derivatives: dict[str, dict[str, float]] = {}
    while iline < nlines:
        line = lines[iline]
        if line.strip().startswith("| IDVAR"):
            # skip header row, sub-header row, then dashed separator
            iline += 3
            break
        iline += 1

    # now read data rows until the closing dashed line
    while iline < nlines:
        line = lines[iline].strip()
        if line.startswith("-----"):
            iline += 1
            break
        if line.startswith("|") and "UNITS" not in line:
            # '|     100|ALPHA   | 0.00000| 0.00000| ...'
            cells = [c.strip() for c in line.split("|")]
            cells = [c for c in cells if c]
            if len(cells) >= 14:
                label = cells[1]
                stability_derivatives[label] = {
                    "drag_rigid": float(cells[2]),
                    "drag_flexible": float(cells[3]),
                    "side_force_rigid": float(cells[4]),
                    "side_force_flexible": float(cells[5]),
                    "lift_rigid": float(cells[6]),
                    "lift_flexible": float(cells[7]),
                    "roll_moment_rigid": float(cells[8]),
                    "roll_moment_flexible": float(cells[9]),
                    "pitch_moment_rigid": float(cells[10]),
                    "pitch_moment_flexible": float(cells[11]),
                    "yaw_moment_rigid": float(cells[12]),
                    "yaw_moment_flexible": float(cells[13]),
                }
        iline += 1

    trim["stability_derivatives"] = stability_derivatives

    # --- trim results ---
    trim_results: dict[str, dict] = {}
    while iline < nlines:
        line = lines[iline]
        if "T R I M   R E S U L T S" in line:
            iline += 1
            break
        if "***  Z A E R O   T E R M I N A T E D ***" in line:
            break
        iline += 1

    # look for COMPUTED: lines
    while iline < nlines:
        line = lines[iline]
        if "COMPUTED:" in line:
            # '      COMPUTED:   NZ         SYMMETRIC   2.4624E-01 G   ...'
            tokens = line.split()
            idx = tokens.index("COMPUTED:")
            label = tokens[idx + 1]
            symmetry = tokens[idx + 2]
            # flexible value is next, then units, then rigid, then units
            flexible_val = float(tokens[idx + 3])
            flexible_units = tokens[idx + 4]
            rigid_val = float(tokens[idx + 5])
            rigid_units = tokens[idx + 6]
            trim_results[label] = {
                "symmetry": symmetry,
                "flexible": flexible_val,
                "flexible_units": flexible_units,
                "rigid": rigid_val,
                "rigid_units": rigid_units,
            }
        if "NUMBER OF TRIM VARIABLES" in line:
            iline += 1
            break
        if "***  Z A E R O   T E R M I N A T E D ***" in line:
            break
        iline += 1

    trim["trim_results"] = trim_results

    # --- trim variables ---
    trim_variables: dict[str, dict] = {}
    while iline < nlines:
        line = lines[iline]
        if "USER INPUT:" in line:
            tokens = line.split()
            idx = tokens.index("INPUT:")
            idvar = int(tokens[idx + 1])
            label = tokens[idx + 2]
            symmetry = tokens[idx + 3]
            flexible_val = float(tokens[idx + 4])
            rigid_val = float(tokens[idx + 5])
            units = tokens[idx + 6]
            trim_variables[label] = {
                "idvar": idvar,
                "symmetry": symmetry,
                "flexible": flexible_val,
                "rigid": rigid_val,
                "units": units,
            }
        if "S U M M A R Y   O F   T O T A L" in line:
            iline += 1
            break
        if "***  Z A E R O   T E R M I N A T E D ***" in line:
            break
        iline += 1

    trim["trim_variables"] = trim_variables

    # --- aero forces summary ---
    aero_forces: dict[str, dict[str, float]] = {}
    while iline < nlines:
        line = lines[iline]
        stripped = line.strip()
        if "***  Z A E R O   T E R M I N A T E D ***" in line:
            break
        if "M O D A L   C O O R D I N A T E S" in line:
            break
        if ":" in stripped and stripped and not stripped.startswith("COEFFICIENTS"):
            # 'INDUCED DRAG(CDL):             0.03363             0.34544     FX/REFS/Q'
            parts = stripped.split(":")
            coeff_name = parts[0].strip()
            vals = parts[1].split()
            if len(vals) >= 3:
                try:
                    flexible_val = float(vals[0])
                    rigid_val = float(vals[1])
                    units = vals[2]
                    aero_forces[coeff_name] = {
                        "flexible": flexible_val,
                        "rigid": rigid_val,
                        "units": units,
                    }
                except (ValueError, IndexError):
                    pass
        iline += 1

    trim["aero_forces"] = aero_forces

    log.info(f"parsed trim data: mach={trim.get('mach')}, "
             f"q={trim.get('dynamic_pressure')}, "
             f"{len(stability_derivatives)} derivatives, "
             f"{len(trim_variables)} trim variables, "
             f"{len(aero_forces)} aero forces")
    return trim


def zaero_lines_to_out(log: SimpleLogger,
                       lines: list[str]) -> tuple[dict, dict, list[str], list[str], list[str], dict, Optional[dict]]:
    # print(log)
    out = {}
    nlines = len(lines)
    iline = 0

    exec_control_lines: list[str] = []
    case_control_lines: list[str] = []
    bulk_data_lines: list[str] = []

    # --- parse EXECUTIVE CONTROL SUMMARY ---
    while iline < nlines:
        if "E X E C U T I V E  C O N T R O L  S U M M A R Y" in lines[iline]:
            iline += 1
            log.debug('reading executive control summary')
            break
        iline += 1
    iline = _read_ruler_section(lines, iline, nlines, exec_control_lines, stop_on="CEND")

    # --- parse CASE CONTROL SUMMARY ---
    while iline < nlines:
        if "C A S E  C O N T R O L  S U M M A R Y" in lines[iline]:
            iline += 1
            log.debug('reading case control summary')
            break
        iline += 1
    iline = _read_ruler_section(lines, iline, nlines, case_control_lines, stop_on="BEGIN BULK")
    line = lines[iline].rstrip()
    # log.info(f'AA {iline}: {line}')

    # --- parse SORTED BULK DATA ECHO ---
    while iline < nlines:
        line = lines[iline]
        if "S O R T E D   B U L K   D A T A   E C H O" in line:
            iline += 1
            log.debug('reading sorted bulk data echo')
            break
        if "SYMMETRIC (OR ASYMMETRIC) FINITE ELEMENT MODAL RESULTS ARE SUCCESSFULLY READ IN FROM FILE" in line:
            break
        iline += 1

    # skip CARD / COUNT / ruler header lines
    # "MODE      EXTRACTION      EIGENVALUE"
    while iline < nlines:
        stripped = lines[iline].strip()
        # log.info(f'A {iline}: {stripped}')
        if stripped in ("", "CARD") or stripped.startswith(("COUNT", "|...1...|")):
            iline += 1
            continue
        if "SYMMETRIC (OR ASYMMETRIC) FINITE ELEMENT MODAL RESULTS ARE SUCCESSFULLY READ IN FROM FILE" in line:
            break
        break

    # read bulk cards until ENDDATA or blank section
    # the card content after "N -" is 80-char fixed-width; preserve column alignment
    while iline < nlines:
        line = lines[iline]
        stripped = line.strip()
        # log.info(f'B {iline}: {stripped}')
        if "SYMMETRIC (OR ASYMMETRIC) FINITE ELEMENT MODAL RESULTS ARE SUCCESSFULLY READ IN FROM FILE" in line:
            break
        if stripped == "":
            iline += 1
            break
        marker = " -       "
        idx = line.find(marker)
        if idx >= 0:
            card_line = line[idx + len(marker) :]
        else:
            card_line = stripped
        card_line = card_line.rstrip()
        if card_line.strip() == "ENDDATA":
            iline += 1
            break
        bulk_data_lines.append(card_line)
        iline += 1

    # --- parse modal eigenvalue table ---
    modal_data, iline = _parse_modal_table(lines, iline, nlines, log)
    # log.debug(f'end of modal table: line {iline}')

    # --- continue to flutter results ---
    # the subcase block has leading spaces before the asterisks

    # ****************************************
    # *                                      *
    # *       SUBCASE       =        1       *
    # *       DISCIPLINE    = FLUTTER        *
    # *       BULK ENTRY ID =      100       *
    # *                                      *
    # ****************************************
    # log.debug(f'iline = {iline}')
    while iline < nlines:
        line = lines[iline]
        if "                                             ****************************************" in line:
            log.debug(f'breaking on {line.rstrip()}')
            break
        # log.debug(line.rstrip())
        iline += 1

    iline += 2
    iline = cast(str, iline)
    subcase_sline = lines[iline].strip("\n *").split("=")
    subcase_sline = [val.strip() for val in subcase_sline]
    assert "SUBCASE" in subcase_sline, subcase_sline
    subcase = subcase_sline[1]
    assert subcase == "1", subcase_sline

    discipline_sline = lines[iline+1].strip("\n *").split("=")
    discipline_sline = [val.strip() for val in discipline_sline]
    assert "DISCIPLINE" in discipline_sline, discipline_sline
    discipline = discipline_sline[1]

    bulk_id_sline = lines[iline+2].strip("\n *").split("=")
    bulk_id_sline = [val.strip() for val in bulk_id_sline]
    assert "BULK ENTRY ID" in bulk_id_sline, bulk_id_sline
    bulk_id = bulk_id_sline[1]

    log.info(f"subcase={subcase!r} discipline={discipline!r} bulk_id={bulk_id!r}")
    iline += 3

    if discipline == "TRIM":
        log.info("TRIM discipline found; parsing trim results")
        trim_data = _parse_trim_section(lines, iline, nlines, log)
        trim_data['subcase'] = int(subcase)
        trim_data['bulk_id'] = int(bulk_id)
        return (
            out,
            {},
            exec_control_lines,
            case_control_lines,
            bulk_data_lines,
            modal_data,
            trim_data,
        )

    while "REFERENCE LENGTH (L) = " not in lines[iline]:
        #' MACH NUMBER =  0.9000, REFERENCE LENGTH (L) = 2.4000E+01/2.0 (IN), VREF= 1.0000E+00 (IN/SEC)'
        #' ALTITUDE = 0.0000E+00 (IN), ATMOS TABLE=STANDARD, REFERENCE LENGTH (L) = 1.0000E+01/2.0 (IN), VREF= 1.0000E+04 (IN/SEC)'
        iline += 1

    ref_line = lines[iline]
    ref_dict = split_ref_line(ref_line, log)
    log.debug(f"ref_dict = {ref_dict}")

    #' THE FOLLOWING V-G-F TABLE LISTS   49 NUMBER OF STRUCTURAL MODES (  1 - 49 ) AND    0 NUMBER OF AERODYNAMIC LAG ROOTS (  0 -  0 )'
    while "THE FOLLOWING V-G-F TABLE LISTS" not in lines[iline]:
        # print(f'{iline}: {lines[iline].rstrip()}')
        iline += 1

    units_line = ""
    imode = 0
    while iline < nlines:
        #   UNITS    UNITS    UNITS          MODE NO.  1              MODE NO.  2              MODE NO.  3              MODE NO.  4
        #    NONE   SLIN/     NONE       STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE
        #           IN**3               MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%
        #   V/VREF  DENSITY     MACH         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V
        #   0.0000 0.000+00   0.0000    0.0000   4.502 INFINT    0.0000   4.704 INFINT    0.0000   5.675 INFINT    0.0000   5.970 INFINT
        #   0.0100 1.146-07   0.0100    0.0000   4.502 8.7162    0.0000   4.704 9.1072    0.0000   5.67510.9868   -0.0004   5.96711.5524
        while "UNITS    UNITS    UNITS          MODE NO." not in lines[iline]:
            if "SUBCASE       =" in lines[iline]:
                log.warning("breaking on new subcase")
                return (out, ref_dict,
                        exec_control_lines, case_control_lines, bulk_data_lines,
                        modal_data, None,)
            if "***  Z A E R O   T E R M I N A T E D ***" in lines[iline]:
                log.info("***zero terminated***")
                return (out, ref_dict,
                        exec_control_lines, case_control_lines, bulk_data_lines,
                        modal_data, None, )
            # log.debug(f'Units {iline}: {lines[iline].rstrip()}')
            iline += 1

        modes = get_mode_sline(lines[iline])

        log.debug(f"modes = {modes}")
        for mode in modes:
            out[mode] = []
        iline += 1
        units_sline1 = lines[iline].strip("\n").split()
        units_sline2 = lines[iline + 1].strip("\n").split()
        names_sline = split_flutter_values(
            lines[iline + 2], apply_float=False
        )  # .strip('\n').split()
        log.debug(f"names_sline = {names_sline}")
        assert names_sline[0] == "V/VREF", names_sline

        density_units = ""
        header = "v/vref"
        if names_sline[1] == "V":
            assert units_sline2[0] == "SEC", units_sline2
            header += ",v"
            if units_sline1[1] == "IN/":
                velocity_units = "in/s"
            else:  # pragma: no cover
                raise RuntimeError((units_sline1[1], units_sline2[0]))
        elif names_sline[1] == "DENSITY":
            header += ",density"
            if units_sline1[1] == "SLIN/":
                assert units_sline2[0] == "IN**3", units_sline2
                density_units = "slinch/in^3"
            elif units_sline1[1] == "SLUG/":
                assert units_sline2[0] == "FT**3", units_sline2
                density_units = "slug/ft^3"
            else:  # pragma: no cover
                raise RuntimeError((units_sline1[1], units_sline2[0]))
        else:  # pragma: no cover
            raise NotImplementedError(
                f"names_sline[1]={names_sline[1]!r}; names_sline={names_sline}"
            )

        if names_sline[2] == "MACH":
            header += ",mach"
        elif names_sline[2] == "DYN P":  # 'DYN P'
            header += ",q"
        elif names_sline[2] == "ALTITUDE":  # 'DYN P'
            header += ",alt"
        else:  # pragma: no cover
            raise NotImplementedError(names_sline)

        # log.debug(f'units_sline1 = {units_sline1}')
        # log.debug(f'units_sline2 = {units_sline2}')
        # log.debug(f'names_sline = {names_sline}')
        assert "V/VREF" in names_sline, names_sline
        iline += 3
        if header not in out:
            if density_units:
                out["density_units"] = density_units
            out["header"] = header
            out[header] = []

        while lines[iline].strip() != "":
            values_str = lines[iline].rstrip()
            values = split_flutter_values(values_str)
            if modes[0] == 1:
                vvref_density_mach = values[:3]
                log.debug(f"*** {vvref_density_mach}")
                out[header].append(vvref_density_mach)
                # header
            other = values[3:]
            assert len(modes) * 3 == len(other), (modes, other)

            for imode, mode in enumerate(modes):
                other = values[(imode + 1) * 3 : (imode + 2) * 3]
                out[mode].append(other)
            iline += 1

        log.debug(f"header...{str(out[header])}")

        # print('modes:')
        nvalues = len(values)
        assert nvalues % 3 == 0, values
        dmode = (nvalues - 3) // 3
        assert dmode >= 1, dmode
        # print(f'dmode = {dmode}')

    log.info("------------------------------------------------")
    log.info(f"{iline}: {lines[iline].rstrip()}")

    raise RuntimeError("end of zaero_lines_to_out")

    #  THE FOLLOWING V-G-F TABLE LISTS   49 NUMBER OF STRUCTURAL MODES (  1 - 49 ) AND    0 NUMBER OF AERODYNAMIC LAG ROOTS (  0 -  0 )
    #
    #
    #                                            SUMMARY     MODES  1 -  4
    #
    # UNITS    UNITS    UNITS          MODE NO.  1              MODE NO.  2              MODE NO.  3              MODE NO.  4
    #  NONE   SLIN/     NONE       STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE
    #          IN**3               MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%
    #  V/VREF  DENSITY     MACH         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V
    #  0.0000 0.000+00   0.0000    0.0000   4.502 INFINT    0.0000   4.704 INFINT    0.0000   5.675 INFINT    0.0000   5.970 INFINT
    #  0.0100 1.146-07   0.0100    0.0000   4.502 8.7162    0.0000   4.704 9.1072    0.0000   5.67510.9868   -0.0004   5.96711.5524

    #    UNITS    UNITS    UNITS          MODE NO.145              MODE NO.146              MODE NO.147              MODE NO.148
    #     NONE   SLIN/   SLIN/        AERODYNAMIC LAG          AERODYNAMIC LAG          AERODYNAMIC LAG          AERODYNAMIC LAG
    #            IN**3   IN/S**2
    #    V/VREF  DENSITY    DYN P         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V
    #    0.0000 0.000+00 0.000+00    0.0000   0.000 INFINT    0.0000   0.000 INFINT    0.0000   0.000 INFINT    0.0000   0.000 INFINT
    #  12176.50 5.082-08 3.767+00   -2.4527   0.000 0.0000   -2.4527   0.000 0.0000  -184439.   0.000 0.0000  -24662.2   0.003 0.0000
    #  12702.30 7.284-08 5.876+00   -2.4526   0.000 0.0000   -2.4526   0.000 0.0000  -521511.   0.000 0.0000  0.000000   0.000 0.0000
    #  13151.00 9.788-08 8.464+00  -516766.   0.000 0.0000  ********   0.000 0.0000   -6.8127   0.000 0.0000   -4.3530   0.000 0.0000


def get_mode_sline(line: str) -> list[int]:
    """
    '   UNITS    UNITS    UNITS          MODE NO. 93              MODE NO. 94              MODE NO. 95              MODE NO. 96'
    -> [93, 94, 95, 96]
    '   UNITS    UNITS    UNITS          MODE NO.101              MODE NO.102              MODE NO.103              MODE NO.104'
    -> [101, 102, 103, 104]
    """
    mode_sline = line.strip("\n").split()
    mode_sline2 = []
    for val in mode_sline:
        val = val.rstrip(".")
        if "." in val:
            val_split = val.split(".")
            mode_sline2.extend(val_split)
        else:
            mode_sline2.append(val)

    # print('mode_sline2', mode_sline2)
    modes_strs = mode_sline2[5::3]
    modes = [int(mode) for mode in modes_strs]
    # log.debug(f'mode_sline = {mode_sline!r}')
    return modes


def _split_num_unit(value_str: str, log: SimpleLogger) -> tuple[float, str]:
    """
    Parameters
    ----------
    value_str

    Returns
    -------
    '1.0726E-07 (SLIN/IN**3  )'
    -> (1.0726E-07, 'SLIN/IN**3')

    """
    log.debug(value_str)
    assert "(" in value_str, value_str

    value_str2, unit = value_str.split(" ", 1)
    value_str2 = value_str2.strip()
    value = float(value_str2)
    unit = unit.strip(" ()").lower()
    return value, unit


def split_ref_line(line: str, log: SimpleLogger) -> dict[str, tuple[float, str]]:
    sline = line.strip().split(",")
    out = {}
    for name_value in sline:
        name, value_str = name_value.split("=")
        name = name.strip()
        value_str = value_str.strip()
        if name in {"MACH NUMBER", "MACH"}:
            value = float(value_str)
            unit = ""
        elif name == "VREF":
            value, unit = _split_num_unit(value_str, log)
            if unit == "in/sec":
                unit = "in/s"
            elif unit == "ft/sec":
                unit = "ft/s"
            assert unit in {"in/s", "ft/s", "m/s"}, unit
        elif name == "DENSITY":
            value, unit = _split_num_unit(value_str, log)

            if unit == "slin/in**3":
                unit = "slinch/in^3"
            # elif unit == 'ft/sec':
            #     unit = 'ft/s'
            assert unit in {"slinch/in^3"}, f"unit={unit!r}"
        elif name == "ALTITUDE":
            value, unit = _split_num_unit(value_str, log)
            assert unit in {"ft"}, f"unit={unit!r}"

        elif name == "REFERENCE LENGTH (L)":
            value_str2, unit = value_str.split(" ")
            value_str2 = value_str2.strip()
            assert value_str2.endswith("/2.0"), value_str2
            value_str3 = value_str2.split("/")[0]
            value = float(value_str3) / 2.0
            unit = unit.strip(" ()").lower()
            assert unit in {"in", "ft", "m"}, unit
        elif name == "ATMOS TABLE":
            if isinstance(value, str):
                assert value_str == "STANDARD", value_str
                value = value_str
            else:
                atmos_value = int(value_str)
                value = atmos_value
            unit = ""
        else:  # pragma: no cover
            raise RuntimeError(f"unhandled name; name={name!r} value={value_str!r}")
        out[name] = (value, unit)
    return out


def split_flutter_values(line: str, apply_float: bool = True) -> list[str]:
    """
    '   V/VREF  DENSITY     MACH         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V'
    '   0.0000 0.000+00   0.0000    0.0000   9.738 INFINT    0.0000  10.169 INFINT    0.0000  12.790 INFINT    0.0000  15.084 INFINT'

    Parameters
    ----------
    line

    Returns
    -------

    """
    values = [
        line[:9],  # 0: V/VREF
        line[9:18],  # 1: DENSITY
        line[18:27],  # 2: DYN P
        line[27:37],  # 3: G
        line[37:45],  # 4: F(HZ)
        line[45:52],  # 5: K=WL/V; infinit
        line[52:62],  # 6: G
        line[62:70],  # 7: F(HZ)
        line[70:77],  # 8: K=WL/V; infinit
        line[77:87],  # 9
        line[87:95],  # 10
        line[95:102],  # 11: infinit
        line[102:112],  # 12: G
        line[112:120],  # 13: F(HZ)
        line[120:],  # 14: K=WL/V
    ]
    # print(values)
    values2 = [value.strip() for value in values]
    if not apply_float:
        return values2

    # print(values2)
    values_out = []
    is_blank = False
    for i, value in enumerate(values2):
        assert " " not in value, f"{i}: {value}; values={values}"
        value = value.strip()
        if value == "" and is_blank:
            continue
        elif value == "" and not is_blank:
            assert is_blank is False, values2
            is_blank = True
            continue
        if value in {"INFINT", "+INFINT"}:
            value_out = np.inf  #'INFINT'
        elif "*" in value:
                             # '*******'
            assert value in {"*******", "********"}, f"i={i} value_out={value_out!r}; values2={values2}"
            value_out = np.nan
        else:
            try:
                value_out = cast_float(value)
            except RuntimeError:
                raise RuntimeError(f"{i}: {value}; values={values}")
            # value_out = float(value)
        values_out.append(value_out)
    return values_out


def cast_float(value_str: str):
    try:
        value = float(value_str)
    except ValueError:
        value_str = value_str.strip()
        if "+" in value_str[1:]:
            value_str2 = value_str[0] + value_str[1:].replace("+", "e+")
        elif "-" in value_str[1:]:
            value_str2 = value_str[0] + value_str[1:].replace("-", "e-")
        else:  # pragma: no cover
            raise RuntimeError(value_str)
        value = float(value_str2)
    return value
