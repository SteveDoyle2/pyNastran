import os
from typing import Optional
import numpy as np
from cpylog import SimpleLogger
from pyNastran.utils import PathLike, print_bad_path
from pyNastran.f06.flutter_response import FlutterResponse

def debug_to_level(debug: Optional[str | bool]) -> str:
    if debug is None:
        return 'warning'
    if isinstance(debug, str):
        return debug
    return 'debug' if debug else 'info'

def debug_to_log(debug: Optional[str|bool],
                 log: Optional[SimpleLogger]) -> SimpleLogger:
    level = debug_to_level(debug)
    log = SimpleLogger(level=level) if log is None else log
    return log

def read_zona_out(zona_out_filename: PathLike,
                  log: Optional[SimpleLogger]=None,
                  debug: Optional[str|bool]=True) -> tuple[dict, dict]:
    log = debug_to_log(debug, log)
    assert os.path.exists(zona_out_filename), print_bad_path(zona_out_filename)
    with open(zona_out_filename, 'r') as zona_out_file:
        lines = zona_out_file.readlines()

    case_dict, ref_dict = zona_lines_to_out(log, lines)
    for key, data in case_dict.items():
        log.debug(f'{key}:\n{str(data)}')

    # data = out_to_data()
    modes, result, in_units = out_dict_to_results(case_dict, ref_dict, log)

    # FlutterResponse.from_nx()
    assert isinstance(in_units, dict), in_units
    resp = FlutterResponse.from_zona(
        modes, result, in_units, zona_out_filename)
    responses = {
        1: resp,
    }
    mass = {}
    return responses, mass

def out_dict_to_results(out_dict: dict,
                        ref_dict: dict[str, tuple[float, str]],
                        log: SimpleLogger) -> tuple[list[int],
                                   np.ndarray,
                                   dict[str, str]]:
    header_name = out_dict['header']
    header_data_strs = out_dict[header_name]  # vvref,density,mach

    # v/vref, density, mach
    header_data = np.array(header_data_strs, dtype='float64')
    nrows = len(header_data)

    if header_name == 'v/vref,density,mach':
        vref, velocity_units = ref_dict['VREF']
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        mach = header_data[:, 2]
        density_units = out_dict['density_units']
    elif header_name == 'v/vref,density,q':
        vref, velocity_units = ref_dict['VREF']
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        q = header_data[:, 2]
        mach = q
        density_units = out_dict['density_units']
    elif header_name == 'v/vref,density,alt':
        # ref_dict = {'MACH': (0.8, ''), 'ATMOS TABLE': ('STANDARD', ''),
        #             'REFERENCE LENGTH (L)': (1.0, 'in'),
        #             'VREF': (2.0, 'in/s')}
        #print('ref_dict = ', ref_dict)
        vref, velocity_units = ref_dict['VREF']
        atmos_table = ref_dict['ATMOS TABLE'][0]
        assert atmos_table == 'STANDARD', atmos_table
        mach = ref_dict['MACH'][0]
        #alt, altitude_units = ref_dict['ALT']
        #print(f'vref={vref}; velocity_units={velocity_units}')
        velocity = header_data[:, 0] * vref
        density = header_data[:, 1]
        q = 0.5 * density * velocity ** 2
        density_units = out_dict['density_units']
    elif header_name == 'v/vref,v,q':
        #print('ref_dict = ', ref_dict)
        #velocity0, velocity_units = ref_dict['VELOCITY']
        #atmos_table = ref_dict['ATMOS TABLE'][0]
        #assert atmos_table == 'STANDARD', atmos_table
        #mach = ref_dict['MACH'][0]
        # alt, altitude_units = ref_dict['ALT']
        #print(f'vref={vref}; velocity_units={velocity_units}')
        v_vref = header_data[:, 0]
        #print((velocity - velocity0).max())
        q = header_data[:, 1]
        if 'DENSITY' in ref_dict and 'VREF' in ref_dict and 'MACH' in ref_dict:
            vref, velocity_units = ref_dict['VREF']
            velocity = v_vref * vref
            density, density_units = ref_dict['DENSITY']
            mach = ref_dict['MACH'][0]
        else:
            raise RuntimeError(f'V,q; {ref_dict}')
        log.debug(f'density = {density}')
    else:  # pragma: no cover
        raise NotImplementedError(header_name)

    if velocity_units in 'in/s':
        altitude_units = 'ft'
        dynamic_pressure_units = 'psi'
    elif velocity_units in 'ft/s':
        altitude_units = 'ft'
        dynamic_pressure_units = 'psf'
    else:
        assert velocity_units == 'm/s', velocity_units
        altitude_units = 'm'
        dynamic_pressure_units = 'Pa'

    #print(f'out_dict = {list(out_dict.keys())}')
    #print(f'out_dict[header] = {out_dict['header']}')
    in_units = {
        'velocity': velocity_units,
        'density': density_units,
        'eas': velocity_units,
        'altitude': altitude_units,
        'dynamic_pressure': dynamic_pressure_units,
    }

    modes = [mode for mode in out_dict if isinstance(mode, int)]
    nmodes = len(modes)
    log.info(f'modes = {modes}')
    results_list = []
    result = np.zeros((nmodes, nrows, 9), dtype='float64')
    log.debug(f'nmodes={nmodes}; nrows={nrows}')
    for imode, mode in enumerate(modes):
        # g, f(Hz), k
        data = out_dict[mode]
        g_freq_k = np.array(data, dtype='float64')
        damping = g_freq_k[:, 0]
        freq = g_freq_k[:, 1]
        kfreq = g_freq_k[:, 2]

        omega = 2 * np.pi * freq
        eigr = omega * damping / 2
        eigi = omega

        # fdata   [                rho,     mach, velocity, damping, freq]
        # results [kfreq, 1/kfreq, density, mach, velocity, damping, freq, eigr, eigi]
        #log.debug(len(kfreq))
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

def zona_lines_to_out(log: SimpleLogger, lines: list[str]) -> tuple[dict, dict]:
    out = {}
    nlines = len(lines)
    iline = 0
    #     SYMMETRIC (OR ASYMMETRIC) FINITE ELEMENT MODAL RESULTS ARE SUCCESSFULLY READ IN FROM FILE junk.f06
    # RIGID BODY DEGREES OF FREEDOM (DEFINED IN THE FEM BASIC COORDINATE SYSTEM) =       0
    #
    #        MODE      EXTRACTION      EIGENVALUE              FREQUENCY                      GENERALIZED
    #                    ORDER         (RAD/S)**2         (RAD/S)            (HZ)            MASS       STIFFNESS
    #           1               1    -0.12395E-05     0.11133E-02     0.17719E-03     0.10000E+01     0.12395E-05
    #           2               2    -0.34338E-06     0.58599E-03     0.93263E-04     0.10000E+01     0.34338E-06

    # ****************************************
    # *                                      *
    # *       SUBCASE       =        1       *
    # *       DISCIPLINE    = FLUTTER        *
    # *       BULK ENTRY ID =  1000001       *
    # *                                      *
    # ****************************************

    while '                                             ****************************************' not in lines[iline]:
        #print(f'{iline}: {lines[iline].rstrip()}')
        iline += 1

    iline += 2
    subcase_sline = lines[iline].strip('\n *').split('=')
    subcase_sline = [val.strip() for val in subcase_sline]
    assert 'SUBCASE' in subcase_sline, subcase_sline
    subcase = subcase_sline[1]
    assert subcase == '1', subcase_sline

    discipline_sline = lines[iline+1].strip('\n *').split('=')
    discipline_sline = [val.strip() for val in discipline_sline]
    assert 'DISCIPLINE' in discipline_sline, discipline_sline
    discipline = discipline_sline[1]

    bulk_id_sline = lines[iline+2].strip('\n *').split('=')
    bulk_id_sline = [val.strip() for val in bulk_id_sline]
    assert 'BULK ENTRY ID' in bulk_id_sline, bulk_id_sline
    bulk_id = bulk_id_sline[1]

    log.info(f'subcase={subcase!r} discipline={discipline!r} bulk_id={bulk_id!r}')
    iline += 3

    while 'REFERENCE LENGTH (L) = ' not in lines[iline]:
        #' MACH NUMBER =  0.9000, REFERENCE LENGTH (L) = 2.4000E+01/2.0 (IN), VREF= 1.0000E+00 (IN/SEC)'
        #' ALTITUDE = 0.0000E+00 (IN), ATMOS TABLE=STANDARD, REFERENCE LENGTH (L) = 1.0000E+01/2.0 (IN), VREF= 1.0000E+04 (IN/SEC)'
        iline += 1

    ref_line = lines[iline]
    ref_dict = split_ref_line(ref_line)
    log.debug(f'ref_dict = {ref_dict}')

    #' THE FOLLOWING V-G-F TABLE LISTS   49 NUMBER OF STRUCTURAL MODES (  1 - 49 ) AND    0 NUMBER OF AERODYNAMIC LAG ROOTS (  0 -  0 )'
    while 'THE FOLLOWING V-G-F TABLE LISTS' not in lines[iline]:
        #print(f'{iline}: {lines[iline].rstrip()}')
        iline += 1

    units_line = ''
    imode = 0
    while iline < nlines:
        #   UNITS    UNITS    UNITS          MODE NO.  1              MODE NO.  2              MODE NO.  3              MODE NO.  4
        #    NONE   SLIN/     NONE       STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE          STRUCTURAL MODE
        #           IN**3               MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%     MODAL DAMPING= 0.00%
        #   V/VREF  DENSITY     MACH         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V         G   F(HZ) K=WL/V
        #   0.0000 0.000+00   0.0000    0.0000   4.502 INFINT    0.0000   4.704 INFINT    0.0000   5.675 INFINT    0.0000   5.970 INFINT
        #   0.0100 1.146-07   0.0100    0.0000   4.502 8.7162    0.0000   4.704 9.1072    0.0000   5.67510.9868   -0.0004   5.96711.5524
        while 'UNITS    UNITS    UNITS          MODE NO.' not in lines[iline]:
            if 'SUBCASE       =' in lines[iline]:
                log.warning('breaking on new subcase')
                return out, ref_dict
            if '***  Z A E R O   T E R M I N A T E D ***' in lines[iline]:
                log.info('***zero terminated***')
                return out, ref_dict
            #log.debug(f'Units {iline}: {lines[iline].rstrip()}')
            iline += 1

        modes = get_mode_sline(lines[iline])

        log.debug(f'modes = {modes}')
        for mode in modes:
            out[mode] = []
        iline += 1
        units_sline1 = lines[iline].strip('\n').split()
        units_sline2 = lines[iline+1].strip('\n').split()
        names_sline = split_flutter_values(
            lines[iline+2], apply_float=False)  #.strip('\n').split()
        log.debug(f'names_sline = {names_sline}')
        assert names_sline[0] == 'V/VREF', names_sline

        density_units = ''
        header = 'v/vref'
        if names_sline[1] == 'V':
            assert units_sline2[0] == 'SEC', units_sline2
            header += ',v'
            if units_sline1[1] == 'IN/':
                velocity_units = 'in/s'
            else:  # pragma: no cover
                raise RuntimeError((units_sline1[1], units_sline2[0]))
        elif names_sline[1] == 'DENSITY':
            header += ',density'
            if units_sline1[1] == 'SLIN/':
                assert units_sline2[0] == 'IN**3', units_sline2
                density_units = 'slinch/in^3'
            elif units_sline1[1] == 'SLUG/':
                assert units_sline2[0] == 'FT**3', units_sline2
                density_units = 'slug/ft^3'
            else:  # pragma: no cover
                raise RuntimeError((units_sline1[1], units_sline2[0]))
        else:  # pragma: no cover
            raise NotImplementedError(f'names_sline[1]={names_sline[1]!r}; names_sline={names_sline}')

        if names_sline[2] == 'MACH':
            header += ',mach'
        elif names_sline[2] == 'DYN P':  # 'DYN P'
            header += ',q'
        elif names_sline[2] == 'ALTITUDE':  # 'DYN P'
            header += ',alt'
        else:  # pragma: no cover
            raise NotImplementedError(names_sline)

        #log.debug(f'units_sline1 = {units_sline1}')
        #log.debug(f'units_sline2 = {units_sline2}')
        #log.debug(f'names_sline = {names_sline}')
        assert 'V/VREF' in names_sline, names_sline
        iline += 3
        if header not in out:
            if density_units:
                out['density_units'] = density_units
            out['header'] = header
            out[header] = []

        while lines[iline].strip() != '':
            values_str = lines[iline].rstrip()
            values = split_flutter_values(values_str)
            if modes[0] == 1:
                vvref_density_mach = values[:3]
                log.debug(f'*** {vvref_density_mach}')
                out[header].append(vvref_density_mach)
                # header
            other = values[3:]
            assert len(modes)*3 == len(other), (modes, other)

            for imode, mode in enumerate(modes):
                other = values[(imode+1)*3:(imode+2)*3]
                out[mode].append(other)
            iline += 1

        log.debug(f'header...{str(out[header])}')

        #print('modes:')
        nvalues = len(values)
        assert nvalues % 3 == 0, values
        dmode = (nvalues - 3) // 3
        assert dmode >= 1, dmode
        #print(f'dmode = {dmode}')

    log.info('------------------------------------------------')
    log.info(f'{iline}: {lines[iline].rstrip()}')

    raise RuntimeError('end of zona_lines_to_out')

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

def get_mode_sline(line: str) -> list[int]:
    """
    '   UNITS    UNITS    UNITS          MODE NO. 93              MODE NO. 94              MODE NO. 95              MODE NO. 96'
    -> [93, 94, 95, 96]
    '   UNITS    UNITS    UNITS          MODE NO.101              MODE NO.102              MODE NO.103              MODE NO.104'
    -> [101, 102, 103, 104]
    """
    mode_sline = line.strip('\n').split()
    mode_sline2 = []
    for val in mode_sline:
        val = val.rstrip('.')
        if '.' in val:
            val_split = val.split('.')
            mode_sline2.extend(val_split)
        else:
            mode_sline2.append(val)

    #print('mode_sline2', mode_sline2)
    modes_strs = mode_sline2[5::3]
    modes = [int(mode) for mode in modes_strs]
    #log.debug(f'mode_sline = {mode_sline!r}')
    return modes


def _split_num_unit(value_str: str) -> tuple[float, str]:
    """
    Parameters
    ----------
    value_str

    Returns
    -------

    '1.0726E-07 (SLIN/IN**3  )'
    -> (1.0726E-07, 'SLIN/IN**3')

    """
    print(value_str)
    assert '(' in value_str, value_str

    value_str2, unit = value_str.split(' ', 1)
    value_str2 = value_str2.strip()
    value = float(value_str2)
    unit = unit.strip(' ()').lower()
    return value, unit

def split_ref_line(line: str) -> dict[str, tuple[float, str]]:
    sline = line.strip().split(',')
    out = {}
    for name_value in sline:
        name, value_str = name_value.split('=')
        name = name.strip()
        value_str = value_str.strip()
        if name in {'MACH NUMBER', 'MACH'}:
            value = float(value_str)
            unit = ''
        elif name == 'VREF':
            value, unit = _split_num_unit(value_str)
            if unit == 'in/sec':
                unit = 'in/s'
            elif unit == 'ft/sec':
                unit = 'ft/s'
            assert unit in {'in/s', 'ft/s', 'm/s'}, unit
        elif name == 'DENSITY':
            value, unit = _split_num_unit(value_str)

            if unit == 'slin/in**3':
                unit = 'slinch/in^3'
            # elif unit == 'ft/sec':
            #     unit = 'ft/s'
            assert unit in {'slinch/in^3'}, f'unit={unit!r}'
        elif name == 'ALTITUDE':
            value, unit = _split_num_unit(value_str)
            assert unit in {'ft'}, f'unit={unit!r}'

        elif name == 'REFERENCE LENGTH (L)':
            value_str2, unit = value_str.split(' ')
            value_str2 = value_str2.strip()
            assert value_str2.endswith('/2.0'), value_str2
            value_str3 = value_str2.split('/')[0]
            value = float(value_str3) / 2.0
            unit = unit.strip(' ()').lower()
            assert unit in {'in', 'ft', 'm'}, unit
        elif name == 'ATMOS TABLE':
            assert value_str == 'STANDARD', value_str
            value = value_str
            unit = ''
        else:  # pragma: no cover
            raise RuntimeError(f'unhandled name; name={name!r} value={value_str!r}')
        out[name] = (value, unit)
    return out

def split_flutter_values(line: str,
                         apply_float: bool=True) -> list[str]:
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
        line[:9],     # 0: V/VREF
        line[9:18],   # 1: DENSITY
        line[18:27],  # 2: DYN P
        line[27:37],  # 3: G
        line[37:45],  # 4: F(HZ)
        line[45:52],  # 5: K=WL/V; infinit

        line[52:62],  # 6: G
        line[62:70],  # 7: F(HZ)
        line[70:77],  # 8: K=WL/V; infinit

        line[77:87],  # 9
        line[87:95],  # 10
        line[95:102],   # 11: infinit
        line[102:112],  # 12: G
        line[112:120],  # 13: F(HZ)
        line[120:],     # 14: K=WL/V
    ]
    #print(values)
    values2 = [value.strip() for value in values]
    if not apply_float:
        return values2

    #print(values2)
    values_out = []
    is_blank = False
    for i, value in enumerate(values2):
        assert ' ' not in value, f'{i}: {value}; values={values}'
        value = value.strip()
        if value == '' and is_blank:
            continue
        elif value == '' and not is_blank:
            assert is_blank is False, values2
            is_blank = True
            continue
        if value in {'INFINT', '+INFINT'}:
            value_out = np.inf #'INFINT'
        else:
            try:
                value_out = cast_float(value)
            except RuntimeError:
                raise RuntimeError(f'{i}: {value}; values={values}')
            #value_out = float(value)
        values_out.append(value_out)
    return values_out

def cast_float(value_str: str):
    try:
        value = float(value_str)
    except ValueError:
        value_str = value_str.strip()
        if '+' in value_str[1:]:
            value_str2 = value_str[0] + value_str[1:].replace('+', 'e+')
        elif '-' in value_str[1:]:
            value_str2 = value_str[0] + value_str[1:].replace('-', 'e-')
        else:  # pragma: no cover
            raise RuntimeError(value_str)
        value = float(value_str2)
    return value
