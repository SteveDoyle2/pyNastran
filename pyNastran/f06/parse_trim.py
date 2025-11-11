from typing import Optional, TextIO, Any
import os
import numpy as np
#import scipy.sparse
from cpylog import SimpleLogger  # get_logger
from pyNastran.utils import print_bad_path
from pyNastran.f06.f06_tables.trim import (
    MonitorLoads, TrimResults, ControllerState,
    AeroPressure, AeroForce, TrimVariables, TrimVariable)

from cpylog import __version__ as CPYLOG_VERSION
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger


#'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'

MONTHS = {'JANUARY', 'FEBRUARY', 'MARCH', 'APRIL', 'MAY', 'JUNE',
          'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER'}
SKIP_FLAGS = [
    'This software and related documentation are',
    'LIMITATIONS TO U.S. GOVERNMENT RIGHTS. UNPUBLISHED',
    'TOLERANCE LIMITS ARE',
    'SUPPORT PT.NO.             EPSILON             STRAIN   ENERGY',
    '*** USER WARNING MESSAGE 2020 (SDRCOMP)',

    'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O',
    'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O',
    'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N',
    'M O D E L   S U M M A R Y',
    'C A S E    C O N T R O L    E C H O',
    'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y',
    'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R',

    # SOL 144
    'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S',
    #'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S',
    'A E R O D Y N A M I C   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S',
    #'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S',
    #'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (READ MODULE)',
    #'R E A L   E I G E N V A L U E S',
    'MAXIMUM  DISPLACEMENTS',
    'FLUTTER  SUMMARY',
    #'* * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *',  # causes a crash
]


def read_f06_trim(f06_filename: str,
                  log: Optional[SimpleLogger]=None,
                  nlines_max: int=1_000_000,
                  debug: bool=False) -> dict[str, TrimResults]:
    """TODO: doesn't handle extra PAGE headers; requires LINE=1000000"""
    log = get_logger(log, debug, encoding='utf-8')
    dirname = os.path.dirname(os.path.abspath(f06_filename))
    assert os.path.exists(f06_filename), print_bad_path(f06_filename)
    log.info(f'reading {f06_filename!r}')
    with open(f06_filename, 'r') as f06_file:
        trim_results, tables, matrices = _read_f06_trim(
            f06_file, log, nlines_max, dirname,
            debug=debug)

    stacked_aero_pressure = _stack_data_dict(trim_results.aero_pressure, 2)
    stacked_aero_force = _stack_data_dict(trim_results.aero_force, 3)

    for subcase, aero_pressure in stacked_aero_pressure.items():
        metadata = trim_results.metadata[subcase]
        trim_results.aero_pressure[subcase] = AeroPressure.from_f06(
            subcase, *aero_pressure, metadata)

    for subcase, aero_force in stacked_aero_force.items():
        metadata = trim_results.metadata[subcase]
        trim_results.aero_force[subcase] = AeroForce.from_f06(
            subcase, *aero_force, metadata)

    out = {
        'trim_results': trim_results,
    }
    if len(tables):
        out['tables'] = tables
        log.info('found the following tables in the f06: %s' % (list(tables)))
    if len(matrices):
        out['matrices'] = matrices
        log.info('found the following matrices in the f06: %s' % (list(matrices)))
    str(trim_results)
    return out


def _stack_data_dict(results: dict[np.ndarray, np.ndarray],
                     narray: int) -> dict[np.ndarray, np.ndarray]:
    results_out = {}
    for subcase, press_data in results.items():
        all_nids = []
        all_eids = []
        all_data = []
        all_labels = []
        if narray == 3:
            for nid, datai, labeli in press_data:
                all_nids.append(nid)
                all_data.append(datai)
                all_labels.extend(labeli)
            nids2 = np.hstack(all_nids)
            data2 = np.vstack(all_data)
            labels = np.array(all_labels)
            results_out[subcase] = (nids2, data2, labels)
        else:
            for nid, eid, datai in press_data:
                all_nids.append(nid)
                all_eids.append(eid)
                all_data.append(datai)
            nids2 = np.hstack(all_nids)
            eids2 = np.hstack(all_eids)
            data2 = np.vstack(all_data)
            results_out[subcase] = (nids2, eids2, data2)
    return results_out


def _skip_to_page_stamp_and_rewind(f06_file: TextIO, line: str, i: int,
                                   nlines_max: int) -> tuple[str, int, int]:
    seek0 = f06_file.tell()
    line_end, iend = _skip_to_page_stamp(f06_file, line, i, nlines_max)
    seek1 = f06_file.tell()

    f06_file.seek(seek0)
    return line_end, iend, seek1


def _skip_to_page_stamp(f06_file: TextIO, line: str, i: int,
                        nlines_max: int) -> tuple[str, int]:
    line = f06_file.readline()
    i += 1
    # JANUARY  26, 2012  SIMCENTER NASTRAN  3/12/20   PAGE     2
    while 'NASTRAN' not in line and 'PAGE' not in line:
        line = f06_file.readline()
        i += 1
        if i > nlines_max:
            raise RuntimeError(f'{nlines_max:d} lines in file is max?...\n'
                               'this will be removed once the parser is better tested')
    #log.debug(f'line = {line.rstrip()}')
    return line, i


def _read_f06_trim(f06_file: TextIO, log: SimpleLogger,
                   nlines_max: int, dirname: str,
                   debug: bool=False) -> tuple[TrimResults,
                                               dict[str, np.ndarray],
                                               dict[str, np.ndarray]]:
    i = 0
    #debug = True
    tables: dict[str, np.ndarray] = {}
    matrices: dict[str, np.ndarray] = {}
    trim_results = TrimResults()
    iblank_count = 0

    ipressure = 0
    iforce = 0
    while True:
        line = f06_file.readline()
        i += 1
        if debug:
            log.debug(f'i={i} {line.strip()}')
        if '* * * END OF JOB * * *' in line:
            #print("****done****")
            break
        iflags = [datai in line for datai in SKIP_FLAGS]
        if any(iflags):
            iblank_count = 0
            flag = SKIP_FLAGS[iflags.index(True)]
            if debug:
                log.info(f'******* found skip flag: {flag}')
            #print('skip', line)
            line, i = _skip_to_page_stamp(f06_file, line, i, nlines_max)
            if 'trademark' not in line:
                line, i, title, subtitle, subcase = _get_title_subtitle_subcase(f06_file, line, i, nlines_max)
            continue
        elif line.startswith('0      MATRIX '):
            iblank_count = 0
            line, i = _skip_to_page_stamp(f06_file, line, i, nlines_max)
            line, i, title, subtitle, subcase = _get_title_subtitle_subcase(f06_file, line, i, nlines_max)
            #matrix_name, matrix, line, i = _read_matrix(f06_file, line, i, log, debug)
            #matrices[matrix_name] = matrix
            #del matrix_name, matrix
        elif 'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S' in line:
            log.debug('reading aero static data recovery tables')
            iblank_count = 0
            line, i, ipressure, iforce, metadata = _read_aerostatic_data_recovery_output_table(
                f06_file, line, i, nlines_max,
                trim_results,
                title, subtitle, subcase,
                ipressure, iforce, log)
        elif 'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' in line:
            log.debug('reading aero static data recovery tables')
            iblank_count = 0
            line, i = _read_structural_monitor_point_integrated_loads(
                f06_file, line, i, nlines_max, trim_results,
                title, subtitle, subcase, log)
        elif 'PAGE' in line and any(month in line for month in MONTHS):
            line, i, title, subtitle, subcase = _get_title_subtitle_subcase(
                f06_file, line, i, nlines_max)
            #log.info(f'title={title!r} subtitle={subtitle!r}')
        else:
            #log.debug(f'else: i={i} {line.strip()}')
            line_strip = line.strip()
            if len(line_strip) == 0:
                iblank_count += 1
            else:
                iblank_count = 0
            #print(line)
        #print('----')

        if iblank_count == 1000:
            log.warning('breaking because 1000 blank lines were found; assuming theres an error (or incomplete deck)')
            break
        if i > nlines_max:
            raise RuntimeError(f'{nlines_max:d} lines in file is max?...\n'
                               'this will be removed once the parser is better tested')
        if i % 1000 == 0:
            log.debug(f'i={i}')
    return trim_results, tables, matrices


def _get_title_subtitle_subcase(f06_file: TextIO,
                                line: str, i: int,
                                nlines_max: int) -> tuple[str, int,
                                                          str, str, str]:
    """
    1    144                                                                   FEBRUARY   1, 2021  SIMCENTER NASTRAN  3/12/20   PAGE    10
      SUBTITLE
    0                                                                                                            SUBCASE 42
    """
    n = f06_file.tell()
    title = ''
    subtitle_line = f06_file.readline()
    subcase_line = f06_file.readline()

    #'SUBCASE 1'
    if 'SUBCASE' not in subcase_line:
        f06_file.seek(n)
        return line, i, title, '', ''

    i += 2
    assert 'SUBCASE' in subcase_line, f'i={i:d} subcase_line={subcase_line!r}'
    subtitle = subtitle_line.strip()
    #'0                                                                                                            SUBCASE 13             '
    subcase = subcase_line[1:].strip()
    subcase_int = subcase.split('SUBCASE ', 1)[1].strip()
    assert ' ' not in subcase_int, f'subcase_int={subcase_int!r}; subcase={subcase!r}'
    assert ',' not in subcase_int, f'subcase_int={subcase_int!r}; subcase={subcase!r}'
    assert len(subcase_int) < 4, f'len(subcase_int)={len(subcase_int)}; subcase_int={subcase_int!r}; subcase={subcase!r}'
    return line, i, title, subtitle, subcase_int


def _read_structural_monitor_point_integrated_loads(f06_file: TextIO,
                                                    line: str, i: int, nlines_max: int,
                                                    trim_results: TrimResults,
                                                    title: str, subtitle: str, subcase: str,
                                                    #ipressure: int, iforce: int,
                                                    log: SimpleLogger) -> tuple[str, int]:
    """
    '                              S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S'
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC'
    '                                           MACH = 1.000000E-01                Q = 1.587000E-01'
    '
    '        CONTROLLER STATE:'
    '        ANGLEA   =   1.0000E-01'
    '
    '        MONITOR POINT NAME = AEROSG2D          COMPONENT =                   CLASS = COEFFICIENT               '
    '        LABEL = Full Vehicle Integrated Loads                           '
    '        CID =      102          X =  0.00000E+00          Y =  0.00000E+00          Z =  0.00000E+00'
    '
    '          AXIS      RIGID AIR       ELASTIC REST.   RIGID APPLIED    REST. APPLIED   '
    '          ----    -------------    -------------    -------------    -------------   '
    '           CX     0.000000E+00     0.000000E+00     0.000000E+00     0.000000E+00'
    '           CY    -1.062477E+01    -1.062477E+01     0.000000E+00     0.000000E+00'
    '           CZ     1.382605E+02     1.382605E+02     0.000000E+00     0.000000E+00'
    '           CMX    3.801091E+03     3.801091E+03     0.000000E+00     0.000000E+00'
    '           CMY   -1.767277E+03    -1.767277E+03     0.000000E+00     0.000000E+00'
    '           CMZ   -1.948738E+02    -1.948738E+02     0.000000E+00     0.000000E+00'
    ''
    '        MONITOR POINT NAME = AE01              COMPONENT = AE01              CLASS = GENERAL                   '
    """
    header_lines = []
    i0 = i

    while 'MONITOR POINT NAME' not in line and ('NASTRAN' not in line and 'PAGE' not in line):
        if i != i0:
            header_lines.append(line)
        line = f06_file.readline()
        i += 1
        #line, i = _read_structural_monitor_point_integrated_loads(
        #f06_file, line, i, nlines_max,
        #title, subtitle, subcase,
        #log)

    if 'MONITOR POINT NAME' not in line:
        raise RuntimeError("'MONITOR POINT NAME' not in line")

    controller_state = _get_controller_state(header_lines)
    line_end, iend, seek1 = _skip_to_page_stamp_and_rewind(f06_file, line, i, nlines_max)

    names = []
    comps = []
    classes = []
    labels = []

    #name_comps_classes_labels = []
    xyzs = []
    cids = []
    all_coeffs = []
    cds = []
    axis_to_index = {
        'CX': 0, 'CY': 1, 'CZ': 2,
        'CMX': 3, 'CMY': 4, 'CMZ': 5, }

    while i < iend and 'MONITOR POINT NAME' in line:
        #'        MONITOR POINT NAME = AEROSG2D          COMPONENT =                   CLASS = COEFFICIENT               '
        #'        LABEL = Full Vehicle Integrated Loads                           '
        #'        CID =      102          X =  0.00000E+00          Y =  0.00000E+00          Z =  0.00000E+00'

        #line = _remove_intermediate_spaces(line)
        #line2 = _remove_intermediate_spaces(line.replace('MONITOR POINT NAME', 'MONITOR_POINT_NAME'))
        #sline2 = line2.split(' ')
        #assert len(line2) == 3, sline2

        name_comp_class = line.split('MONITOR POINT NAME =')[1]
        name_comp, classi = name_comp_class.rsplit('CLASS = ')
        name_comp = name_comp.strip()

        # COEFFICIENT: summation about the CG?
        # GENERAL:     MONPNT1
        classi = classi.strip()
        assert classi in ['COEFFICIENT', 'GENERAL'], classi
        name, comp = name_comp.split('COMPONENT =')
        name = name.strip()
        comp = comp.strip()
        #print(f'name={name!r} comp={comp!r} class={classi!r}')

        #'        LABEL = Full Vehicle Integrated Loads                           '
        line = f06_file.readline()
        i += 1
        label = line.split('LABEL =')[1].strip()

        line = f06_file.readline()
        line = _remove_intermediate_spaces(line)
        i += 1

        sline = line.split(' ')
        if len(sline) == 4:
            # NX
            assert 'CID' in sline[0], sline
            cd = -1
        else:
            # TODO: doesn't support CD
            assert len(sline) == 5, sline
            # MSC
            assert 'CP' in sline[0], sline
            assert 'CD' in sline[4], sline
            #['CP=100', 'X=0.00000E+00', 'Y=0.00000E+00', 'Z=0.00000E+00', 'CD=100']
            sline_cd = sline[4].split('=')
            assert sline_cd[0] == 'CD', sline
            cd = int(sline_cd[1])

        sline_cid_cd = sline[0].split('=')
        sline_x = sline[1].split('=')
        sline_y = sline[2].split('=')
        sline_z = sline[3].split('=')
        assert sline_x[0] == 'X', sline
        assert sline_y[0] == 'Y', sline
        assert sline_z[0] == 'Z', sline

        cid = int(sline_cid_cd[1])
        xyz = [
            float(sline_x[1]),
            float(sline_x[1]),
            float(sline_x[1]),
        ]

        unused_line0a = f06_file.readline()
        line0b = f06_file.readline().strip()
        i += 2

        if line0b == 'AXIS      RIGID AIR       ELASTIC REST.   RIGID APPLIED    REST. APPLIED':
            #AXIS      RIGID AIR       ELASTIC REST.   RIGID APPLIED    REST. APPLIED
            line = f06_file.readline()
            i += 1
            coeffs = np.zeros((6, 4), dtype='float64')

            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CX', axis
            coeffs[0, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CY', axis
            coeffs[1, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CZ', axis
            coeffs[2, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            # -----
            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CMX', axis
            coeffs[3, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CMY', axis
            coeffs[4, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            line = f06_file.readline()
            i += 1
            axis, rigid_air, elastic_rest, rigid_applied, rest_applied = line.split()
            assert axis == 'CMZ', axis
            coeffs[5, :] = [rigid_air, elastic_rest, rigid_applied, rest_applied]

            line = f06_file.readline()
            line = f06_file.readline()
            i += 2

        elif line0b == 'AXIS      RIGID AIR       ELASTIC REST.      INERTIAL      RIGID APPLIED    REST. APPLIED':
            line = f06_file.readline()
            i += 1
            coeffs = np.zeros((6, 5), dtype='float64')

            while i < iend - 1:
                line = f06_file.readline().strip()
                i += 1
                if len(line) == 0:  # or 'PAGE' in line or 'NASTRAN' in line:
                    break
                axis, rigid_air, elastic_rest, inertial, rigid_applied, rest_applied = line.split()
                iaxis = axis_to_index[axis]
                coeffs[iaxis, :] = [rigid_air, elastic_rest, inertial, rigid_applied, rest_applied]
            line = f06_file.readline()
            i += 1

        names.append(name)
        comps.append(comp)
        classes.append(classi)
        labels.append(label)
        #name_comps_classes_labels.append([name, comp, classi, label])
        cids.append(cid)
        xyzs.append(xyz)
        cds.append(cd)
        all_coeffs.append(coeffs)
    all_coeffs_stacked = np.stack(all_coeffs, axis=0)
    nrows = len(names)
    assert all_coeffs_stacked.shape == (nrows, 6, 5), all_coeffs_stacked.shape

    names_array = np.array(names)
    comps_array = np.array(comps)
    classes_array = np.array(classes)
    labels_array = np.array(labels)
    cids_array = np.array(cids, dtype='int32')
    xyzs_array = np.array(xyzs, dtype='float64')
    all_coeffs_array = np.array(all_coeffs_stacked, dtype='float64')
    cds_array = np.array(cds, dtype='int32')

    isubcase = int(subcase)
    trim_results.structural_monitor_loads[isubcase] = MonitorLoads(
        names_array, comps_array, classes_array, labels_array,
        cids_array, xyzs_array, all_coeffs_array, cds_array)
    trim_results.controller_state[isubcase] = controller_state
    f06_file.seek(seek1)
    return line_end, iend


def _remove_intermediate_spaces(line_in: str) -> str:
    """
    Simplifies parsing of the following line

    Parameters
    ----------
    line_in : str
        a line that looks like:
        '        CID =      102          X =  0.00000E+00          Y =  0.00000E+00          Z =  0.00000E+00'

    Returns
    -------
    line : str
        a simple to split line
        'CID=102 X=0.00000E+00 Y=0.00000E+00 Z=0.00000E+00'
    """
    line = line_in.strip()
    while ' =' in line or '= ' in line or '  ' in line:
        line = line.replace(' =', '=').replace('= ', '=').replace('  ', ' ')
    return line


def _get_controller_state(header_lines: list[str]) -> ControllerState:
    controller_state: ControllerState = {}
    controller_lines: list[str] = []
    for i, line in enumerate(header_lines):
        if 'CONTROLLER STATE:' in line:
            controller_lines = header_lines[i+1:]
            break
    for line in controller_lines:
        line2 = line.strip()
        if not line2:
            continue

        #name, value = line2.split('=')
        sline = line2.replace('=', '').strip().split()
        nsline = len(sline)
        assert nsline % 2 == 0, sline
        for j in range(0, nsline, 2):
            name = sline[j]
            value = sline[j+1]
            controller_state[name] = float(value)
    return controller_state


def _read_aeroelastic_trim_variables(f06_file: TextIO,
                                     line: str, i: int, nlines_max: int,
                                     trim_results: TrimResults, isubcase: int,
                                     metadata: dict[str, Any]) -> tuple[str, int]:
    """
    '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = SYMMETRIC      XZ-SYMMETRY = SYMMETRIC'
    '                                           MACH = 1.000000E-01                Q = 1.000000E-01'
    '                         CHORD = 1.0000E+01           SPAN = 1.0000E+02            AREA = 1.0000E+03'
    ''
    '          TRIM ALGORITHM USED: LINEAR TRIM SOLUTION WITHOUT REDUNDANT CONTROL SURFACES.'       <---------- you are here
    ''
    ''
    '                                                      AEROELASTIC TRIM VARIABLES'
    ''
    '                                  ID     LABEL                 TYPE        TRIM STATUS      VALUE OF UX'
    ''
    '                                         INTERCEPT          RIGID BODY           FIXED      1.000000E+00'
    '                                 101     ANGLEA             RIGID BODY           FIXED      1.000000E-01  RADIANS'
    '1    144                                                                      MARCH  12, 2020  SIMCENTER NASTRAN  3/12/20   PAGE    42'
    """
    while 'AEROELASTIC TRIM VARIABLES' not in line:
        line = f06_file.readline()
        i += 1
    #print(i, line)

    line = f06_file.readline()
    line = f06_file.readline()
    i += 2
    # ID     LABEL                 TYPE        TRIM STATUS      VALUE OF UX

    line = f06_file.readline()
    line = f06_file.readline()
    i += 2

    trim_variables: dict[str, TrimVariable] = {}
    assert 'INTERCEPT' in line, line
    idi, name, trim_type, trim_status, ux, ux_unit = _split_trim_variable(line)
    trim_variables[name] = (idi, trim_type, trim_status, ux, ux_unit)

    line_end, iend, seek1 = _skip_to_page_stamp_and_rewind(f06_file, line, i, nlines_max)

    while i < iend:
        line = f06_file.readline()
        if 'NASTRAN' in line and 'PAGE' in line or 'CONTROL SURFACE POSITION AND HINGE MOMENT RESULTS' in line:
            break
        line2 = line.rstrip()
        if len(line2) == 0:
            continue
        idi, name, trim_type, trim_status, ux, ux_unit = _split_trim_variable(line2)
        trim_variables[name] = (idi, trim_type, trim_status, ux, ux_unit)
    trim_results.trim_variables[isubcase] = TrimVariables.from_f06(
        trim_variables, metadata, isubcase)
    f06_file.seek(seek1)
    return line_end, iend


def _split_trim_variable(line: str) -> tuple[int, str, str, str, float, str]:
    """101     ANGLEA             RIGID BODY           FIXED      1.000000E-01  RADIANS'"""
    line2 = line.rstrip()  # s.split()

    # old
    #id_str = line2[30:40].strip()
    #name = line2[40:50].strip()
    #Type = line2[50:70].strip()
    #trim_status = line2[70:90].strip()
    #ux = line2[90:106]
    #ux_unit = line2[106:130]

    id_str = line2[30:40].strip()
    name = line2[40:50].strip()
    trim_type = line2[50:70].strip()
    trim_status = line2[70:90].strip()
    ux_str = line2[90:106]
    ux_unit = line2[106:130]
    assert line2[130:].strip() == '', line2[130:]
    #idi, name, type, trim_status, ux, ux_unit

    if id_str:
        int_id = int(id_str)
    else:
        int_id = 0

    #print('%r %r %r %r ux=%r %r' % (int_id, name, Type, trim_status, ux_str, ux_unit))
    ux = float(ux_str)
    assert trim_type in {'RIGID BODY', 'CONTROL SURFACE'}, trim_type
    assert trim_status in {'FIXED', 'FREE', 'LINKED'}, trim_status
    assert ux_unit in {'', 'LOAD FACTOR', 'RADIANS', 'NONDIMEN. RATE', 'RAD/S/S PER G'}, ux_unit

    return int_id, name, trim_type, trim_status, ux, ux_unit


def _read_metadata_header(f06_file: TextIO,
                          i: int) -> tuple[int, str, str, str]:
    #'CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC'
    line = f06_file.readline()
    i += 1
    sline = line.strip().split()
    assert sline[0] == 'CONFIGURATION', sline
    aero_config = sline[2]  # 'AEROSG2D'
    assert sline[3] == 'XY-SYMMETRY', sline
    xy_symmetry = sline[5]  # ASYMMETRIC
    assert sline[6] == 'XZ-SYMMETRY', sline
    xz_symmetry = sline[8]  # ASYMMETRIC
    assert xz_symmetry in {'ANTISYMMETRIC', 'ASYMMETRIC', 'SYMMETRIC'}, xz_symmetry
    assert xy_symmetry in {'ANTISYMMETRIC', 'ASYMMETRIC', 'SYMMETRIC'}, xy_symmetry
    return i, aero_config, xz_symmetry, xy_symmetry


def _read_aerostatic_data_recovery_output_table(f06_file: TextIO,
                                                line: str, i: int, nlines_max: int,
                                                trim_results: TrimResults,
                                                title: str, subtitle: str, subcase: str,
                                                ipressure: int, iforce: int,
                                                log: SimpleLogger) -> tuple[str, int, int, int, dict[str, Any]]:
    """
    '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'      <----- you are here
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC'
    '                                           MACH = 0.000000E+00                Q = 1.000000E+00'
    '                         CHORD = 1.0000E+00           SPAN = 1.0000E+01            AREA = 5.0000E+00'
    ''
    ''
    '                                             AERODYNAMIC FORCES ON THE AERODYNAMIC ELEMENTS'
    ''
    '    GROUP  GRID ID  LABEL        T1                T2                T3                R1                R2                R3'
    '        1        1   LS     0.000000E+00      0.000000E+00      6.090764E-03      0.000000E+00      7.613455E-05      0.000000E+00'
    '        1        2   LS     0.000000E+00      0.000000E+00      2.964140E-03      0.000000E+00      3.705175E-05      0.000000E+00'

    """
    isubcase = int(subcase)
    del subcase

    out = _read_metadata_header(f06_file, i)
    (i, aero_config, xz_symmetry, xy_symmetry,
     #mach, q,
     ) = out

    line = f06_file.readline()
    i += 1
    assert 'MACH' in line, line.strip()
    sline = line.strip().split()
    # ['MACH', '=', '7.890000E-01', 'Q', '=', '1.500000E+00']
    assert sline[0] == 'MACH', sline
    assert sline[3] == 'Q', sline
    mach = float(sline[2])
    q = float(sline[5])

    line = f06_file.readline()
    i += 1
    assert 'CHORD' in line, line.strip()
    sline = line.strip().split()
    #['CHORD', '=', '1.3110E+02', 'SPAN', '=', '2.5564E+03', 'AREA', '=', '7.3400E+05']
    assert sline[0] == 'CHORD', sline
    assert sline[3] == 'SPAN', sline
    assert sline[6] == 'AREA', sline
    cref = float(sline[2])
    bref = float(sline[5])
    sref = float(sline[8])
    label = ''
    metadata = {
        'mach': mach, 'q': q,
        'cref': cref, 'bref': bref, 'sref': sref,
        'title': title, 'subtitle': subtitle, 'label': label,
    }
    trim_results.metadata[isubcase] = metadata

    unused_line1 = f06_file.readline()
    line2 = f06_file.readline()
    i += 2
    if 'TRIM ALGORITHM USED: LINEAR TRIM SOLUTION WITHOUT REDUNDANT CONTROL SURFACES.' in line2:
        line, i = _read_aeroelastic_trim_variables(
            f06_file, line, i, nlines_max, trim_results, isubcase, metadata)
        return line, i, ipressure, iforce, metadata

    line3 = f06_file.readline()
    i += 1
    if 'TRANSFORMATION FROM REFERENCE TO WIND AXES:' in line3:
        line, i = _skip_to_page_stamp(f06_file, line, i, nlines_max)
        #line, i, title, subtitle, subcase = _get_title_subtitle_subcase(f06_file, line, i, nlines_max)
        return line, i, ipressure, iforce, metadata

    if 'AERODYNAMIC PRESSURES ON THE AERODYNAMIC ELEMENTS' in line3:
        line, i, grid_id, element_id, Cp_pressure = _read_aerostatic_data_recover_output_table_pressure(
            f06_file, line3, i, nlines_max, log)
        nelement, nresult = Cp_pressure.shape

        if subtitle:
            trim_results.metadata[isubcase]['subtitle'] = subtitle
        if isubcase not in trim_results.aero_pressure:
            trim_results.aero_pressure[isubcase] = []

        trim_results.aero_pressure[isubcase].append((grid_id, element_id, Cp_pressure))
        ipressure += 1

    elif 'AERODYNAMIC FORCES ON THE AERODYNAMIC ELEMENTS' in line3:
        line, i, grid_id, loads, labels = _read_aerostatic_data_recover_output_table_force(
            f06_file, line3, i, nlines_max, log)

        if subtitle:
            trim_results.metadata[isubcase]['subtitle'] = subtitle
        if isubcase not in trim_results.aero_force:
            trim_results.aero_force[isubcase] = []

        trim_results.aero_force[isubcase].append((grid_id, loads, labels))
        iforce += 1
    else:
        raise NotImplementedError(line3.strip())
    return line, i, ipressure, iforce, metadata


def _read_aerostatic_data_recover_output_table_pressure(f06_file: TextIO,
                                                        line: str, i: int, nlines_max: int, log: SimpleLogger) -> tuple[
                                                            str, int, np.ndarray, np.ndarray, np.ndarray]:
    """
    # MSC
    '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC'
    '                                           MACH = 0.000000E+00                Q = 1.000000E+00'
    '                         CHORD = 1.0000E+00           SPAN = 1.0000E+01            AREA = 5.0000E+00'
    ''
    ''
    '                                            AERODYNAMIC PRESSURES ON THE AERODYNAMIC ELEMENTS'  <--- you are here
    ''
    '                                                             AERODYNAMIC PRES.       AERODYNAMIC'
    '                                         GRID   LABEL          COEFFICIENTS           PRESSURES' <---here
    '                                          32     LS            5.390633E-02         5.390633E-02'
    '                                          33     LS            4.857536E-02         4.857536E-02'


    # NX
    '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = ASYMMETRIC'
    '                                           MACH = 8.000000E-01                Q = 1.000000E+00'
    '                         CHORD = 6.9500E+00           SPAN = 7.5000E+01            AREA = 5.2100E+02'
    ''
    ''
    '                                            AERODYNAMIC PRESSURES ON THE AERODYNAMIC ELEMENTS'
    ''
    '                                                             AERODYNAMIC PRES.       AERODYNAMIC'
    '                                         GRID   LABEL          COEFFICIENTS           PRESSURES             ELEMENT'
    '                                        3318     LS            4.405547E-01         4.405547E-01           6020022'
    '                                        3319     LS            5.405547E-01         5.405547E-01           6020023'
    """
    log.debug(' - reading aero pressure')
    unused_line1 = f06_file.readline()
    unused_line2 = f06_file.readline()
    unused_line3 = f06_file.readline()
    line_strip = f06_file.readline().strip()
    i += 4

    data_lines = []
    while line_strip != '':
        data_lines.append(line_strip)
        line = f06_file.readline()
        i += 1
        line_strip = line.strip()
        ndata = len(data_lines)

    ndata = len(data_lines)
    grid_id = np.zeros(ndata, dtype='int32')
    element_id = -np.ones(ndata, dtype='int32')
    cp_pressure = np.zeros((ndata, 2), dtype='float64')
    for i, line in enumerate(data_lines):
        sline = line.split()
        #print(sline)
        if len(sline) == 4:  # MSC
            grid_str, label, cp_str, pressure_str = sline
        else:  # NX
            grid_str, label, cp_str, pressure_str, eid_str = sline
            element_id[i] = int(eid_str)
        grid = int(grid_str)
        cp = float(cp_str)
        pressure = float(pressure_str)
        grid_id[i] = grid
        cp_pressure[i, :] = [cp, pressure]
    return line, i, grid_id, element_id, cp_pressure


def _read_aerostatic_data_recover_output_table_force(f06_file: TextIO,
                                                     line: str, i: int, nlines_max: int,
                                                     log: SimpleLogger):
    """
    '                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S'
    '                         CONFIGURATION = AEROSG2D     XY-SYMMETRY = ASYMMETRIC     XZ-SYMMETRY = SYMMETRIC'
    '                                           MACH = 0.000000E+00                Q = 1.000000E+00'
    '                         CHORD = 1.0000E+00           SPAN = 1.0000E+01            AREA = 5.0000E+00'
    ''
    ''
    '                                             AERODYNAMIC FORCES ON THE AERODYNAMIC ELEMENTS'      <----- you are here
    ''
    '    GROUP  GRID ID  LABEL        T1                T2                T3                R1                R2                R3'
    '        1        1   LS     0.000000E+00      0.000000E+00      6.090764E-03      0.000000E+00      7.613455E-05      0.000000E+00'
    '        1        2   LS     0.000000E+00      0.000000E+00      2.964140E-03      0.000000E+00      3.705175E-05      0.000000E+00'

    """
    log.debug(' - reading aero force')
    assert 'AERODYNAMIC FORCES ON THE AERODYNAMIC ELEMENTS' in line, line.strip()

    line = f06_file.readline()
    line = f06_file.readline()
    i += 2
    assert 'GROUP  GRID ID  LABEL' in line, line.strip()

    line_strip = f06_file.readline().strip()
    i += 1
    data_lines = []
    while line_strip != '':
        data_lines.append(line_strip)

        line = f06_file.readline()
        i += 1
        line_strip = line.strip()

    # parse
    #'    GROUP  GRID ID  LABEL        T1                T2                T3                R1                R2                R3'
    #'        1        1   LS     0.000000E+00      0.000000E+00      6.090764E-03      0.000000E+00      7.613455E-05      0.000000E+00'
    #'        1        2   LS     0.000000E+00      0.000000E+00      2.964140E-03      0.000000E+00      3.705175E-05      0.000000E+00'
    ndata = len(data_lines)
    grid_id = np.zeros(ndata, dtype='int32')
    loads = np.zeros((ndata, 6), dtype='float64')
    labels = []
    for i, line in enumerate(data_lines):
        group_str, grid_str, label, *force_moment = line.split()
        unused_group = int(group_str)
        grid = int(grid_str)
        grid_id[i] = grid
        force_moment_float = [float(val) for val in force_moment]
        loads[i] = force_moment_float
        labels.append(label)
    return line, i, grid_id, loads, labels
