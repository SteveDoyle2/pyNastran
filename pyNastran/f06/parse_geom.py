from typing import Optional
from pyNastran.utils import PathLike
from cpylog import SimpleLogger, get_logger


def parse_f06_geom(f06_filename: PathLike,
                   log: Optional[SimpleLogger]=None) -> tuple[list[str], list[str],
                                                              list[str], list[str]]:
    """C A S E    C O N T R O L   D E C K   E C H O"""
    log: SimpleLogger = get_logger(log)
    iline = 0
    nblank = 0
    is_section = False
    system_lines = []
    exec_lines = []
    case_lines = []
    bulk_lines = []
    with open(f06_filename, 'r') as f06_file:
        while nblank < 50:
            if not is_section:
                line = f06_file.readline().rstrip(); iline += 1
            if line == '':
                nblank += 1
                continue
            nblank = 0

            # log.debug(f'A {iline}: {line!r}')
            if 'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O' in line:
                # log.info(f'B found FILE AND SYSTEM PARAMETER: {line!r}')
                iline, line = _get_geom_section(f06_file, iline, system_lines, log)
                is_section = True
            elif 'N A S T R A N    E X E C U T I V E    C O N T R O L    D E C K    E C H O' in line:
                # log.info(f'B found EXECUTIVE CONTROL DECK: {line!r}')
                iline, line = _get_geom_section(f06_file, iline, exec_lines, log)
                is_section = True
            elif 'C A S E    C O N T R O L   D E C K   E C H O' in line:
                # log.info(f'B found CASE CONTROL DECK: {line!r}')
                iline, line = _get_geom_section(f06_file, iline, case_lines, log)
                # log.info(f'end CASE CONTROL DECK {iline}: {line!r}')
                is_section = True
            elif 'INPUT BULK DATA CARD COUNT' in line:
                # INPUT BULK DATA CARD COUNT
                #  - sometimes the f06 doesn't write the 'B U L K   D A T A   E C H O' header
                #  - sometimes it writes both and this section is empty...
                iline, line = _get_geom_section(f06_file, iline, bulk_lines, log, require_lines=False)
                is_section = True
            elif 'S O R T E D   B U L K   D A T A   E C H O' in line or 'INPUT BULK DATA CARD COUNT' in line:
                # log.info(f'B found BULK DATA DECK: {line!r}')
                iline, line = _get_geom_section(f06_file, iline, bulk_lines, log)
                is_section = True

            elif ('O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R' in line or
                  'R E A L   E I G E N V A L U E S' in line or
                  'R E A L   E I G E N V E C T O R   N O .' in line):
                break
            else:
                is_section = False
                #print('???')
                continue

                # while nblank < 50:
                #     line = f06_file.readline().rstrip(); iline += 1
                #     print('B', iline, line)
                #     if line == '':
                #         nblank += 1
                #         continue
                #     nblank = 0

    system_lines = [line.lstrip() for line in system_lines]
    exec_lines = [line.lstrip() for line in exec_lines]
    case_lines = clean_case_lines(case_lines)
    bulk_lines = clean_bulk_lines(bulk_lines)

    # for i, linei in enumerate(system_lines):
    #     print(f'system {i}: {linei}')
    # for i, linei in enumerate(exec_lines):
    #     print(f'exec {i}: {linei}')
    # for i, linei in enumerate(case_lines):
    #     print(f'case {i}: {linei}')
    # for i, linei in enumerate(bulk_lines):
    #     print(f'bulk {i}: {linei}')
    #asdf
    return system_lines, exec_lines, case_lines, bulk_lines

def clean_case_lines(lines: list[str]) -> list[str]:
    case_lines1 = [line.lstrip() for line in lines]
    case_lines = [line.split(' ', 1)[1].lstrip() for line in case_lines1 if ' ' in line]
    return case_lines


def clean_bulk_lines(lines: list[str]) -> list[str]:
    try:
        out_lines = [line.split('-        ', 1)[1].lstrip() for line in lines if '-' in line]
    except:
        out_lines = []
        for line in lines:
            if '-        ' not in line:
                continue
            # print(f'{line!r}')
            sline = line.split('- ', 1)
            line2 = sline[1].lstrip()
            out_lines.append(line2)
    return out_lines

def _get_geom_section(f06_file, iline: int, lines: list[str],
                      log: SimpleLogger, require_lines: bool=True) -> tuple[int, str]:
    debug = True
    line = ''
    nblank = 0
    #lines = []
    while not line.startswith('0'):
        line = f06_file.readline().rstrip(); iline += 1
        # log.debug(f'  section: {line!r}')
        if len(line) > 1 and line.startswith('0'):
            #print(f'*** {line}')
            break
        if line == '':
            nblank += 1
            if nblank == 50:
                raise RuntimeError('nblank={nblank:d}; stopping')
            continue
        nblank = 0
        if line == '0':
            line = ' 0'
            continue

        #print(f'aaa: {line.lstrip()!r}')
        strip_line = line.lstrip()
        if strip_line == 'CARD' or strip_line.startswith(('COUNT ', 'INPUT BULK DATA CARD COUNT')):
            continue
        lines.append(line)
    #print(lines)
    if require_lines:
        assert len(lines), lines
    return iline, line
