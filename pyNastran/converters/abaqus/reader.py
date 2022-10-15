from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np
from .abaqus_cards import ShellSection, SolidSection, Boundary
if TYPE_CHECKING:
    from cpylog import SimpleLogger


def get_param_map(iline: int, word: str, required_keys: Optional[list[str]]=None) -> dict[str, str]:
    """
    get the optional arguments on a line

    Examples
    --------
    >>> iline = 0
    >>> word = 'elset,instance=dummy2,generate'
    >>> params = get_param_map(iline, word, required_keys=['instance'])
    params = {
        'elset' : None,
        'instance' : 'dummy2,
        'generate' : None,
    }
    """
    if required_keys is None:
        required_keys = []
    words = word.split(',')
    param_map = {}
    for wordi in words:
        if '=' not in wordi:
            key = wordi.strip()
            value = None
        else:
            sword = wordi.split('=', 1)
            assert len(sword) == 2, sword
            key = sword[0].strip()
            value = sword[1].strip()
        param_map[key] = value

    msg = ''
    for key in required_keys:
        if key not in param_map:
            msg += f'line {iline:d}: {key!r} not found in {word!r}\n'
    if msg:
        raise RuntimeError(msg)
    return param_map


def read_cload(line0, lines, iline, log: SimpleLogger) -> tuple[int, str, Any]:
    """
    First line
    ----------
     1. Node number or node set label.
     2. Concentrated load type label, TSB.
     3. Magnitude factor, M. The default value is 1.0. This factor will
        be scaled by any *AMPLITUDE specification associated with this *CLOAD option.
     4. Exposed area.

    Give the following direction cosines in the local coordinate system
    if the *TRANSFORM option was used at this node:
     5. X-direction cosine of the outward normal to the exposed area,
        pointing into the fluid, in the initial configuration.
     6. Y-direction cosine of the outward normal to the exposed area,
        pointing into the fluid, in the initial configuration.
     7. Z-direction cosine of the outward normal to the exposed area,
        pointing into the fluid, in the initial configuration.

    The following data should be provided only if it is necessary to change
    the fluid properties specified under the *AQUA option:
     8. Density of the fluid outside the element. This value will override
        the fluid density given on the data line of the *AQUA option.
     9. Free surface elevation of the fluid outside the element. This value
        will override the fluid surface elevation given on the data line
        of the *AQUA option.
     10. Constant pressure, added to the hydrostatic pressure outside the element.
    Repeat this data line as often as necessary to define concentrated
    buoyancy at various nodes or node sets.
    """
    log.debug(f'read_cload {line0!r}')
    cload = []
    while '*' not in line0:
        sline = line0.split(',')
        assert len(sline) == 3, sline
        #cload += sline
        iline += 1
        line0 = lines[iline].strip().lower()
        #print(line0)
        try:
            nid = int(sline[0])
        except ValueError:
            nid = sline[0]
        dof = int(sline[1])
        mag = float(sline[2])
        cloadi = (nid, dof, mag)
        cload.append(cloadi)

    return iline, line0, cload

def read_node(lines, iline, log, skip_star=False):
    """reads *node"""
    if skip_star:
        iline += 1

    nids = []
    nodes = []
    #print('  Node iline=%s' % iline)
    line0 = lines[iline].strip().lower()
    assert '*' not in line0, line0
    #print('  node line0 =', line0)
    is_failed = False
    #if len(nids) > 0:
        #nids0 = copy.deepcopy(nids)
        #nids = []
        #is_failed = False

    while not line0.startswith('*'):
        sline = line0.split(',')
        nids.append(sline[0])
        nsline = len(sline)
        if nsline == 3:
            sline.append(0.)
            nodes.append(sline[1:])
        elif nsline == 4:
            nodes.append(sline[1:])
        else:
            raise NotImplementedError(sline)
        iline += 1
        line0 = lines[iline].strip().lower()
    unused_nnodes = len(nids)

    if is_failed:
        msg = 'nids will overwrite nids0!\n'
        #msg += 'nids0 = %s\n' % nids0
        msg += 'nids = %s\n' % nids
        raise RuntimeError(msg)
    return iline, line0, nids, nodes

def read_elset(lines: list[str], iline: int, word: str, log: SimpleLogger,
               is_instance: bool=True):
    """reads *elset"""
    log.debug('word=%r' % word)
    assert '*' in word, f'word={word!r} param_map={param_map}'
    params_map = get_param_map(iline, word, required_keys=['elset'])
    # TODO: skips header parsing
    iline += 1
    #print('elset: ', line0)
    params_map = get_param_map(iline, word, required_keys=['elset'])
    set_name = params_map['elset']
    line0 = lines[iline].strip().lower()
    set_ids, iline, line0 = read_set(lines, iline, line0, params_map)
    return iline, line0, set_name, set_ids

def read_nset(lines, iline, word, log, is_instance=True):
    """reads *nset"""
    log.debug('word=%r' % word)
    assert 'nset' in word, word
    # TODO: skips header parsing
    required_keys = ['instance'] if is_instance else []
    params_map = get_param_map(iline, word, required_keys=required_keys)
    #print('params_map =', params_map)
    set_name = params_map['nset']
    iline += 1
    line0 = lines[iline].strip().lower()
    set_ids, iline, line0 = read_set(lines, iline, line0, params_map)
    return iline, line0, set_name, set_ids

def read_boundary(lines: list[str], line0: str, iline: int) -> tuple[Boundary, int, str]:
    boundary_lines = []
    line0 = lines[iline]
    assert '*' not in line0, line0
    while '*' not in line0:
        # nid, dof1, dof2, disp
        #1,1,,0
        sline = line0.strip().split(',')
        boundary_lines.append(sline)
        iline += 1
        line0 = lines[iline].strip().lower()
    boundary = Boundary.from_lines(boundary_lines)
    return boundary, iline, line0

def read_solid_section(line0, lines, iline, log):
    """reads *solid section"""
    # TODO: skips header parsing
    #iline += 1
    assert '*solid' in line0, line0
    word2 = line0.strip('*').lower()
    params_map = get_param_map(iline, word2, required_keys=['material'])
    log.debug('    param_map = %s' % params_map)
    #line0 = lines[iline].strip().lower()
    data_lines, iline, line0 = read_star_block2(lines, iline, line0, log)
    #print('line0 =', iline, line0)
    #print(f'lines[{iline}] = {lines[iline]!r}')
    #print('lines[iline+1] =', lines[iline+1])
    #print('data_lines =', data_lines)
    #for line in data_lines:
        #print(line)
    solid_section = SolidSection.add_from_data_lines(params_map, data_lines, log)
    return iline, solid_section

def read_shell_section(line0: str, lines: list[str], iline: int,
                       log: SimpleLogger) -> ShellSection:
    """reads *shell section"""
    assert '*shell' in line0, line0
    # TODO: skips header parsing
    #iline += 1
    word2 = line0.strip('*').lower()
    params_map = get_param_map(iline, word2, required_keys=['material'])
    log.debug('    param_map = %s' % params_map)

    iline += 1
    line0 = lines[iline].strip().lower()
    data_lines, iline, line0 = read_star_block2(lines, iline, line0, log)
    log.info(f'params_map = {params_map}')
    log.info(f'data_lines = {data_lines}')
    #for line in data_lines:
        #print(line)
    assert len(data_lines) > 0, data_lines
    shell_section = ShellSection.add_from_data_lines(params_map, data_lines, log)
    #print(lines[iline])
    return iline, shell_section

def read_hourglass_stiffness(line0: str, lines: list[str], iline: int,
                             log: SimpleLogger) -> None:
    """reads *hourglass stiffness"""
    # TODO: skips header parsing
    #iline += 1
    word2 = line0.strip('*').lower()
    iline += 1
    #params_map = get_param_map(iline, word2, required_keys=['material'])
    #log.debug('    param_map = %s' % params_map)
    #line0 = lines[iline].strip().lower()
    data_lines, iline, line0 = read_star_block2(lines, iline, line0, log)
    assert len(data_lines) == 1, data_lines
    #for line in data_lines:
        #print(line)
    #solid_section = SolidSection(params_map, data_lines, log)
    hourglass_stiffness = None
    return iline, hourglass_stiffness


def read_star_block(lines, iline, line0, log, debug=False):
    """
    because this uses file streaming, there are 30,000 places where a try except
    block is needed, so this should probably be used all over.
    """
    data_lines = []
    try:
        iline += 1
        line0 = lines[iline].strip().lower()
        while not line0.startswith('*'):
            data_lines.append(line0.split(','))
            iline += 1
            line0 = lines[iline].strip().lower()
            #log.debug('line = %r' % line0)
        iline -= 1
        line0 = lines[iline].strip().lower()
    except IndexError:
        pass
    if debug:
        for line in data_lines:
            log.debug(line)
    return data_lines, iline, line0


def read_star_block2(lines, iline, line0, log, debug=False):
    """
    because this uses file streaming, there are 30,000 places where a try except
    block is needed, so this should probably be used all over.
    """
    line0 = lines[iline].strip().lower()
    data_lines = []
    while not line0.startswith('*'):
        sline = line0.strip(', ').split(',')
        data_lines.append(sline)
        iline += 1
        line0 = lines[iline].strip().lower()
    if debug:
        for line in data_lines:
            log.debug(line)
    return data_lines, iline, line0

def read_set(lines, iline, line0, params_map):
    """reads a set"""
    set_ids = []
    while not line0.startswith('*'):
        set_ids += line0.strip(', ').split(',')
        iline += 1
        line0 = lines[iline].strip().lower()
    if 'generate' in params_map:
        assert len(set_ids) == 3, set_ids
        set_ids = np.arange(int(set_ids[0]), int(set_ids[1]), int(set_ids[2]))
    else:
        try:
            set_ids = np.unique(np.array(set_ids, dtype='int32'))
        except ValueError:
            print(set_ids)
            raise
    return set_ids, iline, line0

def read_orientation(line0: str, lines: list[str], iline: int, log: SimpleLogger):
    orientation_fields = []
    iline += 1
    line0 = lines[iline].strip().lower()
    while '*' not in line0:
        sline = line0.split(',')
        iline += 1
        line0 = lines[iline].strip().lower()
        orientation_fields.append(sline)
    log.debug(line0)
    iline -= 1
    line0 = lines[iline].strip().lower()
    return iline, line0, orientation_fields

def read_system(line0: str, lines: list[str], iline: int, log: SimpleLogger):
    coordinate_system_fields = []
    iline += 1
    line0 = lines[iline].strip().lower()
    while '*' not in line0:
        sline = line0.split(',')
        iline += 1
        line0 = lines[iline].strip().lower()
        coordinate_system_fields.extend(sline)
    log.debug(line0)
    iline -= 1
    line0 = lines[iline].strip().lower()
    assert len(coordinate_system_fields) == 9, coordinate_system_fields
    return iline, line0, coordinate_system_fields
