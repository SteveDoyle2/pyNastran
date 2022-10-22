from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from .abaqus_cards import (
    ShellSection, SolidSection, Boundary, Material, Transform)
from .elements import allowed_element_types
from .reader_utils import split_by_equals

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

def read_element(lines: list[str], line0: str, iline: int, log: SimpleLogger, debug: bool):
    """
    '*element, type=mass, elset=topc_inertia-2_mass_'
    """
    #print('------------------')
    assert '*' in line0, line0
    sline = line0.split(',')[1:]
    if len(sline) < 1:
        raise RuntimeError("looking for element_type (e.g., '*Element, type=R2D2')\n"
                           "line0=%r\nsline=%s; allowed:\n[%s]" % (
                               line0, sline, ', '.join(allowed_element_types)))

    etype_sline = sline[0]
    assert 'type' in etype_sline, etype_sline
    etype = etype_sline.split('=')[1].strip()
    if etype not in allowed_element_types:
        msg = 'etype=%r allowed=[%s]' % (etype, ','.join(allowed_element_types))
        raise RuntimeError(msg)

    if debug:
        log.debug('    etype = %r' % etype)

    #iline += 1
    line1 = lines[iline].strip().lower()
    log.debug('    line1 = %r' % line1)

    elements = []
    #print(line1)
    assert '*' not in line1, line1
    while not line1.startswith('*'):
        #print(line1)
        elements.append(line1.strip('\n\t ,').split(','))
        iline += 1
        line1 = lines[iline].strip().lower()
    #log.debug('elements = %s' % elements)
    return line1, iline, etype, elements

def read_spring(lines: list[str], iline: int, word: str, log: SimpleLogger) -> None:
    """
    *SPRING,ELSET=Eall
    blank line
    10.

    Defines a linear spring constant with value 10. for all elements in
    element set Eall and all temperatures.
    Example:
    *SPRING,ELSET=Eall,NONLINEAR
    0.,0.,293.
    10.,1.,293.
    100.,2.,293.
    0.,0.,393.
    5.,1.,393.
    25.,2.,393.
    """
    log.warning('spring is unsupported')
    param_map = get_param_map(iline, word, required_keys=['elset'])
    print(param_map)
    #name = param_map['name']

    iline += 1
    word_line = lines[iline].strip().lower()
    word = word_line.strip('*').lower()
    unused_allowed_words = ['elastic']
    unallowed_words = [
        'shell section', 'solid section',
        'material', 'step', 'boundary', 'amplitude', 'surface interaction',
        'assembly', 'spring']
    iline += 1
    line0 = lines[iline].strip('\n\r\t, ').lower()

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

def read_material(lines: list[str], iline: int, word: str, log: SimpleLogger) -> Material:
    """reads a Material card"""
    param_map = get_param_map(iline, word, required_keys=['name'])
    #print(param_map)
    name = param_map['name']

    iline += 1
    word_line = lines[iline].strip().lower()
    word = word_line.strip('*').lower()
    unused_allowed_words = ['elastic']
    unallowed_words = [
        'shell section', 'solid section',
        'material', 'step', 'boundary', 'amplitude', 'surface interaction',
        'assembly', 'spring']
    iline += 1
    line0 = lines[iline].strip('\n\r\t, ').lower()
    #print('  wordA =', word)
    #while word in allowed_words:
    sections = {}
    density = None
    ndelete = None
    ndepvars = None
    while word not in unallowed_words:
        data_lines = []
        #log.info('  mat_word = %r' % word)
        #print(sections)
        if word.startswith('elastic'):
            key = 'elastic'
            sword = word.split(',')
            #print(key, sword)

            #log.debug('  matword = %s' % sword)
            if len(sword) == 1:
                # elastic
                assert len(sword) in [1, 2], sword
            else:
                mat_type = sword[1]
                assert 'type' in mat_type, sword
                mat_type = mat_type.split('=')[1]

                sline = line0.split(',')
                if mat_type == 'traction':
                    assert len(sline) == 3, sline
                    log.debug('  traction material')
                elif mat_type in ['iso', 'isotropic']:
                    #, TYPE=ISO
                    #1.00000E+07, 3.00000E-01
                    assert len(sline) == 2, sline
                    e, nu = sline
                    e = float(e)
                    nu = float(nu)
                    sections['elastic'] = [e, nu]
                    #print(sections)
                else:
                    raise NotImplementedError(f'mat_type={mat_type!r}')
            iline += 1
        elif word.startswith('plastic'):
            key = 'plastic'
            sword = word.split(',')
            log.debug('  matword = %s' % sword)
            if len(sword) == 1:
                # elastic
                assert len(sline) in [1, 2], sline
            else:
                raise NotImplementedError(sline)
            data_lines, iline, line0 = read_star_block2(lines, iline, line0, log, debug=False)
            #print(data_lines)
        elif word == 'density':
            key = 'density'
            sline = line0.split(',')
            assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
            density = float(sline[0])
            iline += 1
        elif word.startswith('damage initiation'):
            key = 'damage initiation'
            #log.debug('  damage0 %s' % line0)
            sline = line0.split(',')
            log.debug(sline)
            assert len(sline) == 3, sline
            iline += 1
        elif word.startswith('damage evolution'):
            key = 'damage evolution'
            #self.log.debug('  damage_e %s' % line0)
            unused_data = []
            while '*' not in line0:
                sline = line0.split(',')
                assert len(sline) == 3, sline
                iline += 1
                line0 = lines[iline].strip().lower()
            log.debug(line0)
        elif word == 'damage stabilization':
            key = 'damage stabilization'
            sline = line0.split(',')
            assert len(sline) == 1, sline
            iline += 1

        #elif word.startswith('surface interaction'):
            #key = 'surface interaction'
            #data = []
            #while '*' not in line0:
                #sline = line0.split(',')
                #iline += 1
                #line0 = lines[iline].strip().lower()
            #log.debug(line0)
        #elif word.startswith('friction'):
            #key = 'friction'
            #data = []
            #while '*' not in line0:
                #sline = line0.split(',')
                #iline += 1
                #line0 = lines[iline].strip().lower()
            #log.debug(line0)
        #elif word.startswith('surface behavior'):
            #key = 'surface behavior'
            #data = []
            #while '*' not in line0:
                #sline = line0.split(',')
                #iline += 1
                #line0 = lines[iline].strip().lower()
            #log.debug(line0)
        #elif word.startswith('contact damping'):
            #key = 'contact damping'
            #data = []
            #while '*' not in line0:
                #sline = line0.split(',')
                #iline += 1
                #line0 = lines[iline].strip().lower()
            #log.debug(line0)

        elif word.startswith('depvar'):
            key = 'depvar'
            sline = word_line.split()
            if len(sline) > 1:
                assert len(sline) == 2, sline
                sline2 = sline[1].split('=')
                assert  len(sline2) == 2, sline
                assert sline2[0].lower() == 'delete', sline
                ndelete = int(sline2[1])

            sline = line0.split(',')
            assert len(sline) == 1, sline
            ndepvars = int(sline[0])
            iline += 1
        elif word.startswith('user material'):
            key = 'user material'
            words = word.split(',')[1:]

            is_constants = False
            for wordi in words:
                mat_word, value = split_by_equals(wordi, lines, iline-1)
                mat_word = mat_word.strip()
                if mat_word == 'constants':
                    nconstants = int(value)
                    is_constants = True
                elif mat_word == 'type':
                    mat_type = value.strip()
                    allowed_types = ['mechanical']
                    if not mat_type in allowed_types:
                        msg = 'mat_type=%r; allowed_types=[%s]'  % (
                            mat_type, ', '.join(allowed_types))
                        raise NotImplementedError(msg)
                else:
                    raise NotImplementedError('mat_word=%r' % mat_word)

            if not is_constants:
                msg = "line %i: 'constants' was not defined on %r" % (
                    iline, lines[iline-1].rstrip())
                raise RuntimeError(msg)

            #nconstants = 111
            nlines_full = nconstants // 8
            nleftover = nconstants % 8
            mat_data = []
            for unused_iiline in range(nlines_full):
                sline = line0.split(',')
                assert len(sline) == 8, 'len(sline)=%s; sline=%s' % (len(sline), sline)
                mat_data += sline
                iline += 1
                line0 = lines[iline].strip('\n\r\t, ').lower()
            if nleftover:
                sline = line0.split(',')
                iline += 1
                line0 = lines[iline].strip('\n\r\t, ').lower()
        elif word.startswith('initial conditions'):
            # TODO: skips header parsing
            #iline += 1
            #line0 = lines[iline].strip().lower()
            unused_data = []
            while '*' not in line0:
                sline = line0.split(',')
                iline += 1
                line0 = lines[iline].strip().lower()
            log.debug(line0)
        elif word.lower().startswith('hyperelastic, mooney-rivlin'):
            key = 'hyperelastic, mooney-rivlin'
            while '*' not in line0:
                sline = line0.split(',')
                iline += 1
                line0 = lines[iline].strip().lower()
            log.debug(line0)
        elif word.lower().startswith('expansion'):
            #*Expansion, zero=20.
            #80.,
            key = 'expansion'
            while '*' not in line0:
                sline = line0.split(',')
                iline += 1
                line0 = lines[iline].strip().lower()
            #iline += 1
            log.debug(line0)
        else:
            msg = print_data(lines, iline, word, 'is this an unallowed word for *Material?\n')
            raise NotImplementedError(msg)

        if key in sections:
            msg = f'key={key!r} already defined for Material name={name!r}'
            log.warning(msg)
        else:
            #raise RuntimeError(msg)
            sections[key] = data_lines

        try:
            line = lines[iline]
        except IndexError:
            is_broken = True
            log.debug('  breaking on end of file')
            break
        word_line = line.strip('\n\r\t, ').lower()
        del line
        word = word_line.strip('*').lower()

        iline += 1
        line0 = lines[iline].strip('\n\r\t, ').lower()
        #log.debug('  lineB = %r' % line0)
        #log.debug('  wordB = %r' % word)

        is_broken = False
        for unallowed_word in unallowed_words:
            if word.startswith(unallowed_word):
                log.debug('  breaking on %r' % unallowed_word)
                is_broken = True
                break
        if is_broken:
            iline -= 1
            break
    #print(name, sections)
    material = Material(name, sections=sections,
                        is_elastic=True, density=density,
                        ndepvars=ndepvars, ndelete=ndelete)
    iline -= 1
    return iline, line0, material

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
    boundary = Boundary.from_data(boundary_lines)
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
    """
    *SYSTEM
    6414.0    , 0.0       , -678.0    , 6414.03642, 0.81906834, -677.42746
    6414.05199, -0.573696 , -677.18258
    *NODE, SYSTEM=C
        241,  1724.0602494441,  54.991639723866,  -0.280338873236
    """
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

def read_transform(line0: str, lines: list[str], iline: int, log: SimpleLogger):
    """
    *TRANSFORM, TYPE=C, NSET=HM_auto_transform_3
    6414.0    , 0.0       , -678.0    ,  6513.79832,  -4.661E-06,  -684.34787
    *NSET, NSET=HM_auto_transform_3
             3,       241

    First (and only) line:
      Global X-coordinate of point a
      Global Y-coordinate of point a
      Global Z-coordinate of point a
      Global X-coordinate of point b
      Global Y-coordinate of point b
      Global Z-coordinate of point b
    """
    params = get_param_map(iline, line0, required_keys=['type', 'nset'])
    transform_type = params['type'].upper()
    nset = params['nset']
    assert transform_type in ['R', 'C', 'S'], f'transform_type={transform_type!r}'

    transform_fields = []
    iline += 1
    line0 = lines[iline].strip().lower()
    while '*' not in line0:
        sline = line0.split(',')
        iline += 1
        line0 = lines[iline].strip().lower()
        transform_fields.extend(sline)
    log.debug(line0)
    iline -= 1
    line0 = lines[iline].strip().lower()
    assert len(transform_fields) == 6, transform_fields
    transform = Transform.from_data(transform_type, nset, transform_fields)
    return iline, line0, transform
