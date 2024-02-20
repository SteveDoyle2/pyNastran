from __future__ import annotations
from typing import Union, Optional, Any, TYPE_CHECKING
import numpy as np

from .abaqus_cards import (
    Mass, BeamSection, ShellSection, SolidSection,
    Boundary, Material, Transform,
    Surface, Tie, Orientation,
    Frequency, Step, Part,
)
from .elements import allowed_element_types
from .reader_utils import split_by_equals, print_data

if TYPE_CHECKING:
    from cpylog import SimpleLogger


def get_param_map(iline: int, word: str,
                  required_keys: Optional[list[str]]=None) -> dict[str, str]:
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


def read_cload(iline: int, line0: str,
               lines: list[str],
               log: SimpleLogger) -> tuple[int, str, Any]:
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
    cload: list[tuple[Union[int, str], int, float]] = []
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

def read_dload(iline: int, line0: str,
               lines: list[str],
               log: SimpleLogger) -> tuple[int, str, Any]:
    """
    First line:

    *DLOAD
    Enter any needed parameters and their value
    Following line for surface loading:

    Element number or element set label.
    Distributed load type label.
    Actual magnitude of the load (for Px type labels) or fluid node number (for PxNU type labels)
    Repeat this line if needed.
    Example:

    *DLOAD,AMPLITUDE=A1
    element_set_name,P3,10.
    assigns a pressure loading with magnitude 10. times the amplitude curve of
    amplitude A1 to face number three of all elements belonging to set element_set_name.

    for hexahedral elements:
    face 1: 1-2-3-4
    face 2: 5-8-7-6
    face 3: 1-5-6-2
    face 4: 2-6-7-3
    face 5: 3-7-8-4
    face 6: 4-8-5-1

    for tetrahedral elements:
    Face 1: 1-2-3
    Face 2: 1-4-2
    Face 3: 2-4-3
    Face 4: 3-4-1
    for wedge elements:

    Face 1: 1-2-3
    Face 2: 4-5-6
    Face 3: 1-2-5-4
    Face 4: 2-3-6-5
    Face 5: 3-1-4-6
    for quadrilateral plane stress, plane strain and axisymmetric elements:
    Face 1: 1-2
    Face 2: 2-3
    Face 3: 3-4
    Face 4: 4-1
    for triangular plane stress, plane strain and axisymmetric elements:

    Face 1: 1-2
    Face 2: 2-3
    Face 3: 3-1
    for beam elements:

    Face 1: pressure in 1-direction
    Face 2: pressure in 2-direction

    Pressure
    --------
    ['Internal-1_Internal_Selection-1_Uniform_Pressure-1_S1', 'P1', '-1']

    Gravity
    -------
    ['Internal_Selection-1_Gravity-1', ' Grav', ' 9810', ' 0', ' 0', ' -1']

    Following line for gravity loading with known gravity vector:
    Element number or element set label.
    GRAV
    Actual magnitude of the gravity vector.
    Coordinate 1 of the normalized gravity vector
    Coordinate 2 of the normalized gravity vector
    Coordinate 3 of the normalized gravity vector
    Repeat this line if needed. Here "gravity" really stands for any acceleration vector.
    Example:

    *DLOAD
    Eall,GRAV,9810.,0.,0.,-1.
    assigns gravity loading in the negative z-direction with magnitude 9810. to all elements.
    """
    log.debug(f'read_dload {line0!r}')
    line0 = lines[iline]
    if len(line0) and line0[0] == '*':
        dload = None
        return iline, line0, dload

    dload = []
    while '*' not in line0:
        sline = line0.split(',')
        tag = sline[1].upper().strip()
        if tag in {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'}:  # 6 faces on a CHEXA
            assert len(sline) == 3, sline
            element_set_name, face, mag_str = sline
            try:
                element_set_name = int(element_set_name)
            except ValueError:
                element_set_name = sline[0]
            face = face.strip().upper()
            assert face in {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'}, face
            mag = float(mag_str)
            dloadi = (element_set_name, face, mag)
        elif tag == 'GRAV':
            element_set_name, unused_tag, mag_str, gx_str, gy_str, gz_str = sline
            try:
                element_set_name = int(element_set_name)
            except ValueError:
                element_set_name = sline[0]
            mag = float(mag_str)
            gx = float(gx_str)
            gy = float(gy_str)
            gz = float(gz_str)
            dloadi = (element_set_name, tag, mag, gx, gy, gz)
        else:
            raise NotImplementedError(tag)

        #cload += sline
        iline += 1
        line0 = lines[iline].strip().lower()
        #print(line0)
        dload.append(dloadi)
    return iline, line0, dload

def read_surface(iline: int, line0: str, lines: list[str],
                 log: SimpleLogger) -> tuple[int, str, Surface]:
    """
    *Surface, Name=Internal_Selection-1_Uniform_Pressure-1, Type=Element
    Internal-1_Internal_Selection-1_Uniform_Pressure-1_S1, S1

    *Surface, Name=Internal_Selection-1_Uniform_Pressure-1, Type=Element
    Internal-1_Internal_Selection-1_Uniform_Pressure-1_S1, S1
    Internal-1_Internal_Selection-1_Uniform_Pressure-1_S2, S2

    First line:

    *SURFACE
    Enter the parameter NAME and its value, and, if necessary, the TYPE parameter.
    Following line for nodal surfaces:

    Node or node set to be assigned to this surface (maximum 1 entry per line).
    Repeat this line if needed.
    Following line for element face surfaces:

    Element or element set (maximum 1 entry per line).
    Surface label (maximum 1 entry per line).
    https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node247.html
    """
    log.debug(f'read_surface {line0!r}')
    iline += 1
    iline, line, flags, lines_out = read_generic_section(iline, line0, lines, log)

    surface_name = ''
    surface_type = ''
    for key, value in split_strict_flags(flags):
        if key == 'name':
            surface_name = value
        elif key == 'type':
            surface_type = value
        else:  # pragma: no cover
            raise NotImplementedError((key, value))
    assert surface_type in {'element'}, surface_type

    set_names = []
    faces = []
    for line in lines_out:
        sline = line.strip().split(',')
        assert len(sline) == 2, sline
        set_name, face = sline
        set_names.append(set_name.strip().lower())
        faces.append(face.strip().lower())
    surface = Surface(surface_name, surface_type, set_names, faces)
    return iline, line0, surface


def read_frequency(iline: int, line0: str, lines: list[str],
                   log: SimpleLogger) -> tuple[int, str, Frequency]:
    """
    *Frequency, solver=Pardiso
    10
    """
    log.debug(f'read_frequency {line0!r}')
    iline += 1
    iline, line, flags, lines_out = read_generic_section(iline, line0, lines, log)

    solver = ''
    for key, value in split_strict_flags(flags):
        if key == 'solver':
            solver = value
        else:  # pragma: no cover
            raise NotImplementedError((key, value))

    assert len(lines_out) == 1, lines_out
    line = lines_out[0]
    sline = line.strip().split(',')
    assert len(sline) == 1, sline
    nmodes = int(sline[0])
    frequency = Frequency(solver, nmodes)
    return iline, line0, frequency

def split_strict_flags(flags: list[str]) -> list[tuple[str, str]]:
    for key_value in flags:
        key, value = key_value.split('=')
        key = key.strip().lower()
        value = value.strip().lower()
        yield (key, value)

def read_node(iline: int, lines: list[str], log: SimpleLogger,
              skip_star: bool=False) -> tuple[int, str,
                                             list[str],
                                             list[str]]:
    """reads *node"""
    if skip_star:
        iline += 1

    nids: list[str] = []
    nodes: list[str] = []
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
            sline.append('0.')
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

def read_element(iline: int, line0: str, lines: list[str],
                 log: SimpleLogger, debug: bool) -> tuple[int, str,
                                                          str, str, list[tuple[str, ...]]]:
    """
    '*element, type=mass, elset=topc_inertia-2_mass_'
    """
    assert isinstance(iline, int), iline
    assert isinstance(line0, str), line0
    assert '*' in line0, line0
    iline, line_out, flags, lines_out = read_generic_section(iline, line0, lines, log)
    assert len(lines_out), lines_out

    if len(flags) < 1:
        raise RuntimeError("looking for element_type (e.g., '*Element, type=R2D2')\n"
                           "line0=%r\nsline=%s; allowed:\n[%s]" % (
                               line0, flags, ', '.join(allowed_element_types)))

    etype = ''
    elset = ''
    for flag_value in flags:
        flag, value = flag_value.split('=')
        flag = flag.strip().lower()
        value = value.strip()
        if flag == 'type':
            assert etype == '', etype
            etype = value
        elif flag == 'elset':
            elset = value
        else:
            raise RuntimeError(flag_value)

    if etype not in allowed_element_types:
        msg = 'etype=%r allowed=[%s]' % (etype, ','.join(allowed_element_types))
        raise RuntimeError(msg)

    elements = []
    for line in lines_out:
        element = line.strip('\n\t ,').split(',')
        elements.append(element)

    return iline, line_out, etype, elset, elements

def read_spring(iline: int, lines: list[str], word: str, log: SimpleLogger) -> None:
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
    #unallowed_words = [
        #'shell section', 'solid section',
        #'material', 'step', 'boundary', 'amplitude', 'surface interaction',
        #'assembly', 'spring']
    iline += 1
    line0 = lines[iline].strip('\n\r\t, ').lower()

def read_nset(iline: int, line0: str, lines: list[str],
              log: SimpleLogger,
              is_instance: bool) -> tuple[int, str, str, np.ndarray]:
    #line0_backup = line0
    iline, line0, flags, lines_out = read_generic_section(iline, line0, lines, log)

    generate = False
    nset = ''
    instance_name = ''
    for key_value in flags:
        if key_value == 'generate':
            generate = True
            continue
        key, value = key_value.split('=', 1)
        key = key.strip().lower()
        value = value.strip().lower()

        instance_name = ''
        if key == 'instance':
            instance_name = value
        elif key == 'nset':
            nset = value
        else:  # pragma: no cover
            raise RuntimeError((key, value))

    if is_instance:
        set_name = instance_name
    else:
        set_name = nset
    #print(flags)

    set_ids = _generate_set(lines_out, generate)
    assert isinstance(iline, int), iline
    assert len(set_name) > 0, flags
    return iline, line0, set_name, set_ids

def read_elset(iline: int, line0: str, lines: list[str],
               log: SimpleLogger,
               is_instance: bool) -> tuple[int, str, str, np.ndarray]:
    #line0_backup = line0
    iline, line0, flags, lines_out = read_generic_section(iline, line0, lines, log)

    generate = False
    instance_name = ''
    elset = ''
    for key_value in flags:
        if key_value == 'generate':
            generate = True
            continue
        key, value = key_value.split('=', 1)
        key = key.strip().lower()
        value = value.strip().lower()

        if key == 'instance':
            assert elset == '', elset
            instance_name = value
        elif key == 'elset':
            assert instance_name == '', instance_name
            elset = value
        else:  # pramga: no cover
            raise RuntimeError((key, value))

    if is_instance:
        assert elset == '', elset
        set_name = instance_name
    else:
        assert instance_name == '', instance_name
        set_name = elset

    set_ids = _generate_set(lines_out, generate)
    assert isinstance(iline, int), iline
    assert len(set_name) > 0, flags
    assert set_name == set_name.lower()
    return iline, line0, set_name, set_ids

def read_material(iline: int, word: str,
                  lines: list[str],
                  log: SimpleLogger) -> tuple[int, str, str, Material]:
    """reads a Material card"""
    lines2 = lines[iline:]
    assert isinstance(iline, int), iline
    assert isinstance(word, str), word

    param_map = get_param_map(iline, word, required_keys=['name'])
    #print(param_map)
    name = param_map['name']

    iline += 1
    word_line = lines[iline].strip().lower()
    word = word_line.strip('*').lower()
    unused_allowed_words = ['elastic']
    unallowed_words = [
        'shell section', 'solid section', 'beam section',
        'material', 'step', 'boundary', 'amplitude', 'surface interaction',
        'assembly', 'spring', 'orientation']
    #orientations = {}
    iline += 1
    line0 = lines[iline].strip('\n\r\t, ').lower()
    #print('  wordA =', word)
    #while word in allowed_words:
    sections = {}
    density = 0.0
    ndelete = None
    ndepvars = None
    while word not in unallowed_words:
        data_lines = []
        log.info('  mat_word = %r' % word)
        #print(sections)
        word_lower = word.lower()
        if word.startswith('elastic'):
            iline, key = _read_material_elastic(
                iline, line0, lines,
                word, sections, log)
            iline += 1

        elif word.startswith('plastic'):
            key = 'plastic'
            sline = word.split(',')
            log.debug('  matword = %s' % sline)
            if len(sline) == 1:
                # elastic
                assert len(sline) in [1, 2], sline
            else:
                raise NotImplementedError(sline)
            iline, line0, data_lines = read_star_block2(iline, line0, lines, log, debug=False)
            #print(data_lines)
        elif word == 'density':
            key = 'density'
            sline = line0.split(',')
            assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
            density = float(sline[0])
            sections['density'] = [density]
            iline += 1
        elif word == 'conductivity':
            key = 'conductivity'
            sline = line0.split(',')
            assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
            conductivity = float(sline[0])
            sections['conductivity'] = [conductivity]
            iline += 1
        elif word == 'specific heat':
            key = 'specific_heat'
            sline = line0.split(',')
            assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
            specific_heat = float(sline[0])
            sections['specific_heat'] = [specific_heat]
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
            asdf
            # TODO: skips header parsing
            #iline += 1
            #line0 = lines[iline].strip().lower()
            unused_data = []
            while '*' not in line0:
                sline = line0.split(',')
                iline += 1
                line0 = lines[iline].strip().lower()
            log.debug(line0)
        elif word_lower.startswith('hyperelastic, mooney-rivlin'):
            key = 'hyperelastic, mooney-rivlin'
            iline, line0, flags, lines_out = read_generic_section(iline, word_line, lines, log)
            iline += 1
            log.debug(str((iline, line0)))

        elif word_lower.startswith('expansion'):
            iline, line0 = _read_material_expansion(iline, word_line, lines, sections, log)
        else:
            msg = print_data(lines, iline, word, 'is this an unallowed word for *Material?\n')
            raise NotImplementedError(msg)

        if key in sections:
            msg = f'  key={key!r} already defined for Material name={name!r}'
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
        log.debug(f'{iline}: word_line={word_line!r}; word={word!r}')

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
    #assert 'elastic' in sections or 'engineering constants' in sections, sections
    material = Material(name, sections=sections,
                        #is_elastic=True,
                        density=density,
                        #conductivity=conductivity, specific_heat=specific_heat,
                        ndepvars=ndepvars, ndelete=ndelete)
    iline -= 1
    return iline, line0, word, material

def _read_material_elastic(iline: int,
                           line0: str,
                           lines: list[str],
                           word: str,
                           sections: dict[str, Any],
                           log: SimpleLogger) -> tuple[int, str]:
    """reads a *Material/*Elastic block"""
    key = 'elastic'
    sword = word.split(',')
    #print(key, sword)

    #log.debug('  matword = %s' % sword)
    if len(sword) == 1:
        # elastic
        mat_type = 'isotropic'
        assert len(sword) in [1, 2], sword
    else:
        mat_type = sword[1]
        assert 'type' in mat_type, sword
        mat_type = mat_type.split('=')[1]

    sline = line0.split(',')
    if mat_type == 'traction':
        assert len(sline) == 3, sline
        log.debug('  traction material')
    elif mat_type in ('iso', 'isotropic'):
        #, TYPE=ISO
        #1.00000E+07, 3.00000E-01
        assert len(sline) == 2, sline
        e_str, nu_str = sline
        e = float(e_str)
        nu = float(nu_str)
        sections['elastic'] = [e, nu]
        #print(sections)
    elif mat_type == 'engineering constants':
        # https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node193.html
        #  two lines:
        # - E_1, E_2, E_3, nu12, nu13, nu23, G12, G13,
        # - G23, Temperature
        e1, e2, e3, nu12, nu13, n23, g12, g13 = [float(val) if val else 0. for val in sline]
        iline += 1
        sline = lines[iline].split(',')
        g23, tref = [float(val) if val else 0. for val in sline]
        key = 'engineering constants'
        sections['engineering constants'] = [e1, e2, e3, nu12, nu13, n23, g12, g13,
                                             g23, tref]
    else:
        raise NotImplementedError(f'mat_type={mat_type!r}')
    return iline, key

def read_boundary(iline: int, line0: str, lines: list[str]) -> tuple[int, str, Boundary]:
    assert isinstance(iline, int), iline
    assert isinstance(line0, str), line0
    boundary_lines = []
    line0 = lines[iline]
    if len(line0) and line0[0] == '*':
        boundary = None
        return iline, line0, boundary

    assert '*' not in line0, line0
    while '*' not in line0:
        # nid, dof1, dof2, disp
        #1,1,,0
        sline = line0.strip().split(',')
        boundary_lines.append(sline)
        iline += 1
        line0 = lines[iline].strip().lower()
    boundary = Boundary.from_data(boundary_lines)
    return iline, line0, boundary

def read_solid_section(iline: int, line0: str, lines: list[str],
                       log: SimpleLogger) -> tuple[int, SolidSection]:
    """reads *solid section"""
    assert isinstance(iline, int), iline
    assert isinstance(line0, str), line0

    #print(line0)
    iline, line0, flags, lines_out = read_generic_section(
        iline, line0, lines, log, require_lines_out=False)
    params_map = {}
    for key, value in split_strict_flags(flags):
        if key == 'material':
            params_map[key] = value.lower()
        elif key == 'elset':
            params_map[key] = value.lower()
        elif key == 'controls':
            params_map[key] = value.lower()
        #elif key == 'offset':
            #params_map[key] = float(value)
        else:  # pragma: no cover
            raise RuntimeError((key, value))
    solid_section = SolidSection.add_from_data_lines(params_map, lines_out, log)
    return iline, solid_section

def read_shell_section(iline: int, line0: str, lines: list[str],
                       log: SimpleLogger) -> tuple[int, ShellSection]:
    """reads *shell section"""
    assert isinstance(iline, int), iline
    assert '*shell' in line0, line0

    iline, line0, flags, lines_out = read_generic_section(iline, line0, lines, log)
    is_composite = 'composite' in flags
    if is_composite:
        flags.remove('composite')

    params_map = {
        'is_composite': is_composite,
        'orientation': '',
    }
    for key, value in split_strict_flags(flags):
        if key == 'material':
            params_map[key] = value.lower()
        elif key == 'elset':
            params_map[key] = value.lower()
        elif key == 'orientation':
            params_map[key] = value.lower()
        elif key == 'offset':
            params_map[key] = float(value)
        else:  # pragma: no cover
            raise RuntimeError(key)

    shell_section = ShellSection.add_from_data_lines(params_map, lines_out, log)
    return iline, shell_section

def read_hourglass_stiffness(iline: int, line0: str, lines: list[str],
                             log: SimpleLogger) -> tuple[int, Any]:
    """reads *hourglass stiffness"""
    # TODO: skips header parsing
    #iline += 1
    word2 = line0.strip('*').lower()
    iline += 1
    #params_map = get_param_map(iline, word2, required_keys=['material'])
    #log.debug('    param_map = %s' % params_map)
    #line0 = lines[iline].strip().lower()
    iline, line0, data_lines = read_star_block2(iline, line0, lines, log)
    assert len(data_lines) == 1, data_lines
    #for line in data_lines:
        #print(line)
    #solid_section = SolidSection(params_map, data_lines, log)
    hourglass_stiffness = None
    return iline, hourglass_stiffness


def read_star_block(iline: int, line0: str, lines: list[str],
                    log: SimpleLogger, debug: bool=False) -> tuple[int, str,
                                                                   list[tuple[str, ...]]]:
    """
    because this uses file streaming, there are 30,000 places where a
    try except block is needed, so this should probably be used all
    over.
    """
    assert isinstance(iline, int), iline
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
    return iline, line0, data_lines


def read_star_block2(iline: int, line0: str,
                     lines: list[str], log: SimpleLogger,
                     debug: bool=False) -> tuple[int, str, list[str]]:
    """
    because this uses file streaming, there are 30,000 places where a
    try except block is needed, so this should probably be used all
    over.
    """
    assert isinstance(iline, int), iline
    assert isinstance(line0, str), line0
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
    return iline, line0, data_lines

def _generate_set(lines_out: list[str], generate: bool) -> np.ndarray:
    set_ids_list = []
    for line in lines_out:
        set_idsi = line.strip(', \n').split(',')
        set_ids_list.extend(set_idsi)

    if generate:
        assert len(set_ids_list) == 3, set_ids_list
        # +1 is becausee it's an inclusive set
        set_ids = np.arange(int(set_ids_list[0]), int(set_ids_list[1])+1, int(set_ids_list[2]))
    elif len(set_ids_list) == 1 and not set_ids_list[0].isdigit():
        set_ids = set_ids_list[0]
    else:
        try:
            set_ids = np.unique(np.array(set_ids_list, dtype='int32'))
        except ValueError:
            print(set_ids_list)
            raise
    return set_ids

def read_orientation(iline: int, line0: str, lines: list[str],
                     log: SimpleLogger) -> tuple[int, str, Orientation]:
    assert isinstance(iline, int)

    iline, line, flags, lines_out = read_generic_section(iline, line0, lines, log)
    #assert len(flags) == 2, flags

    system = 'rectangular'
    definition = 'coordinates'
    name = ''
    for key, value in split_strict_flags(flags):
        if key == 'system':
            system = value
            assert system in ('cylindrical'), system
        elif key == 'name':
            name = value
        elif key == 'name':
            definition = value
        else:  # pramga: no cover
            raise RuntimeError((key, value))

    if system in {'rectangular', 'r'}:
        system = 'rectangular'

    assert definition == 'coordinates', definition
    if system == 'rectangular':
        assert len(lines_out) in {1, 2}, lines_out
        sline1 = lines_out[0].split(',')
        if len(sline1) == 6:
            i1, i2, i3, j1, j2, j3 = [float(val) for val in sline1]
            origin = np.zeros(3)
        elif len(sline1) == 9 and system in {'rectangular', 'z rectangular'}:
            i1, i2, i3, j1, j2, j3, origin1, origin2, origin3 = [float(val) for val in sline1]
            origin = np.array([origin1, origin2, origin3])
        else:
            raise RuntimeError(sline1)
        x_axis = np.array([i1, i2, i3])
        xy_plane = np.array([j1, j2, j3])

        axis = None
        alpha = None
        if len(lines_out) == 2:
            sline2 = lines_out[1].split(',')
            assert len(sline2) == 2, sline2
            axis, alpha = [float(val) for val in sline2]
    else:
        raise RuntimeError(system)

    #if system == ''
    #orientation_fields = []
    #for line in lines_out:
        #sline = line.split(',')
        #orientation_fields.append(sline)
    orientation = Orientation(name, system, origin,
                              x_axis=x_axis, xy_plane=xy_plane,
                              axis=axis, alpha=alpha)
    return iline, line0, orientation

def read_system(iline: int, line0: str, lines: list[str],
                log: SimpleLogger) -> tuple[int, str, list[str]]:
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

def read_transform(iline: int, line0: str, lines: list[str],
                   log: SimpleLogger) -> tuple[int, str, Transform]:
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
    iline += 1
    iline, line_out, flags, lines_out = read_generic_section(iline, line0, lines, log)

    transform_type = ''
    nset = ''
    for key, value in split_strict_flags(flags):
        if key == 'type':
            transform_type = value.upper()
        elif key == 'nset':
            nset = value
        else:  # pragma: no cover
            raise RuntimeError((key, value))

    #params = get_param_map(iline, line0, required_keys=['type', 'nset'])
    #transform_type = params['type'].upper()
    #nset = params['nset']
    assert transform_type in ['R', 'C', 'S'], f'transform_type={transform_type!r}'

    transform_fields = []
    for line in lines_out:
        sline = line.split(',')
        transform_fields.extend(sline)

    assert len(transform_fields) == 6, transform_fields
    transform = Transform.from_data(transform_type, nset, transform_fields)
    return iline, line0, transform

def _read_material_expansion(iline: int, word_line: str, lines: list[str],
                             sections, log: SimpleLogger) -> tuple[int, str]:
    """
    *Expansion, zero=20.
    80.,
    """
    iline, line0, flags, lines_out = read_generic_section(iline, word_line, lines, log)
    iline += 1

    tref = 0.
    for key, value in split_strict_flags(flags):
        if key == 'zero':
            tref = float(value)
        else:  # pragma: no cover
            raise NotImplementedError(key)

    assert len(lines_out) == 1, lines_out
    for line in lines_out:
        sline = line.split(',')
        alpha = float(sline[0])

    sections['expansion'] = [tref, alpha]
    del tref, alpha
    #log.debug(line0)
    return iline, line0

def read_tie(iline: int, line0: str, lines: list[str],
             log: SimpleLogger) -> tuple[int, str, Tie]:
    """
    '*tie, name=top_to_hub'
    ['Internal_Selection-1_Top_to_Hub_Slave, Internal_Selection-1_Top_to_Hub_Master\n']
    https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node251.html
    """
    iline, line_out, flags, lines_out = read_generic_section(iline, line0, lines, log)
    assert len(lines_out), lines_out

    if len(flags) == 0:
        raise RuntimeError("looking for element_type (e.g., '*Element, type=R2D2')\n"
                           "line0=%r\nsline=%s; allowed:\n[%s]" % (
                               line0, flags, ', '.join(allowed_element_types)))

    name = ''
    position_tolerance = np.nan
    #elset = ''
    for flag_value in flags:
        flag, value = flag_value.split('=')
        flag = flag.strip().lower()
        value = value.strip()
        if flag == 'name':
            #assert etype == '', etype
            name = value
        #elif flag == 'elset':
            #elset = value
        else:
            raise RuntimeError(flag_value)

    #elements = []
    assert len(lines_out) == 1, lines_out
    for line in lines_out:
        slave, master = line.strip('\n\t ,').split(',')
        master = master.strip().lower()
        slave = slave.strip().lower()
        #elements.append(element)

    tie = Tie(name, master, slave, position_tolerance)
    return iline, line0, tie

def read_beam_section(iline: int, line0: str, lines: list[str],
                      log: SimpleLogger) -> tuple[int, str, BeamSection]:
    """
    *BEAM SECTION, ELSET=M0B0RstdD0, MATERIAL=MaterialSolid, SECTION=RECT
    10,25
    6.1E-17, 1, 0
    https://docs.software.vt.edu/abaqusv2022/English/SIMACAEKEYRefMap/simakey-r-beamsection.htm

    First line
    ----------
    1. Beam section geometric data.  Values should be given as specified in Beam
       cross-section library for the chosen section type.  For solid circular
       sections only, you can specify a distribution name to define a spatial
       distribution for the beam radius.
    2. Etc.

    Second line (optional; enter a blank line if the default values are to be used)
    -----------
    1. First direction cosine of the first beam section axis.
    2. Second direction cosine of the first beam section axis.
    3. Third direction cosine of the first beam section axis.

    The entries on this line must be (0, 0, −1) for planar beams. The default
    for beams in space is (0, 0, −1) if the first beam section axis is not
    defined by an additional node in the element's connectivity. See Beam
    element cross-section orientation for details.

    Third line (optional)
    ---------------------
    Number of integration points in the first direction or branch.
    This number must be an odd number (for Simpson's integration),
    unless noted otherwise in Beam cross-section library.

    Number of integration points in the second direction or branch.
    This number must be an odd number (for Simpson's integration),
    unless noted otherwise in Beam cross-section library.
    This entry is needed for the THICK PIPE section, as well as for beams in space.

    Number of integration points in the third direction or branch.
    This number must be an odd number (for Simpson's integration),
    unless noted otherwise in Beam cross-section library. This entry is needed only for I-beams.

    Sections
    --------
    - ARBITRARY, for an arbitrary section.
    - BOX, for a rectangular, hollow box section.
    - CIRC, for a solid circular section.
    - HEX, for a hollow hexagonal section.
    - I, for an I-beam section.
    - L, for an L-beam section.
    - PIPE, for a thin-walled circular section.
    - RECT, for a solid, rectangular section.
    - THICK PIPE, for a thick-walled circular section (Abaqus/Standard only).
    - TRAPEZOID, for a trapezoidal section.
    """
    iline, line_out, flags, lines_out = read_generic_section(iline, line0, lines, log)
    assert len(lines_out), lines_out

    allowed_beam_section = ['ELSET', 'MATERIAL', 'SECTION']
    if len(flags) == 0:
        #ELSET=M0B0RstdD0, MATERIAL=MaterialSolid, SECTION=RECT
        raise RuntimeError("looking for beam section info (e.g., '*Beam Section, ELSET=set_name, MATERIAL=mat_name, SECTION=RECT')\n"
                           "line0=%r\nsline=%s; allowed:\n[%s]" % (
                               line0, flags, ', '.join(allowed_beam_section)))

    elset = ''
    material = ''
    section = ''
    for flag_value in flags:
        flag, value = flag_value.split('=')
        flag = flag.strip().lower()
        value = value.strip()
        #if flag == 'name':
            #assert etype == '', etype
            #name = value
        if flag == 'elset':
            assert elset == '', elset
            elset = value
        elif flag == 'material':
            assert material == '', material
            material = value
        elif flag == 'section':
            assert section == '', section
            section = value.upper()
        else:  # pragma: no cover
            raise RuntimeError(flag_value)
    assert elset is not None, elset
    assert material is not None, material
    assert section in {'ARBITRARY', 'BOX', 'CIRC', 'HEX', 'I', 'L',
                       'PIPE', 'RECT', 'THICK PIPE', 'TRAPEZOID',}, section
    assert len(lines_out) == 2, lines_out
    sline1 = lines_out[0].strip('\n\t ,').split(',')
    sline2 = lines_out[1].strip('\n\t ,').split(',')

    dimensions = np.array(sline1, dtype='float64')
    x_vector = np.array(sline2, dtype='float64')
    assert len(x_vector) == 3, x_vector

    beam_section = BeamSection(elset, material, section, dimensions, x_vector)
    return iline, line0, beam_section

def read_mass(iline: int, line0: str, lines: list[str],
              log: SimpleLogger) -> tuple[int, str, Mass]:
    """
    *Mass, elset=mass_element
    0.1
    *Mass, elset=mass_element
    0.1
    https://docs.software.vt.edu/abaqusv2022/English/SIMACAEELMRefMap/simaelm-c-masses.htm

    """
    iline, line_out, flags, lines_out = read_generic_section(iline, line0, lines, log)
    assert len(lines_out), lines_out

    allowed_mass_flags = ['ELSET']
    if len(flags) == 0:
        pass
        #ELSET=M0B0RstdD0, MATERIAL=MaterialSolid, SECTION=RECT
        raise RuntimeError("looking for beam section info (e.g., '*Mass, ELSET=set_name')\n"
                           "line0=%r\nsline=%s; allowed:\n[%s]" % (
                               line0, flags, ', '.join(allowed_mass_flags)))

    elset = ''
    for flag_value in flags:
        flag, value = flag_value.split('=')
        flag = flag.strip().lower()
        value = value.strip()
        if flag == 'elset':
            assert elset == '', elset
            elset = value
        else:  # pragma: no cover
            raise RuntimeError(flag_value)
    assert elset is not None, elset
    assert len(lines_out) == 1, lines_out

    sline = lines_out[0].strip('\n\t ,').split(',')
    assert len(sline) == 1, sline

    mass_value = float(sline[0])
    mass = Mass(elset, mass_value)

    return iline, line0, mass

def read_generic_section(iline: int, line0: str, lines: list[str],
                         log: SimpleLogger,
                         require_lines_out: bool=True) -> tuple[int, str, list[str], list[str]]:
    """
    Parameters
    ----------
    iline: int
        points to first line (not header line)
    line0: str
        the header line

    Returns
    -------
    iline: int
        pointer to line
    line : str
        is the start of the next header
    flags : list[str]
        [key, value]
    lines_out : list[str]
        the lines in the main block

    """
    assert isinstance(iline, int)
    assert isinstance(line0, str)
    iline0 = iline
    assert '*' in line0, line0
    #'*element, type=s8, elset=shell_structure' to ['type=s8', 'elset=shell_structure']
    flags = [val.strip() for val in line0.split(',')[1:]]

    line = ''
    lines_out = []
    line = lines[iline]
    while '*' not in line:
        lines_out.append(line)
        iline += 1
        line = lines[iline]

    if require_lines_out:
        assert len(lines_out), line0
    iline -= 1
    return iline, line, flags, lines_out

def read_heading(iline: int, line0: str, lines: list[str],
                         log: SimpleLogger) -> tuple[int, str, list[str]]:
    heading: list[str] = []
    if not lines[iline+1].startswith('*'):
        iline += 1
        iline, line0, flags, heading = read_generic_section(iline, line0, lines, log)
        #iline -= 1
    else:
        log.debug('empty header section')
    return iline, line0, heading

def read_part(lines: list[str], iline: int, line0: str,
              word: str,
              log: SimpleLogger,
              debug: bool) -> tuple[int, str, str, Part]:
    """reads a Part object"""
    sline2 = word.split(',', 1)[1:]

    assert len(sline2) == 1, f'looking for part_name; word={word!r} sline2={sline2}'
    name_slot = sline2[0]
    assert 'name' in name_slot, name_slot
    part_name = name_slot.split('=', 1)[1]
    log.debug(f'part_name = {part_name!r}')
    #self.part_name = part_name

    iline += 1
    line0 = lines[iline].strip().lower()
    assert line0.startswith('*node'), line0


    #iline += 1
    #line0 = lines[iline].strip().lower()

    #iline += 1
    #line0 = lines[iline].strip().lower()
    #print('line0 * = ', line0)
    element_types = {}
    node_sets: dict[str, np.ndarray] = {}
    element_sets: dict[str, np.ndarray] = {}
    #print('resetting nids...')
    nids = []
    nodes = []
    unused_is_start = True
    beam_sections: dict[str, BeamSection] = {}
    solid_sections: list[SolidSection]= []
    shell_sections: list[ShellSection] = []
    masses: dict[str, Mass] = {}
    orientations: dict[str, Orientation]= {}
    while not line0.startswith('*end part'):
        #if is_start:
        iline += 1 # skips over the header line
        log.debug('  ' + line0)
        iword = line0.strip('*').lower()
        log.info(f'part: {iword:s}')
        if '*node' in line0:
            assert len(nids) == 0, nids
            iline, line0, nids, nodes = read_node(iline, lines, log)

        elif '*element' in line0:
            #print(line0)
            iline, line0, etype, elset, elements = read_element(
                iline, line0, lines, log, debug)
            element_types[etype] = (elements, elset)
            iline += 1

        elif '*nset' in line0:
            iline, line0, set_name, set_ids = read_nset(
                iline, line0, lines, log, is_instance=False)
            node_sets[set_name] = set_ids
            iline += 1

        elif '*elset' in line0:
            iline, line0, set_name, set_ids = read_elset(
                iline, line0, lines, log, is_instance=False)
            element_sets[set_name] = set_ids
            iline += 1

        elif '*surface' in line0:
            raise RuntimeError('surface part')
            # TODO: skips header parsing
            #iline += 1
            line0 = lines[iline].strip().lower()
            data_lines = []
            while not line0.startswith('*'):
                data_lines.append(line0.split(','))
                iline += 1
                line0 = lines[iline].strip().lower()

        elif '*solid section' in line0:
            iline, solid_section = read_solid_section(iline, line0, lines, log)
            iline += 1
            solid_sections.append(solid_section)
        elif '*shell section' in line0:
            iline, shell_section = read_shell_section(iline, line0, lines, log)
            iline += 1
            shell_sections.append(shell_section)

        elif '*cohesive section' in line0:
            # TODO: skips header parsing
            #iline += 1
            #cohesive_section
            line0 = lines[iline].strip().lower()
            data_lines = []
            while not line0.startswith('*'):
                data_lines.append(line0.split(','))
                iline += 1
                line0 = lines[iline].strip().lower()
        elif '*mass' in line0:
            iline, line0, mass = read_mass(iline, line0, lines, log)
            masses[mass.elset] = mass
            # TODO: skips header parsing
            #iline, line0, flags, data_lines = reader.read_generic_section(iline, line0, lines, log)
            iline += 1
        elif '*rotary inertia' in line0:
            # TODO: skips header parsing
            iline, line0, flags, data_lines = read_generic_section(iline, line0, lines, log)
            iline += 1
        elif '*orientation' in line0:
            iline, line0, orientation = read_orientation(iline, line0, lines, log)
            orientations[orientation.name] = orientation
        else:
            msg = f'line={line0!r}\n'
            allowed = ['*node', '*element', '*nset', '*elset', '*surface',
                       '*solid section', '*cohesive section']
            msg += 'expected=[%r]' % ', '.join(allowed)
            raise NotImplementedError(msg)

        line0 = lines[iline].strip().lower()
        unused_is_start = False

        #print(line0)
    #node_sets = []
    #element_sets = []

    if debug:
        log.debug('part_name = %r' % part_name)
    #print('part.shell_sections =', shell_sections)

    del masses, orientations
    part = Part(part_name, nids, nodes, element_types,
                node_sets, element_sets,
                beam_sections, solid_sections, shell_sections, log)
    return iline, line0, part_name, part

def read_step(lines: list[str], iline: int, line0: str, istep: int,
              log: SimpleLogger) -> tuple[int, str, Step]:
    """reads a step object"""
    log.debug(f'  start of step {istep:d}...')

    boundaries = []
    node_output = []
    element_output = []
    cloads = []
    dloads = []
    surfaces = []
    frequencies = []
    # case 1
    # ------
    # *Step, name=Step-1, nlgeom=NO, inc=10000
    # *Static
    # 0.01, 1., 1e-05, 0.01
    #
    # case 2
    # ------
    #*STEP, NLGEOM=YES, AMPLITUDE=RAMP, INC=10000
    # Increase from T=117.0C to T=122.0C over 300.0 seconds  (1C/min)
    # *Static
    # 0.01, 1., 1e-05, 0.01
    #
    # case 3
    # ------
    # *STEP
    # *STATIC
    # *CLOAD
    # LOAD,2,-25

    iline += 1
    line0 = lines[iline].strip().lower()
    step_name = ''
    if not line0.startswith('*'):
        step_name = lines[iline].strip()
        iline += 1
        line0 = lines[iline].strip().lower()
    word = line0.strip('*').lower()


    #allowed_words = ['static', 'boundary', 'dsload', 'restart', 'output', 'node',
                     #'element output']
    #print('  word =', word)
    #print('  lineA =', line0)
    while word != 'end step':
        log.debug('    step_word = %r' % word)
        iline += 1
        line0 = lines[iline].strip().lower()
        #print('word =', word)
        #print('active_line =', line0)
        unused_data_lines = []
        if word.startswith('static'):
            sword = word.split(',')
            for word in sword[1:]:
                word = word.strip()
                assert word in {'solver=pardiso'}, f'words={sword}; word={word!r}'

            sline = line0.split(',')
            _line = lines[iline].strip().lower()
            if _line.startswith('*'):
                pass
            else:
                # 0.01, 1., 1e-05, 0.01
                #print('sline', sline, line0)
                assert len(sline) == 4, sline
                iline += 1
        elif word.startswith('restart'):
            line0 = lines[iline].strip().lower()
            word = line0.strip('*').lower()
            continue
            #print('  line_sline =', line0)
            #iline -= 1
            #line0 = lines[iline].strip().lower()
            #sline = line0.split(',')
            #assert len(sline) == 3, sline
            #iline += 1

        elif word.startswith('dsload'):
            #iline += 1
            #line0 = lines[iline].strip().lower()
            #print('  line_sline =', line0)
            sline = line0.split(',')
            assert len(sline) == 3, sline
            iline += 1
        elif word.startswith('dynamic'):
            log.debug(f'    line_sline = {line0!r}')
            #iline += 1
            #line0 = lines[iline].strip().lower()
            sline = line0.split(',')
            assert len(sline) >= 2, sline
            iline += 1
        elif word.startswith('visco'):
            iline += 1
        elif word.startswith('temperature'):
            iline -= 1
            line0 = lines[iline].strip().lower()
            iline, line0, unused_data_lines = read_star_block(
                iline, line0, lines, log, debug=True)
            iline += 1
        elif word.startswith('controls'):
            #log.debug('      controls')
            iline, line0, unused_data_lines = read_star_block(
                iline, line0, lines, log)
            iline += 1
            line0 = lines[iline].strip().lower()
            #for line in data_lines:
                #print(line)

        elif word.startswith('output'):
            line0 = lines[iline].strip().lower()
            word = line0.strip('*').lower()
            continue
        elif word == 'node output':
            node_output = []
            while '*' not in line0:
                sline = line0.split(',')
                node_output += [val.strip() for val in sline]
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('element output'):
            element_output = []
            while '*' not in line0:
                sline = line0.split(',')
                element_output += [val.strip() for val in sline]
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('contact output'):
            unused_contact_output = []
            while '*' not in line0:
                sline = line0.split(',')
                element_output += [val.strip() for val in sline]
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('boundary'):
            iline, line0, boundary = read_boundary(iline, line0, lines)
            if boundary:
                boundaries.append(boundary)
        elif word.startswith('buckle'):
            node_output = []
            while '*' not in line0:
                sline = line0.split(',')
                node_output += sline
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('cload'):
            iline, line0, cload = read_cload(iline, line0, lines, log)
            cloads.append(cload)
        elif word.startswith('dload'):
            iline, line0, dload = read_dload(iline, line0, lines, log)
            if dload:
                dloads.append(dload)

        elif word.startswith('surface'):
            iline, line0, surface = read_surface(iline, line0, lines, log)
            surfaces.append(surface)
        elif word.startswith('node print'):
            node_output = []
            while '*' not in line0:
                sline = line0.split(',')
                node_output += sline
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('node file'):
            node_output = []
            while '*' not in line0:
                sline = line0.split(',')
                node_output += [val.strip() for val in sline]
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('el file'):
            element_output = []
            while '*' not in line0:
                sline = line0.strip().split(',')
                element_output += [val.strip() for val in sline]
                iline += 1
                line0 = lines[iline].strip().lower()
        elif word.startswith('frequency'):
            iline -= 1
            line0 = lines[iline]
            iline, line0, frequency = read_frequency(iline, line0, lines, log)
            iline += 1
            frequencies.append(frequency)
        else:
            msg = print_data(lines, iline, word, 'is this an unallowed word for *Step?\n')
            raise NotImplementedError(msg)
        line0 = lines[iline].strip().lower()
        word = line0.strip('*').lower()
        #print('  lineB =', line0)
        #print('  word2 =', word)
    #iline += 1
    #iline -= 1
    step = Step(step_name, boundaries,
                node_output, element_output,
                cloads, dloads, surfaces,
                frequencies,
                is_nlgeom=False)
    log.debug('  end of step %i...' % istep)
    return iline, line0, step
