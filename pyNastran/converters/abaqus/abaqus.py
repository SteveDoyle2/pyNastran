"""Defines the Abaqus class"""
from typing import Tuple, List, Dict, Union, Optional, Any
import numpy as np
from cpylog import SimpleLogger, get_logger2
from pyNastran.converters.abaqus.abaqus_cards import (
    Assembly, Material, Part, Elements,
    Step, SolidSection, ShellSection,
    cast_nodes, allowed_element_types)


def read_abaqus(abaqus_inp_filename, log=None, debug=False):
    """reads an abaqus model"""
    model = Abaqus(log=log, debug=debug)
    model.read_abaqus_inp(abaqus_inp_filename)
    return model

def _clean_lines(lines):
    """removes comment lines and concatenates include files"""
    lines2 = []
    for line in lines:
        line2 = line.strip().split('**', 1)[0]
        #print(line2)
        if line2:
            if 'include' in line2.lower():
                sline = line2.split(',')
                assert len(sline) == 2, sline
                assert '=' in sline[1], sline
                sline2 = sline[1].split('=')
                assert len(sline2) == 2, sline2
                base, inc_filename = sline2
                base = base.strip()
                inc_filename = inc_filename.strip()
                assert base.lower() == 'input', 'base=%r' % base.lower()

                with open(inc_filename, 'r') as inc_file:
                    inc_lines = inc_file.readlines()
                inc_lines = _clean_lines(inc_lines)
                lines2 += inc_lines
                continue

            lines2.append(line)
    return lines2


class Abaqus:
    """defines the abaqus reader"""
    def __init__(self, log=None, debug: Union[str, bool, None]=True):
        self.debug = debug
        self.parts = {}
        self.boundaries = {}
        self.materials = {}
        self.amplitudes = {}
        self.assembly = None
        self.initial_conditions = {}
        self.steps = []
        self.heading = None
        self.preprint = None

        self.shell_sections = []  # List[ShellSection]
        self.solid_sections = []  # List[SolidSection]
        self.log = get_logger2(log, debug)

    def read_abaqus_inp(self, abaqus_inp_filename):
        """reads an abaqus model"""
        if isinstance(abaqus_inp_filename, str):
            with open(abaqus_inp_filename, 'r') as abaqus_inp:
                lines = abaqus_inp.readlines()
        elif isinstance(abaqus_inp_filename, list):
            lines = abaqus_inp_filename
        else:
            msg = 'abaqus_inp_filename=%s type=%r' % (
                abaqus_inp_filename, type(abaqus_inp_filename))
            raise NotImplementedError(msg)

        lines = _clean_lines(lines)

        unused_ilines = []
        iline = 0
        nlines = len(lines)
        nassembly = 0
        istep = 1

        nids = []
        nodes = []
        node_sets = {}
        element_types = {}
        element_sets = {}
        solid_sections = []
        shell_sections = []
        steps = []

        log = self.log
        while iline < nlines:
            # not handling comments right now
            line0 = lines[iline].strip().lower()
            #self.log.debug('%s, %s' % (iline, line0))
            #sline = line.split('**', 1)
            #if len(sline) == 1:
                #line0 = sline[0]
                #comment = ''
            #else:
                #line0, comment = sline
                #if not line0:
                    #iline += 1
                    #continue

            if '*' in line0[0]:
                word = line0.strip('*').lower()
                log.debug('main: word = %r' % word)
                if word == 'heading':
                    pass
                elif word.startswith('preprint'):
                    pass
                elif word == 'boundary':
                    iline += 1
                    boundary, iline, line0 = read_boundary(lines, line0, iline)
                    iline -= 1
                    line0 = lines[iline].strip().lower()

                elif word.startswith('assembly'):
                    if nassembly != 0:
                        raise RuntimeError('only one assembly can be defined...')
                    iline, line0, assembly = self.read_assembly(lines, iline, line0, word)
                    self.assembly = assembly
                    nassembly += 1

                elif word.startswith('part'):
                    iline, line0, part_name, part = self.read_part(lines, iline, line0, word)
                    self.parts[part_name] = part
                    #print('part_name', part_name)
                    if self.debug:
                        self.log.debug('-------------------------------------')
                elif 'section controls' in word:
                    # TODO: skips header parsing
                    data_lines, iline, line0 = _read_star_block(lines, iline, line0, self.log, )

                elif word.startswith('amplitude'):
                    param_map = get_param_map(iline, word)
                    name = param_map['name']
                    if name in self.amplitudes:
                        raise RuntimeError('name=%r is already defined...' % name)
                    # TODO: skips header parsing
                    iline += 1
                    line0 = lines[iline].strip().lower()
                    data_lines = []
                    while not line0.startswith('*'):
                        data_lines.append(line0.split(','))
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    amplitude = []
                    for sline in data_lines[:-1]:
                        assert len(sline) == 8, sline
                        amplitude += sline
                    assert len(data_lines[-1]) <= 8, sline
                    amplitude += data_lines[-1]
                    self.amplitudes[name] = np.array(amplitude)
                    continue
                    #iline -= 1
                    #line0 = lines[iline].strip().lower()

                #elif 'include' in word:
                    #pass
                elif word.startswith('material'):
                    self.log.debug('start of material...')
                    iline, line0, material = self.read_material(lines, iline, word)
                    if material.name in self.materials:
                        msg = 'material.name=%r is already defined...\n' % material.name
                        msg += 'old %s' % self.materials[material.name]
                        msg += 'new %s' % material
                        raise RuntimeError(msg)
                    self.materials[material.name] = material
                    self.log.debug('end of material')
                #elif word.startswith('spring'):
                    #self.log.debug('start of spring...')
                    #iline, line0, material = self.read_spring(lines, iline, word)
                    #asdf
                    #if material.name in self.materials:
                        #msg = 'material.name=%r is already defined...\n' % material.name
                        #msg += 'old %s' % self.materials[material.name]
                        #msg += 'new %s' % material
                        #raise RuntimeError(msg)
                    #self.materials[material.name] = material
                    #self.log.debug('end of spring')

                elif word.startswith('step'):
                    #print('step!!!!!!!')
                    iline, line0, step = self.read_step(lines, iline, line0, istep)
                    steps.append(step)
                    istep += 1
                elif word.startswith('initial conditions'):
                    data_lines, iline, line0 = _read_star_block(lines, iline, line0, self.log, )
                    for line in data_lines:
                        self.log.debug(line)
                    self.log.debug('line_end_of_IC = %s' % line0)
                elif word.startswith('surface interaction'):
                    unused_key = 'surface interaction'
                    unused_data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('friction'):
                    unused_key = 'friction'
                    unused_data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('surface behavior'):
                    unused_key = 'surface behavior'
                    unused_data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('contact damping'):
                    unused_key = 'contact damping'
                    unused_data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('contact pair'):
                    unused_key = 'contact pair'
                    unused_data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                #elif word.startswith('contact output'):
                    #key = 'contact output'
                    #data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)

                #  part...
                elif word.startswith('node'):
                    iline, line0, nidsi, nodesi = read_node(lines, iline, self.log, skip_star=True)
                    nids.append(nidsi)
                    nodes.append(nodesi)
                    #print(f'end of node; iline={iline}')
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    #print(line0)
                elif '*element' in line0:
                    # line0: *ELEMENT,TYPE=C3D4
                    # iline: doesn't start on *element line
                    # 1,263,288,298,265
                    #print(f'start of element; iline={iline}')
                    iline0 = iline
                    line0 = lines[iline].strip().lower()
                    line0, iline, etype, elements = self._read_elements(lines, line0, iline+1)
                    element_types[etype] = elements
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    #print(f'end of element; iline={iline}')
                    #print(line0)
                    assert iline > iline0
                    #print(line0)

                elif word.startswith('nset'):
                    self.log.debug('reading nset')
                    iline0 = iline
                    self.log.debug(line0)
                    iline, line0, set_name, set_ids = read_nset(lines, iline, line0, self.log,
                                                                is_instance=False)
                    node_sets[set_name] = set_ids
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    self.log.info(f'end of nset; line={line0} iline={iline}')
                    assert iline > iline0
                elif word.startswith('elset'):
                    self.log.debug('reading elset')
                    iline0 = iline
                    self.log.debug(line0)
                    iline, line0, set_name, set_ids = read_elset(lines, iline, line0, self.log,
                                                                is_instance=False)
                    node_sets[set_name] = set_ids
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    self.log.info(f'end of elset; line={line0} iline={iline}')
                    assert iline > iline0
                elif '*solid section' in line0:
                    iline, solid_section = read_solid_section(line0, lines, iline, self.log)
                    self.log.debug(f'solid_section = {solid_section}')
                    solid_sections.append(solid_section)
                    #iline -= 1
                    #line0 = lines[iline].strip().lower()
                elif '*shell section' in line0:
                    iline, shell_section = read_shell_section(line0, lines, iline, self.log)
                    #print(shell_section)
                    shell_sections.append(shell_section)
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                elif '*hourglass stiffness' in line0:
                    iline, hourglass_stiffness = read_hourglass_stiffness(line0, lines, iline, self.log)
                elif '*orientation' in line0:
                    unused_key = 'orientation'
                    unused_data = []
                    iline += 1
                    line0 = lines[iline].strip().lower()
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                    iline -= 1
                    line0 = lines[iline].strip().lower()

                else:
                    raise NotImplementedError(f'word={word!r} line0={line0!r}')
                wordi = word.split(',')[0]
                self.log.info(f'end of main {wordi!r}; line={line0!r} iline={iline}')
            else:
                # pass
                raise NotImplementedError(f'this should not happen; last_word={word!r} line={line0!r}')
            iline += 1

            #if self.debug:
                #self.log.debug('')

        self.nids = None
        self.nodes = None
        if nids or nodes:
            self.nids, self.nodes = cast_nodes(nids[0], nodes[0], self.log)
        self.elements = Elements(element_types, self.log)
        self.element_sets = element_sets
        self.shell_sections = shell_sections
        self.solid_sections = solid_sections
        self.steps = steps
        self.log.debug('nassembly = %s' % nassembly)
        for part_name, part in sorted(self.parts.items()):
            self.log.info(str(part))
            part.check_materials(self.materials)
        for unused_mat_name, mat in sorted(self.materials.items()):
            self.log.debug(str(mat))

    def read_spring(self, lines: List[str], iline: int, word: str) -> Material:
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

    def read_material(self, lines: List[str], iline: int, word: str) -> Material:
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
            #self.log.info('  mat_word = %r' % word)
            #print(sections)
            if word.startswith('elastic'):
                key = 'elastic'
                sword = word.split(',')
                #print(key, sword)

                #self.log.debug('  matword = %s' % sword)
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
                        self.log.debug('  traction material')
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
                self.log.debug('  matword = %s' % sword)
                if len(sword) == 1:
                    # elastic
                    assert len(sline) in [1, 2], sline
                else:
                    raise NotImplementedError(sline)
                data_lines, iline, line0 = _read_star_block2(lines, iline, line0, self.log, debug=False)
                #print(data_lines)
            elif word == 'density':
                key = 'density'
                sline = line0.split(',')
                assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
                density = float(sline[0])
                iline += 1
            elif word.startswith('damage initiation'):
                key = 'damage initiation'
                #self.log.debug('  damage0 %s' % line0)
                sline = line0.split(',')
                self.log.debug(sline)
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
                self.log.debug(line0)
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
                #self.log.debug(line0)
            #elif word.startswith('friction'):
                #key = 'friction'
                #data = []
                #while '*' not in line0:
                    #sline = line0.split(',')
                    #iline += 1
                    #line0 = lines[iline].strip().lower()
                #self.log.debug(line0)
            #elif word.startswith('surface behavior'):
                #key = 'surface behavior'
                #data = []
                #while '*' not in line0:
                    #sline = line0.split(',')
                    #iline += 1
                    #line0 = lines[iline].strip().lower()
                #self.log.debug(line0)
            #elif word.startswith('contact damping'):
                #key = 'contact damping'
                #data = []
                #while '*' not in line0:
                    #sline = line0.split(',')
                    #iline += 1
                    #line0 = lines[iline].strip().lower()
                #self.log.debug(line0)

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
                self.log.debug(line0)
            elif word.lower().startswith('hyperelastic, mooney-rivlin'):
                key = 'hyperelastic, mooney-rivlin'
                while '*' not in line0:
                    sline = line0.split(',')
                    iline += 1
                    line0 = lines[iline].strip().lower()
                self.log.debug(line0)
            elif word.lower().startswith('expansion'):
                #*Expansion, zero=20.
                #80.,
                key = 'expansion'
                while '*' not in line0:
                    sline = line0.split(',')
                    iline += 1
                    line0 = lines[iline].strip().lower()
                #iline += 1
                self.log.debug(line0)
            else:
                msg = print_data(lines, iline, word, 'is this an unallowed word for *Material?\n')
                raise NotImplementedError(msg)

            if key in sections:
                msg = f'key={key!r} already defined for Material name={name!r}'
                self.log.warning(msg)
            else:
                #raise RuntimeError(msg)
                sections[key] = data_lines

            try:
                line = lines[iline]
            except IndexError:
                is_broken = True
                self.log.debug('  breaking on end of file')
                break
            word_line = line.strip('\n\r\t, ').lower()
            del line
            word = word_line.strip('*').lower()

            iline += 1
            line0 = lines[iline].strip('\n\r\t, ').lower()
            #self.log.debug('  lineB = %r' % line0)
            #self.log.debug('  wordB = %r' % word)

            is_broken = False
            for unallowed_word in unallowed_words:
                if word.startswith(unallowed_word):
                    self.log.debug('  breaking on %r' % unallowed_word)
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


    def read_assembly(self, lines, iline, line0, word):
        """reads an Assembly object"""
        # TODO: skips header parsing

        iline += 1
        nlines = len(lines)
        line0 = lines[iline].strip().lower()
        element_types = {}
        node_sets = {}
        element_sets = {}

        while not line0.startswith('*end assembly') and iline < nlines:
            self.log.debug('line0 assembly = %s' % line0)

            word = line0.strip('*').lower()
            self.log.info('assembly: %s' % word)
            if '*instance' in line0:
                # TODO: skips header parsing
                iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
                assert line0.startswith('*end instance'), line0
                iline += 1
                line0 = lines[iline].strip().lower()
            elif (word.startswith('surface') or word.startswith('rigid body') or
                  word.startswith('mpc') or word.startswith('tie')):
                # TODO: skips header parsing
                iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('nset'):
                iline, line0, set_name, set_ids = read_nset(lines, iline, word, self.log,
                                                            is_instance=True)
                node_sets[set_name] = set_ids
            elif word.startswith('elset'):
                # TODO: skips header parsing
                params_map = get_param_map(iline, word, required_keys=['instance'])
                set_name = params_map['elset']
                iline += 1
                line0 = lines[iline].strip().lower()
                set_ids, iline, line0 = read_set(lines, iline, line0, params_map)
                element_sets[set_name] = set_ids
            elif word == 'node':
                iline, line0, nids, nodes = read_node(lines, iline, self.log, skip_star=True)
            elif '*element' in line0:
                # doesn't actually start on *element line
                # 1,263,288,298,265
                line0, iline, etype, elements = self._read_elements(lines, line0, iline)
                element_types[etype] = elements
                iline += 1
                line0 = lines[iline].strip().lower()
                #print('line_end =', line0)
            else:
                raise NotImplementedError('\nword=%r\nline=%r' % (word, line0))

        assembly = Assembly(element_types, node_sets, element_sets)
        return iline, line0, assembly

    def read_part(self, lines, iline, line0, word):
        """reads a Part object"""
        sline2 = word.split(',', 1)[1:]

        assert len(sline2) == 1, 'looking for part_name; word=%r sline2=%s' % (word, sline2)
        name_slot = sline2[0]
        assert 'name' in name_slot, name_slot
        part_name = name_slot.split('=', 1)[1]
        self.log.debug('part_name = %r' % part_name)
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
        node_sets = {}
        element_sets = {}
        #print('resetting nids...')
        nids = []
        nodes = []
        unused_is_start = True
        solid_sections = []
        shell_sections = []
        while not line0.startswith('*end part'):
            #if is_start:
            iline += 1 # skips over the header line
            self.log.debug('  ' + line0)
            iword = line0.strip('*').lower()
            self.log.info('part: %s' % iword)
            if '*node' in line0:
                assert len(nids) == 0, nids
                iline, line0, nids, nodes = read_node(lines, iline, self.log)

            elif '*element' in line0:
                #print(line0)
                line0, iline, etype, elements = self._read_elements(lines, line0, iline)
                element_types[etype] = elements

            elif '*nset' in line0:
                #print(line0)
                iline, line0, set_name, set_ids = read_nset(lines, iline, line0, self.log, is_instance=False)
                node_sets[set_name] = set_ids

            elif '*elset' in line0:
                # TODO: skips header parsing
                #iline += 1
                #print('elset: ', line0)
                params_map = get_param_map(iline, line0, required_keys=['elset'])
                set_name = params_map['elset']
                line0 = lines[iline].strip().lower()
                set_ids, iline, line0 = read_set(lines, iline, line0, params_map)
                element_sets[set_name] = set_ids

            elif '*surface' in line0:
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()

            elif '*solid section' in line0:
                iline, solid_section = read_solid_section(line0, lines, iline, self.log)
                solid_sections.append(solid_section)
            elif '*shell section' in line0:
                iline, shell_section = read_shell_section(line0, lines, iline, self.log)
                shell_sections.append(shell_section)

            elif '*cohesive section' in line0:
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif '*mass' in line0:
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif '*rotary inertia' in line0:
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif '*orientation' in line0:
                unused_key = 'orientation'
                unused_data = []
                while '*' not in line0:
                    sline = line0.split(',')
                    iline += 1
                    line0 = lines[iline].strip().lower()
                self.log.debug(line0)
            else:
                msg = 'line=%r\n' % line0
                allowed = ['*node', '*element', '*nset', '*elset', '*surface',
                           '*solid section', '*cohesive section']
                msg += 'expected=[%r]' % ', '.join(allowed)
                raise NotImplementedError(msg)

            line0 = lines[iline].strip().lower()
            unused_is_start = False

            #print(line0)
        #node_sets = []
        #element_sets = []

        if self.debug:
            self.log.debug('part_name = %r' % part_name)
        #print('part.shell_sections =', shell_sections)
        part = Part(part_name, nids, nodes, element_types, node_sets, element_sets,
                    solid_sections, shell_sections, self.log)
        return iline, line0, part_name, part

    def _read_elements(self, lines, line0, iline):
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

        if self.debug:
            self.log.debug('    etype = %r' % etype)

        #iline += 1
        line1 = lines[iline].strip().lower()
        self.log.debug('    line1 = %r' % line1)

        elements = []
        #print(line1)
        assert '*' not in line1, line1
        while not line1.startswith('*'):
            #print(line1)
            elements.append(line1.split(','))
            iline += 1
            line1 = lines[iline].strip().lower()
        #self.log.debug('elements = %s' % elements)
        return line1, iline, etype, elements

    def read_step(self, lines, iline, line0, istep):
        """reads a step object"""
        self.log.debug('  start of step %i...' % istep)

        boundaries = []
        outputs = []
        cloads = []
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
            self.log.debug('    step_word = %r' % word)
            iline += 1
            line0 = lines[iline].strip().lower()
            #print('word =', word)
            #print('active_line =', line0)
            unused_data_lines = []
            if word == 'static':
                sline = line0.split(',')
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
                self.log.debug('    line_sline = %r' % line0)
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
                unused_data_lines, iline, line0 = _read_star_block(
                    lines, iline, line0, self.log, debug=True)
                iline += 1
            elif word.startswith('controls'):
                #self.log.debug('      controls')
                unused_data_lines, iline, line0 = _read_star_block(lines, iline, line0, self.log, )
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
                    node_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('element output'):
                element_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    element_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('contact output'):
                unused_contact_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    element_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('boundary'):
                boundary, iline, line0 = read_boundary(lines, line0, iline)
            elif word.startswith('buckle'):
                node_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    node_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('cload'):
                iline, line0, cload = read_cload(line0, lines, iline)
                cloads.append(cload)
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
                    node_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('el file'):
                node_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    node_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            else:
                msg = print_data(lines, iline, word, 'is this an unallowed word for *Step?\n')
                raise NotImplementedError(msg)
            line0 = lines[iline].strip().lower()
            word = line0.strip('*').lower()
            #print('  lineB =', line0)
            #print('  word2 =', word)
        #iline += 1
        #iline -= 1
        step = Step(step_name, boundaries, outputs, cloads, is_nlgeom=False)
        self.log.debug('  end of step %i...' % istep)
        return iline, line0, step

    def write(self, abaqus_filename_out, is_2d=False):
        self.log.info('writing %r' % abaqus_filename_out)
        assert isinstance(self.steps, list), self.steps
        #self.parts = {}
        #self.boundaries = {}
        #self.materials = {}
        #self.amplitudes = {}
        #self.assembly = None
        #self.initial_conditions = {}
        #self.steps = {}
        #self.heading = None
        #self.preprint = None
        with open(abaqus_filename_out, 'w') as abq_file:
            self.log.debug("  nparts = %s" % len(self.parts))
            self.log.debug("  nmaterials = %s" % len(self.materials))
            if self.assembly is not None:
                self.assembly.write(abq_file)
            for unused_part_name, part in self.parts.items():
                part.write(abq_file, is_2d=is_2d)
            for unused_part_name, initial_conditions in self.initial_conditions.items():
                initial_conditions.write(abq_file)
            for unused_part_name, amplitude in self.amplitudes.items():
                amplitude.write(abq_file)
            for unused_mat_name, mat in self.materials.items():
                mat.write(abq_file)
            for step in self.steps:
                #print(step)
                #print(abq_file)
                step.write(abq_file)

def get_nodes_nnodes_nelements(model: Abaqus, stop_for_no_elements: bool=True):
    """helper method"""
    nnodes = 0
    nelements = 0
    nids = []
    all_nodes = []
    #if model.nodes and model.elements:
    if model.nids is not None and len(model.nids):
        nidsi = model.nids
        nodes = model.nodes
        elements = model.elements

        nnodes += nodes.shape[0]
        nelements += elements.nelements
        nids.append(nidsi)
        all_nodes.append(nodes)

    for unused_part_name, part in model.parts.items():
        #unused_nids = part.nids - 1
        nidsi = part.nids
        nodes = part.nodes
        elements = part.elements

        nnodes += nodes.shape[0]
        nelements += elements.nelements
        nids.append(nidsi)
        all_nodes.append(nodes)

    if nelements == 0 and stop_for_no_elements:
        raise RuntimeError('nelements=0')

    if len(all_nodes) == 1:
        nids = nids[0]
        nodes = all_nodes[0]
    else:
        nids = np.vstack(nids)
        nodes = np.vstack(all_nodes)
    return nnodes, nids, nodes, nelements

def read_cload(line0, lines, iline) -> Tuple[int, str, Any]:
    cload = []
    while '*' not in line0:
        sline = line0.split(',')
        assert len(sline) == 3, sline
        #cload += sline
        iline += 1
        line0 = lines[iline].strip().lower()
        nid = int(sline[0])
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

def read_elset(lines: List[str], iline: int, word: str, log: SimpleLogger,
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

def read_boundary(lines: List[str], line0: str, iline: int) -> Tuple[Any, int, str]:
    boundary = []
    line0 = lines[iline]
    assert '*' not in line0, line0
    while '*' not in line0:
        sline = line0.split(',')
        boundary += sline
        iline += 1
        line0 = lines[iline].strip().lower()
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
    data_lines, iline, line0 = _read_star_block2(lines, iline, line0, log)
    #print('line0 =', iline, line0)
    #print(f'lines[{iline}] = {lines[iline]!r}')
    #print('lines[iline+1] =', lines[iline+1])
    #print('data_lines =', data_lines)
    #for line in data_lines:
        #print(line)
    solid_section = SolidSection.add_from_data_lines(params_map, data_lines, log)
    return iline, solid_section

def read_shell_section(line0: str, lines: List[str], iline: int,
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
    data_lines, iline, line0 = _read_star_block2(lines, iline, line0, log)
    log.info(f'params_map = {params_map}')
    log.info(f'data_lines = {data_lines}')
    #for line in data_lines:
        #print(line)
    assert len(data_lines) > 0, data_lines
    shell_section = ShellSection.add_from_data_lines(params_map, data_lines, log)
    #print(lines[iline])
    return iline, shell_section

def read_hourglass_stiffness(line0: str, lines: List[str], iline: int,
                             log: SimpleLogger) -> None:
    """reads *hourglass stiffness"""
    # TODO: skips header parsing
    #iline += 1
    word2 = line0.strip('*').lower()
    iline += 1
    #params_map = get_param_map(iline, word2, required_keys=['material'])
    #log.debug('    param_map = %s' % params_map)
    #line0 = lines[iline].strip().lower()
    data_lines, iline, line0 = _read_star_block2(lines, iline, line0, log)
    assert len(data_lines) == 1, data_lines
    #for line in data_lines:
        #print(line)
    #solid_section = SolidSection(params_map, data_lines, log)
    hourglass_stiffness = None
    return iline, hourglass_stiffness


def _read_star_block(lines, iline, line0, log, debug=False):
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


def _read_star_block2(lines, iline, line0, log, debug=False):
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

def get_param_map(iline: int, word: str, required_keys: Optional[List[str]]=None) -> Dict[str, 'str']:
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

def split_by_equals(word, unused_lines, iline):
    """
    splits 'x = 42'
    into 'x' and '42'
    """
    if '=' not in word:
        msg = 'line %i: %r cannot be split by an equals sign (=)' % (iline, word)
        raise RuntimeError(msg)
    word_out, value = word.split('=')
    return word_out, value

def print_data(lines, iline, word, msg, nlines=20):
    """prints the last N lines"""
    msg = 'word=%r\n%s\n' % (word, msg)
    iline_start = iline - nlines
    iline_start = max(iline_start, 0)
    for iiline in range(iline_start, iline):
        msg += lines[iiline]
    return msg

def main(): # pragma: no cover
    """tests a simple abaqus model"""
    abaqus_inp_filename = 'mesh.inp'
    part_name = 'part-spec'
    eid = 3707

    model = read_abaqus(abaqus_inp_filename)
    part = model.parts[part_name]
    print(part)
    etype, ieid, elem = part.element(eid)
    print('etype=%s ieid=%s elem=%s' % (etype, ieid, elem))
    #return

    unused_nids = part.nids - 1
    nodes = part.nodes
    cohesive_elements = part.coh2d4
    assert cohesive_elements is not None, cohesive_elements
    n1 = cohesive_elements[:, 1] - 1
    n2 = cohesive_elements[:, 2] - 1
    #print('n1 =', n1)
    #print('n2 =', n2)
    #print('nodes =', nodes)


    #ix = np.unique(np.hstack([n2, n1]))
    ix = np.append(n2, n1[-1])
    eids = cohesive_elements[:, 0] #- cohesive_elements[0, 0]
    x = nodes[ix, 0]
    edge_length_21 = np.abs(nodes[n2, 0] - nodes[n1, 0])
    edge_length_max = edge_length_21.max()
    edge_length_min = edge_length_21.min()
    dedge = edge_length_max - edge_length_min
    #print('edge_length_21 =\n%s' % edge_length_21)

    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.suptitle(abaqus_inp_filename)
    plt.plot(eids, edge_length_21 * 1000., 'b-o')
    if dedge < 1e-6:
        plt.ylim(0.98 * edge_length_min * 1000.,
                 1.02 * edge_length_min * 1000.)
    plt.ylabel('edge length (mm)')
    plt.xlabel('element id')
    plt.grid()

    plt.figure(2)
    plt.suptitle(abaqus_inp_filename)
    plt.plot(x[:-1] * 1000., edge_length_21 * 1000., 'b-o')
    if dedge < 1e-6:
        plt.ylim(0.98 * edge_length_min * 1000.,
                 1.02 * edge_length_min * 1000.)
    plt.grid()
    plt.ylabel('edge length (mm)')
    plt.xlabel('x location (mm)')
    plt.show()

if __name__ == '__main__': # pragma: no cover
    main()
