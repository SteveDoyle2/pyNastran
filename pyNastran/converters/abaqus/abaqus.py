"""Defines the Abaqus class"""
import os
from io import StringIO
from typing import Union, Optional

import numpy as np
from cpylog import SimpleLogger, get_logger2
from pyNastran.converters.abaqus.abaqus_cards import (
    Assembly, Part, Elements, Step, cast_nodes,
    ShellSection, SolidSection)
import pyNastran.converters.abaqus.reader as reader
from pyNastran.converters.abaqus.reader_utils import print_data, clean_lines


def read_abaqus(abaqus_inp_filename, encoding=None,
                log=None, debug=False):
    """reads an abaqus model"""
    model = Abaqus(log=log, debug=debug)
    model.read_abaqus_inp(abaqus_inp_filename, encoding=encoding)
    return model


class Abaqus:
    """defines the abaqus reader"""
    def __init__(self, log: Optional[SimpleLogger]=None,
                 debug: Union[str, bool, None]=True):
        self.debug = debug
        self.parts = {}
        self.boundaries = {}
        self.materials = {}
        self.amplitudes = {}
        self.assembly = None
        self.initial_conditions = {}
        self.steps: list[Step] = []
        self.heading = None
        self.preprint = None
        self.node_sets = {}
        self.element_sets = {}

        self.shell_sections: list[ShellSection] = []
        self.solid_sections: list[SolidSection] = []
        self.log = get_logger2(log, debug)

    def read_abaqus_inp(self, abaqus_inp_filename: str, encoding: str=None):
        """reads an abaqus model"""
        if isinstance(abaqus_inp_filename, str):
            with open(abaqus_inp_filename, 'r', encoding=encoding) as abaqus_inp:
                lines = abaqus_inp.readlines()
        elif isinstance(abaqus_inp_filename, list):
            lines = abaqus_inp_filename
        elif isinstance(abaqus_inp_filename, StringIO):
            lines = abaqus_inp_filename.readlines()
        else:
            msg = 'abaqus_inp_filename=%s type=%r' % (
                abaqus_inp_filename, type(abaqus_inp_filename))
            raise NotImplementedError(msg)

        lines = clean_lines(lines)
        write_clean_lines = False
        if write_clean_lines:  # pragma: no cover
            dirname = os.path.dirname(abaqus_inp_filename)
            with open(os.path.join(dirname, 'spike.out'), 'w') as file_obj:
                for line in lines:
                    file_obj.write(line)

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
        boundaries = []
        surfaces = []
        steps = []

        log = self.log
        while iline < nlines:
            # not handling comments right now
            line0 = lines[iline].strip().lower()
            self.log.debug('%s: %s' % (iline, line0))
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
                    boundary, iline, line0 = reader.read_boundary(lines, line0, iline)
                    if boundary:
                        boundaries.append(boundary)
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
                    data_lines, iline, line0 = reader.read_star_block(lines, iline, line0, self.log)

                elif word.startswith('amplitude'):
                    param_map = reader.get_param_map(iline, word)
                    name = param_map['name']
                    if name in self.amplitudes:
                        raise RuntimeError(f'name={name!r} is already defined...')
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
                    log.debug('start of material...')
                    iline, line0, material = reader.read_material(lines, iline, word, log)
                    if material.name in self.materials:
                        msg = 'material.name=%r is already defined...\n' % material.name
                        msg += 'old %s' % self.materials[material.name]
                        msg += 'new %s' % material
                        raise RuntimeError(msg)
                    self.materials[material.name] = material
                    log.debug('end of material')
                #elif word.startswith('spring'):
                    #log.debug('start of spring...')
                    #iline, line0, material = read_spring(lines, iline, word, log)
                    #asdf
                    #if material.name in self.materials:
                        #msg = 'material.name=%r is already defined...\n' % material.name
                        #msg += 'old %s' % self.materials[material.name]
                        #msg += 'new %s' % material
                        #raise RuntimeError(msg)
                    #self.materials[material.name] = material
                    #log.debug('end of spring')

                elif word.startswith('step'):
                    #print('step!!!!!!!')
                    iline, line0, step = self.read_step(lines, iline, line0, istep)
                    steps.append(step)
                    istep += 1
                elif word.startswith('initial conditions'):
                    data_lines, iline, line0 = reader.read_star_block(lines, iline, line0, self.log, )
                    for line in data_lines:
                        self.log.debug(line)
                    self.log.debug(f'line_end_of_IC = {line0!r}')
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
                    iline, line0, nidsi, nodesi = reader.read_node(
                        lines, iline, self.log, skip_star=True)
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
                    line0, iline, etype, elements = reader.read_element(
                        lines, line0, iline+1, self.log, self.debug)
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
                    iline, line0, set_name, set_ids = reader.read_nset(
                        lines, iline, line0, self.log, is_instance=False)
                    node_sets[set_name] = set_ids
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    self.log.info(f'end of nset; line={line0} iline={iline}')
                    assert iline > iline0
                elif word.startswith('elset'):
                    self.log.debug('reading elset')
                    iline0 = iline
                    self.log.debug(line0)
                    iline, line0, set_name, set_ids = reader.read_elset(
                        lines, iline, line0, self.log, is_instance=False)
                    element_sets[set_name] = set_ids
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                    self.log.info(f'end of elset; line={line0} iline={iline}')
                    assert iline > iline0
                elif '*solid section' in line0:
                    iline, solid_section = reader.read_solid_section(line0, lines, iline, self.log)
                    self.log.debug(f'solid_section = {solid_section}')
                    solid_sections.append(solid_section)
                    #iline -= 1
                    #line0 = lines[iline].strip().lower()
                elif '*shell section' in line0:
                    iline, shell_section = reader.read_shell_section(line0, lines, iline, self.log)
                    #print(shell_section)
                    shell_sections.append(shell_section)
                    iline -= 1
                    line0 = lines[iline].strip().lower()
                elif '*surface' in line0:
                    iline, line0, surface = reader.read_surface(line0, lines, iline, self.log)
                    surfaces.extend(surface)
                    iline -= 1
                    line0 = lines[iline].strip().lower()

                elif '*hourglass stiffness' in line0:
                    iline, hourglass_stiffness = reader.read_hourglass_stiffness(line0, lines, iline, self.log)
                elif '*orientation' in line0:
                    iline, line0, orientation = reader.read_orientation(line0, lines, iline, self.log)
                elif '*system' in line0:
                    iline, line0, system = reader.read_system(line0, lines, iline, self.log)
                elif '*transform' in line0:
                    iline, line0, transform = reader.read_transform(line0, lines, iline, self.log)
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
        self.node_sets = node_sets
        self.element_sets = element_sets
        self.boundaries = boundaries
        self.surfaces = surfaces
        self.steps = steps
        self.log.debug('nassembly = %s' % nassembly)
        for part_name, part in sorted(self.parts.items()):
            self.log.info(str(part))
            part.check_materials(self.materials)
        for unused_mat_name, mat in sorted(self.materials.items()):
            self.log.debug(str(mat))

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
                iline, line0, set_name, set_ids = reader.read_nset(
                    lines, iline, word, self.log, is_instance=True)
                node_sets[set_name] = set_ids
            elif word.startswith('elset'):
                # TODO: skips header parsing
                params_map = reader.get_param_map(iline, word, required_keys=['instance'])
                set_name = params_map['elset']
                iline += 1
                line0 = lines[iline].strip().lower()
                set_ids, iline, line0 = reader.read_set(lines, iline, line0, params_map)
                element_sets[set_name] = set_ids
            elif word == 'node':
                iline, line0, nids, nodes = reader.read_node(
                    lines, iline, self.log, skip_star=True)
            elif '*element' in line0:
                # doesn't actually start on *element line
                # 1,263,288,298,265
                line0, iline, etype, elements = reader.read_element(lines, line0, iline)
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

        assert len(sline2) == 1, f'looking for part_name; word={word!r} sline2={sline2}'
        name_slot = sline2[0]
        assert 'name' in name_slot, name_slot
        part_name = name_slot.split('=', 1)[1]
        self.log.debug(f'part_name = {part_name!r}')
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
        log = self.log
        while not line0.startswith('*end part'):
            #if is_start:
            iline += 1 # skips over the header line
            log.debug('  ' + line0)
            iword = line0.strip('*').lower()
            log.info(f'part: {iword:s}')
            if '*node' in line0:
                assert len(nids) == 0, nids
                iline, line0, nids, nodes = reader.read_node(lines, iline, log)

            elif '*element' in line0:
                #print(line0)
                line0, iline, etype, elements = reader.read_element(
                    lines, line0, iline, self.log, self.debug)
                element_types[etype] = elements

            elif '*nset' in line0:
                #print(line0)
                iline, line0, set_name, set_ids = reader.read_nset(
                    lines, iline, line0, log, is_instance=False)
                node_sets[set_name] = set_ids

            elif '*elset' in line0:
                # TODO: skips header parsing
                #iline += 1
                #print('elset: ', line0)
                params_map = reader.get_param_map(iline, line0, required_keys=['elset'])
                set_name = params_map['elset']
                line0 = lines[iline].strip().lower()
                set_ids, iline, line0 = reader.read_set(lines, iline, line0, params_map)
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
                iline, solid_section = reader.read_solid_section(line0, lines, iline, log)
                solid_sections.append(solid_section)
            elif '*shell section' in line0:
                iline, shell_section = reader.read_shell_section(line0, lines, iline, log)
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
                iline, line0, orientation_fields = reader.read_orientation(line0, lines, iline, log)
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

        if self.debug:
            self.log.debug('part_name = %r' % part_name)
        #print('part.shell_sections =', shell_sections)
        part = Part(part_name, nids, nodes, element_types, node_sets, element_sets,
                    solid_sections, shell_sections, self.log)
        return iline, line0, part_name, part

    def read_step(self, lines, iline, line0, istep):
        """reads a step object"""
        log = self.log
        log.debug(f'  start of step {istep:d}...')

        boundaries = []
        node_output = []
        element_output = []
        cloads = []
        dloads = []
        surfaces = []
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
            self.log.debug('    step_word = %r' % word)
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
                self.log.debug(f'    line_sline = {line0!r}')
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
                unused_data_lines, iline, line0 = reader.read_star_block(
                    lines, iline, line0, self.log, debug=True)
                iline += 1
            elif word.startswith('controls'):
                #self.log.debug('      controls')
                unused_data_lines, iline, line0 = reader.read_star_block(lines, iline, line0, self.log, )
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
                boundary, iline, line0 = reader.read_boundary(lines, line0, iline)
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
                iline, line0, cload = reader.read_cload(line0, lines, iline, log)
                cloads.append(cload)
            elif word.startswith('dload'):
                iline, line0, dload = reader.read_dload(line0, lines, iline, log)
                if dload:
                    dloads.append(dload)

            elif word.startswith('surface'):
                iline, line0, surface = reader.read_surface(line0, lines, iline, log)
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
                    cloads, dloads, surfaces, is_nlgeom=False)
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
            self.log.debug(f'  nparts = {len(self.parts):d}')
            self.log.debug(f'  nmaterials = {len(self.materials):d}')
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
            for set_name, seti in self.node_sets.items():
                asdf
            for set_name, seti in self.element_sets.items():
                asdf
            for step in self.steps:
                #print(step)
                #print(abq_file)
                step.write(abq_file)

    def __repr__(self) -> str:
        msg = (
            'Abaqus:\n'
            f'  parts={self.parts}\n'
            f'  boundaries={self.boundaries}\n'
            f'  materials={self.materials}\n'
            f'  amplitudes={self.amplitudes}\n'
            f'  assembly={self.assembly}\n'
            f'  initial_conditions={self.initial_conditions}\n'
            f'  steps={self.steps}\n'
            f'  shell_sections={self.shell_sections}\n'
            f'  solid_sections={self.solid_sections}\n'
        )
        return msg

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
        nidsi = nnodes + part.nids
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
        nids = np.hstack(nids)
        nodes = np.vstack(all_nodes)
        assert len(nodes) == len(nids)
    return nnodes, nids, nodes, nelements

def main(): # pragma: no cover
    """tests a simple abaqus model"""
    abaqus_inp_filename = 'mesh.inp'
    part_name = 'part-spec'
    eid = 3707

    model = read_abaqus(abaqus_inp_filename)
    part = model.parts[part_name]
    print(part)
    etype, ieid, elem = part.element(eid)
    print(f'etype={etype} ieid={ieid:d} elem={elem}')
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
