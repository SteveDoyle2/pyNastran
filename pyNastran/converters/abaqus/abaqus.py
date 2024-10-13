"""Defines the Abaqus class"""
import os
from io import StringIO
from typing import Union, Optional, Any

import numpy as np
from cpylog import SimpleLogger, get_logger2
from pyNastran.converters.abaqus.abaqus_cards import (
    Assembly, Part, Elements, Step, cast_nodes,
    ShellSection, SolidSection, Surface, BeamSection,
    Mass, Boundary, Material)
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
                 debug: str | bool | None=True):
        self.debug = debug
        self.parts: dict[str, Part] = {}
        self.boundaries: dict[str, Boundary] = {}
        self.materials: dict[str, Material] = {}
        self.amplitudes: dict[str, Any] = {}
        self.assembly: Optional[Assembly] = None
        self.initial_conditions: dict[str, Any] = {}
        self.steps: list[Step] = []
        self.heading: list[str] = []
        self.preprint = None
        self.node_sets: dict[str, np.ndarray] = {}
        self.element_sets: dict[str, np.ndarray] = {}

        self.shell_sections: list[ShellSection] = []
        self.solid_sections: list[SolidSection] = []
        self.log = get_logger2(log, debug)

    def read_abaqus_inp(self, abaqus_inp_filename: str, encoding: Optional[str]=None):
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

        heading: list[str] = []
        nids = []
        nodes = []
        node_sets = {}
        element_types = {}
        element_sets = {}
        orientations = {}

        beam_sections: dict[str, BeamSection] = {}
        solid_sections = []
        shell_sections = []
        boundaries = []
        surfaces: dict[str, Surface] = {}
        steps: list[Step] = []
        ties = []
        masses: dict[str, Mass] = {}

        #for ii, linei in enumerate(lines):
            #assert isinstance(linei, str), linei

        log = self.log
        while iline < nlines:
            # not handling comments right now
            line0 = lines[iline].strip().lower()
            log.debug('%s: %r' % (iline, line0))
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
                    assert len(heading) == 0, heading
                    iline, line0, heading = reader.read_heading(iline, line0, lines, log)
                elif word.startswith('preprint'):
                    pass
                elif word == 'boundary':
                    iline += 1
                    iline, line0, boundary = reader.read_boundary(iline, line0, lines)
                    if boundary:
                        boundaries.append(boundary)
                    iline -= 1
                    line0 = lines[iline].strip().lower()

                elif word.startswith('assembly'):
                    if nassembly != 0:
                        raise RuntimeError('only one assembly can be defined...')
                    iline, line0, assembly = self.read_assembly(iline, line0, lines, word)
                    self.assembly = assembly
                    nassembly += 1

                elif word.startswith('part'):
                    iline, line0, part_name, part = reader.read_part(
                        lines, iline, line0, word, self.log, self.debug)
                    self.parts[part_name] = part
                    #print('part_name', part_name)
                    if self.debug:
                        self.log.debug('-------------------------------------')
                elif 'section controls' in word:
                    # TODO: skips header parsing
                    iline, line0, data_lines = reader.read_star_block(iline, line0, lines, log)

                #elif word.startswith('amplitude'):
                    #amplitude
                    #param_map = reader.get_param_map(iline, word)
                    #name = param_map['name']
                    #if name in self.amplitudes:
                        #raise RuntimeError(f'name={name!r} is already defined...')
                    ## TODO: skips header parsing
                    #iline += 1
                    #line0 = lines[iline].strip().lower()
                    #data_lines = []
                    #while not line0.startswith('*'):
                        #data_lines.append(line0.split(','))
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #amplitude = []
                    #for sline in data_lines[:-1]:
                        #assert len(sline) == 8, sline
                        #amplitude += sline
                    #assert len(data_lines[-1]) <= 8, sline
                    #amplitude += data_lines[-1]
                    #self.amplitudes[name] = np.array(amplitude)
                    #continue

                #elif 'include' in word:
                    #pass
                elif word.startswith('material'):
                    log.debug('start of material...')
                    iline, line0, word, material = reader.read_material(iline, word, lines, log)
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
                    iline, line0, step = reader.read_step(lines, iline, line0, istep, log)
                    steps.append(step)
                    istep += 1
                #elif word.startswith('initial conditions'):
                    ##initial_conditions
                    #data_lines, iline, line0 = reader.read_star_block(lines, iline, line0, log)
                    #for line in data_lines:
                        #log.debug(line)
                    #log.debug(f'line_end_of_IC = {line0!r}')

                #elif word.startswith('surface interaction'):
                    #unused_key = 'surface interaction'
                    ##surface_interaction
                    #unused_data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)
                #elif word.startswith('friction'):
                    #unused_key = 'friction'
                    ##friction
                    #unused_data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)
                #elif word.startswith('surface behavior'):
                    #unused_key = 'surface behavior'
                    ##surface_behavior
                    #unused_data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)
                #elif word.startswith('contact damping'):
                    #unused_key = 'contact damping'
                    ##contact_damping
                    #unused_data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)
                #elif word.startswith('contact pair'):
                    #unused_key = 'contact pair'
                    ##contact_pair
                    #unused_data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #self.log.debug(line0)
                #elif word.startswith('contact output'):
                    #key = 'contact output'
                    #data = []
                    #while '*' not in line0:
                        #sline = line0.split(',')
                        #iline += 1
                        #line0 = lines[iline].strip().lower()
                    #log.debug(line0)

                #  part...
                elif word.startswith('node'):
                    iline, line0, nidsi, nodesi = reader.read_node(
                        iline, lines, log, skip_star=True)
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
                    iline, line0, etype, elset, elements = reader.read_element(
                        iline+1, line0, lines, log, self.debug)
                    element_types[etype] = (elements, elset)
                    #iline -= 1
                    #line0 = lines[iline].strip().lower()
                    #print(f'end of element; iline={iline}')
                    #print(line0)
                    assert iline > iline0
                    #print(line0)

                elif word.startswith('nset'):
                    self.log.debug('reading nset')
                    iline += 1
                    self.log.debug(line0)
                    iline, line0, set_name, set_ids = reader.read_nset(
                        iline, line0, lines, log, is_instance=False)
                    node_sets[set_name] = set_ids
                    log.debug(f'{iline}: end of nset; line={line0}')
                    #assert iline > iline0
                elif word.startswith('elset'):
                    self.log.debug('reading elset')
                    iline += 1
                    #iline0 = iline
                    #self.log.debug(line0)
                    iline, line0, set_name, set_ids = reader.read_elset(
                        iline, line0, lines, log, is_instance=False)
                    element_sets[set_name] = set_ids
                    log.debug(f'{iline}: end of elset {set_name!r}; line={line0}')
                    #assert iline > iline0
                elif '*solid section' in line0:
                    iline += 1
                    iline, solid_section = reader.read_solid_section(
                        iline, line0, lines, log)
                    log.debug(f'solid_section = {solid_section}')
                    solid_sections.append(solid_section)
                    line0 = line0.strip().lower()
                elif '*shell section' in line0:
                    iline += 1
                    iline, shell_section = reader.read_shell_section(iline, line0, lines, log)
                    #print(shell_section)
                    shell_sections.append(shell_section)
                    line0 = line0.strip().lower()
                elif '*surface' in line0:
                    iline, line0, surface = reader.read_surface(iline, line0, lines, log)
                    surfaces[surface.name] = surface

                #elif '*hourglass stiffness' in line0:
                    #iline, hourglass_stiffness = reader.read_hourglass_stiffness(iline, line0, lines, log)
                elif '*orientation' in line0:
                    iline += 1
                    iline, line0, orientation = reader.read_orientation(iline, line0, lines, log)
                    orientations[orientation.name] = orientation
                elif '*system' in line0:
                    iline, line0, system = reader.read_system(iline, line0, lines, log)
                elif '*transform' in line0:
                    iline, line0, transform = reader.read_transform(iline, line0, lines, log)
                elif '*tie' in line0:
                    iline += 1
                    iline, line0, tie = reader.read_tie(iline, line0, lines, log)
                    ties.append(tie)

                    #iline += 1
                    #iline, line0, flags, section = reader.read_generic_section(line0, lines, iline, log)
                    #log.warning('skipping tie section')
                elif '*beam section' in line0:
                    iline += 1
                    iline, line0, beam_section = reader.read_beam_section(iline, line0, lines, log)
                    beam_sections[beam_section.elset] = beam_section
                elif '*mass' in line0:
                    iline += 1
                    iline, line0, mass = reader.read_mass(iline, line0, lines, log)
                    masses[mass.elset] = mass
                    del mass
                else:
                    raise NotImplementedError(f'word={word!r} line0={line0!r}')
                assert isinstance(iline, int), word
                wordi = word.split(',')[0]
                log.debug(f'end of main {wordi!r}; line={line0!r} iline={iline}')
            else:
                # pass
                raise NotImplementedError(f'this should not happen; last_word={word!r} line={line0!r}')
            iline += 1

            #if self.debug:
                #log.debug('')

        self.nids = None
        self.nodes = None
        if nids or nodes:
            self.nids, self.nodes = cast_nodes(nids[0], nodes[0], log)

        self.heading = heading
        self.elements = Elements(element_types, self.log)
        for etype, elset_name in self.elements.element_type_to_elset_name.items():
            if elset_name == '':
                continue
            eids = getattr(self.elements, f'{etype}_eids')
            element_sets[elset_name] = eids
            del eids
        self.ties = ties
        self.masses = masses
        self.beam_sections = beam_sections
        self.shell_sections = shell_sections
        self.solid_sections = solid_sections
        self.node_sets = node_sets
        self.element_sets = element_sets
        self.orientations = orientations
        self.boundaries = boundaries
        self.surfaces = surfaces
        self.steps = steps
        log.debug('nassembly = %s' % nassembly)
        for part_name, part in sorted(self.parts.items()):
            log.info(str(part))
            part.check_materials(self.materials)
        for unused_mat_name, mat in sorted(self.materials.items()):
            log.debug(str(mat))

    def read_assembly(self, iline: int, line0: str, lines: list[str],
                      word: str) -> tuple[int, str, Assembly]:
        """reads an Assembly object"""
        assert isinstance(iline, int), iline
        assert isinstance(line0, str), line0
        log = self.log
        # TODO: skips header parsing

        iline += 1
        nlines = len(lines)
        line0 = lines[iline].strip().lower()
        element_types = {}
        node_sets = {}
        element_sets = {}
        debug = self.debug

        while not line0.startswith('*end assembly') and iline < nlines:
            log.debug('line0 assembly = %s' % line0)

            word = line0.strip('*').lower()
            log.info('assembly: %s' % word)
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
            elif (word.startswith('surface') or
                  word.startswith('rigid body') or
                  word.startswith('mpc') or
                  word.startswith('tie')):
                # TODO: skips header parsing
                iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('nset'):
                iline += 1
                iline, line0, set_name, set_ids = reader.read_nset(
                    iline, line0, lines, log, is_instance=True)
                node_sets[set_name] = set_ids
                iline += 1
            elif word.startswith('elset'):
                iline += 1
                iline, line0, set_name, set_ids = reader.read_elset(
                    iline, line0, lines, log, is_instance=True)
                element_sets[set_name] = set_ids
                iline += 1

            elif word == 'node':
                iline, line0, nids, nodes = reader.read_node(
                    iline, lines, log, skip_star=True)
            elif '*element' in line0:
                # doesn't actually start on *element line
                # 1,263,288,298,265
                iline, line0, etype, elset, elements = reader.read_element(iline, line0, lines, log, debug)
                element_types[etype] = (elements, elset)
                iline += 1
                line0 = lines[iline].strip().lower()
                #print('line_end =', line0)
            else:
                raise NotImplementedError('\nword=%r\nline=%r' % (word, line0))
            assert isinstance(iline, int), iline
            assert isinstance(line0, str), line0

        assembly = Assembly(element_types, node_sets, element_sets)
        return iline, line0, assembly

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
                raise NotImplementedError(('node_set', set_name, seti))
            for set_name, seti in self.element_sets.items():
                raise NotImplementedError(('element_set', set_name, seti))
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
            f'  # Sections:\n'
            f'    beam_sections={self.beam_sections}\n'
            f'    shell_sections={self.shell_sections}\n'
            f'    solid_sections={self.solid_sections}\n'
            f'  # Sets:\n'
            f'    node_sets={list(self.node_sets.keys())}\n'
            f'    element_sets={list(self.element_sets.keys())}\n'
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
