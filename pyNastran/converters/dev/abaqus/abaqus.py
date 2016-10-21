"""
Defines the Abaqus class
"""
from __future__ import print_function
import copy
import numpy as np
from pyNastran.utils.log import get_logger2

def read_abaqus(abaqus_inp_filename, log=None, debug=False):
    """reads an abaqus model"""
    model = Abaqus()
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


class Abaqus(object):
    """defines the abaqus reader"""
    def __init__(self, log=None, debug=True):
        self.debug = debug
        self.parts = {}
        self.boundaries = {}
        self.materials = {}
        self.amplitudes = {}
        self.assembly = {}
        self.initial_conditions = {}
        self.steps = {}
        self.heading = None
        self.preprint = None
        self.log = get_logger2(log, debug)

    def read_abaqus_inp(self, abaqus_inp_filename):
        """reads an abaqus model"""
        with open(abaqus_inp_filename, 'r') as abaqus_inp:
            lines = abaqus_inp.readlines()

        lines = _clean_lines(lines)

        ilines = []
        iline = 0
        nlines = len(lines)

        while iline < nlines:
            # not handling comments right now
            line0 = lines[iline].strip().lower()
            self.log.debug('%s, %s' % (iline, line0))
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
                #print('word1 = %r' % word)
                if word == 'heading':
                    pass
                elif word.startswith('preprint'):
                    pass
                elif word == 'boundary':
                    #print('  line_sline =', line0)
                    iline += 1
                    line0 = lines[iline].strip().lower()
                    sline = line0.split(',')
                    assert len(sline) >= 2, sline
                    #iline += 1

                elif word.startswith('assembly'):
                    iline, line0 = self.read_assembly(lines, iline, line0, word)

                elif word.startswith('part'):
                    iline, line0, part_name, part = self.read_part(lines, iline, line0, word)
                    self.parts[part_name] = part
                    if self.debug:
                        self.log.debug('-------------------------------------')
                elif 'section controls' in word:
                    # TODO: skips header parsing
                    data_lines, iline, line0 = self._read_star_block(lines, iline, line0)

                elif word.startswith('amplitude'):
                    param_map = get_param_map(word)
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
                    self.materials[material.name] = material
                    self.log.debug('end of material')

                elif word.startswith('step'):
                    self.log.debug('---------step-----------')
                    iline, line0 = self.read_step(lines, iline, line0)
                elif word.startswith('initial conditions'):
                    # TODO: skips header parsing
                    iline += 1
                    line0 = lines[iline].strip().lower()
                    data_lines = []
                    while not line0.startswith('*'):
                        data_lines.append(line0.split(','))
                        iline += 1
                        line0 = lines[iline].strip().lower()
                elif word.startswith('surface interaction'):
                    key = 'surface interaction'
                    data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('friction'):
                    key = 'friction'
                    data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('surface behavior'):
                    key = 'surface behavior'
                    data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('contact damping'):
                    key = 'contact damping'
                    data = []
                    while '*' not in line0:
                        sline = line0.split(',')
                        iline += 1
                        line0 = lines[iline].strip().lower()
                    self.log.debug(line0)
                elif word.startswith('contact pair'):
                    key = 'contact pair'
                    data = []
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

                else:
                    raise NotImplementedError(word)
            else:
                pass
                #raise NotImplementedError('this shouldnt happen; line=%r' % line0)
            iline += 1

            if self.debug:
                self.log.debug('')

        for part_name, part in sorted(self.parts.items()):
            self.log.info(part)

    def _read_star_block(self, lines, iline, line0):
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
                self.log.debug('line = %r' % line0)
            iline -= 1
            line0 = lines[iline].strip().lower()
        except IndexError:
            pass
        return data_lines, iline, line0

    def read_material(self, lines, iline, word):
        """reads a Material card"""
        param_map = get_param_map(word)
        #print(param_map)
        name = param_map['name']

        iline += 1
        line0 = lines[iline].strip().lower()
        word = line0.strip('*').lower()
        allowed_words = ['elastic']
        unallowed_words = ['material', 'step', 'boundary', 'amplitude', 'surface interaction']
        #print('  line0 =', line0)
        iline += 1
        line0 = lines[iline].strip('\n\r\t, ').lower()
        #print('  wordA =', word)
        #while word in allowed_words:
        sections = {}
        while word not in unallowed_words:
            data_lines = []
            self.log.debug('mat_word = %r' % word)
            if word.startswith('elastic'):
                key = 'elastic'
                sword = word.split(',')

                self.log.debug('  matword = %s' % sword)
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
                    else:
                        raise NotImplementedError(mat_type)
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
                iline += 1
            elif word == 'density':
                key = 'density'
                sline = line0.split(',')
                assert len(sline) == 1, 'sline=%s line0=%r' % (sline, line0)
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
                data = []
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

            elif word == 'depvar':
                key = 'depvar'
                sline = line0.split(',')
                assert len(sline) == 1, sline
                ndepvars = int(sline[0])
                iline += 1
            elif word.startswith('user material'):
                key = 'user material'
                words = word.split(',')[1:]
                for wordi in words:
                    assert '=' in wordi, wordi
                    mat_word, value = wordi.split('=')
                    mat_word = mat_word.strip()
                    if mat_word == 'constants':
                        nconstants = int(value)
                    elif mat_word == 'type':
                        mat_type = value.strip()
                        allowed_types = ['mechanical']
                        if not mat_type in allowed_types:
                            msg = 'mat_type=%r; allowed_types=[%s]'  % (mat_type, ', '.join(allowed_types))
                            raise NotImplementedError(msg)
                    else:
                        raise NotImplementedError('mat_word=%r' % mat_word)

                #nconstants = 111
                nlines_full = nconstants // 8
                nleftover = nconstants % 8
                mat_data = []
                for iiline in range(nlines_full):
                    sline = line0.split(',')
                    assert len(sline) == 8, 'len(sline)=%s; sline=%s' % (len(sline), sline)
                    mat_data += sline
                    iline += 1
                    line0 = lines[iline].strip('\n\r\t, ').lower()
                if nleftover:
                    sline = line0.split(',')
                    iline += 1
                    line0 = lines[iline].strip('\n\r\t, ').lower()
            else:
                raise NotImplementedError(print_data(lines, iline, word, 'is this an unallowed word for *Material?\n'))
            if key in sections:
                msg = 'key=%r already defined for Material name=%r' % (key, name)
                self.log.warning(msg)
                #raise RuntimeError(msg)
            sections[key] = data_lines

            line0 = lines[iline].strip('\n\r\t, ').lower()
            word = line0.strip('*').lower()

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
        material = Material(name, sections=sections)
        iline -= 1
        return iline, line0, material


    def read_assembly(self, lines, iline, line0, word):
        """reads an Assembly object"""
        # TODO: skips header parsing

        iline += 1
        line0 = lines[iline].strip().lower()
        while not line0.startswith('*end assembly'):
            #print('line0 assembly =', line0)

            word = line0.strip('*').lower()
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
                # TODO: skips header parsing
                params_map = get_param_map(word)
                name = params_map['nset']
                iline += 1
                line0 = lines[iline].strip().lower()
                assert 'instance' in params_map, params_map
                set_ids, iline, line0 = read_set(lines, iline, line0, params_map)
            else:
                raise NotImplementedError('\nword=%r\nline=%r' % (word, line0))
        return iline, line0

    def read_part(self, lines, iline, line0, word):
        """reads a Part object"""
        sline2 = word.split(',', 1)[1:]
        #aq
        assert len(sline2) == 1, sline2
        name_slot = sline2[0]
        assert 'name' in name_slot, name_slot
        part_name = name_slot.split('=', 1)[1]
        self.log.info('part_name = %r' % part_name)
        #asdf

        iline += 1
        line0 = lines[iline].strip().lower()
        assert line0 == '*node', line0


        #iline += 1
        #line0 = lines[iline].strip().lower()

        #iline += 1
        #line0 = lines[iline].strip().lower()
        #print('line0 * = ', line0)
        element_types = {}
        #print('resetting nids...')
        nids = []
        nodes = []
        is_start = True
        while not line0.startswith('*end part'):
            #if is_start:
            iline += 1 # skips over the header line
            if '*node' in line0:
                #print('  Node iline=%s' % iline)
                line0 = lines[iline].strip().lower()
                #print('  node line0 =', line0)
                is_failed = False
                #if len(nids) > 0:
                    #nids0 = copy.deepcopy(nids)
                    #nids = []
                    #is_failed = False

                #print('  ', line0)
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
                nnodes = len(nids)
                if is_failed:
                    msg = 'nids will overwrite nids0!\n'
                    #msg += 'nids0 = %s\n' % nids0
                    msg += 'nids = %s\n' % nids
                    raise RuntimeError(msg)

            elif '*element' in line0:
                sline = line0.split(',')[1:]
                assert len(sline) == 1, line0
                etype_sline = sline[0]
                assert 'type' in etype_sline, etype_sline
                etype = etype_sline.split('=')[1]
                assert etype in ['r2d2', 'cpe3', 'cpe4', 'cpe4r', 'coh2d4', 'c3d10h', 'cohax4',
                                 'cax3', 'cax4r'], etype
                if self.debug:
                    self.log.debug('etype = %r' % etype)

                #iline += 1
                line0 = lines[iline].strip().lower()

                elements = []
                while not line0.startswith('*'):
                    elements.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
                element_types[etype] = elements
            elif '*nset' in line0:
                params_map = get_param_map(word)
                name = params_map['name']
                line0 = lines[iline].strip().lower()
                assert 'part' in params_map, params_map
                set_ids, iline, line0 = read_set(lines, iline, line0, params_map)

            elif '*elset' in line0:
                # TODO: skips header parsing
                #iline += 1
                params_map = get_param_map(word)
                name = params_map['name']
                assert 'part' in params_map, params_map
                line0 = lines[iline].strip().lower()
                set_ids, iline, line0 = read_set(lines, iline, line0, params_map)

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
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()

            elif '*cohesive section' in line0:
                # TODO: skips header parsing
                #iline += 1
                line0 = lines[iline].strip().lower()
                data_lines = []
                while not line0.startswith('*'):
                    data_lines.append(line0.split(','))
                    iline += 1
                    line0 = lines[iline].strip().lower()
            else:
                raise NotImplementedError(line0)

            line0 = lines[iline].strip().lower()
            is_start = False

            #print(line0)
            #qqq
        node_sets = []
        element_sets = []

        if self.debug:
            self.log.debug('part_name = %r' % part_name)
        part = Part(part_name, nids, nodes, element_types, node_sets, element_sets, self.log)
        return iline, line0, part_name, part

    def read_step(self, lines, iline, line0):
        """reads a step object"""
        # TODO: skips header parsing

        iline += 2
        line0 = lines[iline].strip().lower()
        word = line0.strip('*').lower()
        #allowed_words = ['static', 'boundary', 'dsload', 'restart', 'output', 'node',
                         #'element output']
        #print('  word =', word)
        #print('  lineA =', line0)
        while word != 'end step':
            iline += 1
            line0 = lines[iline].strip().lower()
            #print('word =', word)
            #print('active_line =', line0)
            data_lines = []
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
                self.log.debug('  line_sline = %r' % line0)
                #iline += 1
                #line0 = lines[iline].strip().lower()
                sline = line0.split(',')
                assert len(sline) >= 2, sline
                iline += 1
            elif word.startswith('visco'):
                iline += 1
            elif word.startswith('temperature'):
                iline += 1
            elif word.startswith('controls'):
                iline += 1
                line0 = lines[iline].strip().lower()
                while not line0.startswith('*'):
                    iline += 1
                    line0 = lines[iline].strip().lower()

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
                contact_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    element_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            elif word.startswith('boundary'):
                node_output = []
                while '*' not in line0:
                    sline = line0.split(',')
                    node_output += sline
                    iline += 1
                    line0 = lines[iline].strip().lower()
            else:
                raise NotImplementedError('word = %r' % word)
            line0 = lines[iline].strip().lower()
            word = line0.strip('*').lower()
            #print('  lineB =', line0)
            #print('  word2 =', word)
        #iline += 1
        #iline -= 1
        return iline, line0

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

class Material(object):
    """a Material object is a series of nodes & elements (of various types)"""
    def __init__(self, name, sections):
        self.name = name
        if 'density' in sections:
            self.density = sections['density']
        if 'depvar' in sections:
            self.depvar = sections['depvar']
        if 'user_material' in sections:
            self.user_material = sections['user_material']


class Part(object):
    """a Part object is a series of nodes & elements (of various types)"""
    def __init__(self, name, nids, nodes, element_types, node_sets, element_sets, log):
        """creates a Part object"""
        self.name = name
        self.log = log

        self.nids = np.array(nids, dtype='int32')
        nnodes = len(self.nids)

        node0 = nodes[0]
        node_shape = len(node0)

        if node_shape == 3:
            self.nodes = np.array(nodes, dtype='float32')
        elif node_shape == 2:
            # abaqus is stupid and can have only x/y coordinates
            self.nodes = np.zeros((nnodes, 3), dtype='float32')
            nodes2 = np.array(nodes, dtype='float32')
            #print(nodes2.shape, self.nodes.shape)
            self.nodes[:, :2] = nodes2
        else:
            raise NotImplementedError(node0)

        # bars
        self.r2d2 = None

        # shells
        self.cpe3 = None
        self.cpe4 = None
        self.cpe4r = None
        self.coh2d4 = None
        self.cohax4 = None
        self.cax3 = None
        self.cax4r = None

        # solids
        self.c3d10h = None

        # bars
        self.r2d2_eids = None

        # shells
        self.cpe3_eids = None
        self.cpe4_eids = None
        self.cpe4r_eids = None
        self.coh2d4_eids = None
        self.cohax4_eids = None
        self.cax3_eids = None
        self.cax4r_eids = None

        # solids
        self.c3d10h_eids = None

        if 'r2d2' in element_types: # similar to CBAR
            elements = element_types['r2d2']
            self.r2d2 = np.array(elements, dtype='int32')
            self.r2d2_eids = self.r2d2[:, 0]

        # shells
        if 'cpe3' in element_types: # similar to CTRIA3
            elements = element_types['cpe3']
            self.cpe3 = np.array(elements, dtype='int32')
            self.cpe3_eids = self.cpe3[:, 0]

        if 'cpe4' in element_types: # similar to CQUAD4
            elements = element_types['cpe4']
            self.cpe4 = np.array(elements, dtype='int32')
            self.cpe4_eids = self.cpe4[:, 0]
            #print('  n_cpe4=%r' % str(self.cpe4.shape))

        if 'cpe4r' in element_types: # similar to CQUAD4
            elements = element_types['cpe4r']
            self.cpe4r = np.array(elements, dtype='int32')
            self.cpe4r_eids = self.cpe4r[:, 0]

        if 'coh2d4' in element_types:
            elements = element_types['coh2d4']
            #print(elements)
            self.coh2d4 = np.array(elements, dtype='int32')
            self.coh2d4_eids = self.coh2d4[:, 0]
            #print('  n_coh2d4=%r' % str(self.coh2d4.shape))

        if 'cohax4' in element_types:
            elements = element_types['cohax4']
            #print(elements)
            self.cohax4 = np.array(elements, dtype='int32')
            self.cohax4_eids = self.cohax4[:, 0]
            #print('  n_cohax4=%r' % str(self.cohax4.shape))

        if 'cax3' in element_types:
            elements = element_types['cax3']
            #print(elements)
            self.cax3 = np.array(elements, dtype='int32')
            self.cax3_eids = self.cax3[:, 0]
            #print('  n_cax3=%r' % str(self.cax3.shape))

        if 'cax4r' in element_types:
            elements = element_types['cax4r']
            #print(elements)
            self.cax4r = np.array(elements, dtype='int32')
            self.cax4r_eids = self.cax4r[:, 0]
            #print('  n_cax4r=%r' % str(self.cax4r.shape))

        # solids
        if 'c3d10h' in element_types: # similar to CTRIA3
            elements = element_types['c3d10h']
            self.c3d10h = np.array(elements, dtype='int32')
            self.c3d10h_eids = self.c3d10h[:, 0]

    def element(self, eid):
        """gets a specific element of the part"""
        elem = None
        # bars
        if self.r2d2_eids is not None:
            ieid = np.where(eid == self.r2d2_eids)[0]
            #print('self.cpe3_eids =', self.cpe3_eids)
            self.log.info('ieid_r2d2 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'r2d2'
                elem = self.r2d2[ieid, :]
                return etype, ieid, elem

         # shells
        if self.cpe3_eids is not None:
            ieid = np.where(eid == self.cpe3_eids)[0]
            #print('self.cpe3_eids =', self.cpe3_eids)
            self.log.debug('ieid_cpe3 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe3'
                elem = self.cpe3[ieid, :]
                return etype, ieid, elem

        if self.cpe4_eids is not None:
            ieid = np.where(eid == self.cpe4_eids)[0]
            #print('self.cpe4_eids =', self.cpe4_eids)
            #print('ieid = %s' % ieid)
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe4'
                elem = self.cpe4[ieid, :]
                return etype, ieid, elem

        if self.cpe4r_eids is not None:
            ieid = np.where(eid == self.cpe4r_eids)[0]
            #print('self.cpe4r_eids =', self.cpe4r_eids)
            #print('ieid = %s' % ieid)
            if len(ieid):
                ieid = ieid[0]
                etype = 'cpe4r'
                elem = self.cpe4r[ieid, :]
                return etype, ieid, elem

        if self.coh2d4_eids is not None:
            ieid = np.where(eid == self.coh2d4_eids)[0]
            #print('self.coh2d4_eids =', self.coh2d4_eids)
            self.log.debug('ieid_coh2d4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'coh2d4'
                elem = self.coh2d4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)

        if self.coh2d4_eids is not None:
            ieid = np.where(eid == self.coh2d4_eids)[0]
            #print('self.coh2d4_eids =', self.coh2d4_eids)
            print('ieid_coh2d4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'coh2d4'
                elem = self.coh2d4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)

        if self.cohax4_eids is not None:
            ieid = np.where(eid == self.cohax4_eids)[0]
            #print('self.cohax4_eids =', self.cohax4_eids)
            print('ieid_cohax4 = %s, %s' % (ieid, len(ieid)))
            if len(ieid):
                ieid = ieid[0]
                etype = 'cohax4'
                elem = self.cohax4[ieid, :]
                return etype, ieid, elem
            else:
                self.log.debug('ieid = %s' % ieid)
        return None, None, None

    def __repr__(self):
        """prints a summary for the part"""
        nnodes = self.nodes.shape[0]
        n_r2d2 = 0
        n_cpe3 = 0
        n_cpe4 = 0
        n_cpe4r = 0
        n_coh2d4 = 0
        n_c3d10h = 0

        n_cohax4 = 0
        n_cax3 = 0
        n_cax4r = 0
        if self.r2d2 is not None:
            n_r2d2 = self.r2d2.shape[0]
        if self.cpe3 is not None:
            n_cpe3 = self.cpe3.shape[0]
        if self.cpe4 is not None:
            n_cpe4 = self.cpe4.shape[0]
        if self.cpe4r is not None:
            n_cpe4r = self.cpe4r.shape[0]
        if self.coh2d4 is not None:
            n_coh2d4 = self.coh2d4.shape[0]
        if self.c3d10h is not None:
            n_c3d10h = self.c3d10h.shape[0]

        if self.cohax4 is not None:
            n_cohax4 = self.cohax4.shape[0]
        if self.cax3 is not None:
            n_cax3 = self.cax3.shape[0]
        if self.cax4r is not None:
            n_cax4r = self.cax4r.shape[0]

        neids = n_r2d2 + n_cpe3 + n_cpe4 + n_cpe4r + n_coh2d4 + n_c3d10h + n_cohax4 + n_cax3 + n_cax4r
        return ('Part(name=%r, nnodes=%i, neids=%i,\n'
                '     n_r2d2=%i, n_cpe3=%i, n_cpe4=%i, n_cpe4r=%i, n_coh2d4=%i,\n'
                '     n_cohax4=%i, n_cax3=%i, n_cax4r=%i,\n'
                '     n_c3d10h=%i)' % (
                    self.name, nnodes, neids,
                    n_r2d2, n_cpe3, n_cpe4, n_cpe4r, n_coh2d4,
                    n_cohax4, n_cax3, n_cax4r,
                    n_c3d10h))

def get_param_map(word):
    """get the optional arguments on a line"""
    words = word.split(',')
    param_map = {}
    for wordi in words:
        if '=' not in wordi:
            key = wordi.strip()
            value = None
        else:
            sword = wordi.split('=')
            assert len(sword) == 2, sword
            key = sword[0].strip()
            value = sword[1].strip()
        param_map[key] = value
    return param_map

def print_data(lines, iline, word, msg, nlines=20):
    msg = 'word=%r\n%s\n' % (word, msg)
    iline_start = iline - nlines
    iline_start = max(iline_start, 0)
    for iiline in range(iline_start, iline):
        msg += lines[iiline]
    raise NotImplementedError(msg)

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

    nids = part.nids - 1
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
    print('dedge =', dedge)
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

