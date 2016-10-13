from __future__ import print_function
import copy
import numpy as np

def read_abaqus(abaqus_inp_filename, log=None, debug=False):
    model = Abaqus()
    model.read_abaqus_inp(abaqus_inp_filename)
    return model

def _clean_lines(lines):
    lines2 = []
    for line in lines:
        line2 = line.strip().split('**', 1)[0]
        #print(line2)
        if line2:
            lines2.append(line)
    return lines2


class Abaqus(object):
    def __init__(self, log=None, debug=False):
        self.debug = debug

    def read_abaqus_inp(self, abaqus_inp_filename):
        with open(abaqus_inp_filename, 'r') as abaqus_inp:
            lines = abaqus_inp.readlines()

        lines = _clean_lines(lines)
        self.parts = {}


        ilines = []
        iline = 0
        nlines = len(lines)

        while iline < nlines:
            # not handling comments right now
            line0 = lines[iline].strip().lower()
            print(iline, line0)
            #sline = line.split('**', 1)
            #if len(sline) == 1:
                #line0 = sline[0]
                #comment = ''
            #else:
                #line0, comment = sline
                #if not line0:
                    #iline += 1
                    #continue

            #print('sline =', sline)
            if '*' in line0[0]:
                word = line0.strip('*').lower()
                print('word = %r' % word)
                if word == 'heading':
                    pass
                elif word.startswith('preprint'):
                    pass

                elif word.startswith('assembly'):
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
                              word.startswith('mpc')):
                            # TODO: skips header parsing
                            iline += 1
                            line0 = lines[iline].strip().lower()
                            data_lines = []
                            while not line0.startswith('*'):
                                data_lines.append(line0.split(','))
                                iline += 1
                                line0 = lines[iline].strip().lower()
                        else:
                            raise NotImplementedError(line0)


                elif word.startswith('part'):
                    sline2 = word.split(',', 1)[1:]
                    #aq
                    assert len(sline2) == 1, sline2
                    name_slot = sline2[0]
                    assert 'name' in name_slot, name_slot
                    part_name = name_slot.split('=', 1)[1]
                    print('part_name = %r' % part_name)
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
                            #print('  nnodes =', nnodes)
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
                            assert etype in ['r2d2', 'cpe3', 'cpe4', 'cpe4r', 'coh2d4'], etype
                            if self.debug:
                                print('etype = %r' % etype)

                            #iline += 1
                            line0 = lines[iline].strip().lower()

                            elements = []
                            while not line0.startswith('*'):
                                elements.append(line0.split(','))
                                iline += 1
                                line0 = lines[iline].strip().lower()
                            element_types[etype] = elements
                        elif '*nset' in line0:
                            # TODO: skips header parsing
                            #iline += 1
                            line0 = lines[iline].strip().lower()
                            set_ids = []
                            while not line0.startswith('*'):
                                set_ids.append(line0.split(','))
                                iline += 1
                                line0 = lines[iline].strip().lower()

                        elif '*elset' in line0:
                            # TODO: skips header parsing
                            #iline += 1
                            line0 = lines[iline].strip().lower()
                            set_ids = []
                            while not line0.startswith('*'):
                                set_ids.append(line0.split(','))
                                iline += 1
                                line0 = lines[iline].strip().lower()

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
                        print('part_name = %r' % part_name)
                    part_obj = Part(part_name, nids, nodes, element_types, node_sets, element_sets)

                    #print(part_obj)
                    self.parts[part_name] = part_obj
                    if self.debug:
                        print('-------------------------------------')
                elif 'section controls' in word:
                    # TODO: skips header parsing
                    iline += 1
                    line0 = lines[iline].strip().lower()
                    data_lines = []
                    while not line0.startswith('*'):
                        data_lines.append(line0.split(','))
                        iline += 1
                        try:
                            line0 = lines[iline].strip().lower()
                        except IndexError:
                            return
                else:
                    raise NotImplementedError(word)
            else:
                raise NotImplementedError('this shouldnt happen; line=%r' % line0)
            iline += 1

            if self.debug:
                print('')

        for part_name, part in sorted(self.parts.items()):
            print(part)


class Part(object):
    def __init__(self, name, nids, nodes, element_types, node_sets, element_sets):
        self.name = name

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

        self.r2d2 = None
        self.cpe3 = None
        self.cpe4 = None
        self.cpe4r = None
        self.coh2d4 = None

        self.r2d2_eids = None
        self.cpe3_eids = None
        self.cpe4_eids = None
        self.cpe4r_eids = None
        self.coh2d4_eids = None

        if 'r2d2' in element_types: # similar to CBAR
            elements = element_types['r2d2']
            self.r2d2 = np.array(elements, dtype='int32')
            self.r2d2_eids = self.r2d2[:, 0]

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

    def element(self, eid):
        elem = None
        if self.cpe3_eids is not None:
            ieid = np.where(eid == self.cpe3_eids)[0]
            #print('self.cpe3_eids =', self.cpe3_eids)
            print('ieid_cpe3 = %s' % ieid, len(ieid))
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
            print('ieid_coh2d4 = %s' % ieid, len(ieid))
            if len(ieid):
                ieid = ieid[0]
                etype = 'coh2d4'
                elem = self.coh2d4[ieid, :]
                return etype, ieid, elem
            else:
                print('ieid = %s' % ieid)
        return None, None, None

    def __repr__(self):
        nnodes = self.nodes.shape[0]
        n_r2d2 = 0
        n_cpe3 = 0
        n_cpe4 = 0
        n_cpe4r = 0
        n_coh2d4 = 0
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
        neids = n_r2d2 + n_cpe3 + n_cpe4 + n_cpe4r + n_coh2d4
        return ('Part(name=%r, nnodes=%i, neids=%i,\n'
                '     n_r2d2=%i, n_cpe3=%i, n_cpe4=%i, n_cpe4r=%i n_coh2d4=%i)' % (
                    self.name, nnodes, neids,
                    n_r2d2, n_cpe3, n_cpe4, n_cpe4r, n_coh2d4))

def main():
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

if __name__ == '__main__':
    main()

