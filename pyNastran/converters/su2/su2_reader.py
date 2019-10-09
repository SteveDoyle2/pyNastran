import numpy as np


def read_su2(su2_filename, log=None, debug=False):
    """reads a 2d/3d SU2 single zone model"""
    model = SU2Reader(log=log, debug=debug)
    unused_ndim, nodes, elements, regions = model.read_su2(su2_filename)
    #model.to_cart3d()
    return model, nodes, elements, regions

class SU2Reader:
    """SU2 reader/writer"""
    etype_nnodes_map = {
        #Line       3
        #Triangle   5
        #Quadrilateral  9
        #Tetrahedral    10
        #Hexahedral 12
        #Wedge      13
        #Pyramid    14
        3 : 2,
        5 : 3,
        9 : 4,
        10 : 4,
        12 : 8,
        13 : 6,
        14 : 5,
    }
    def __init__(self, log=None, debug=False):
        """initializes the SU2Reader object"""
        self.log = log
        self.debug = debug
        self.su2_filename = None

    def read_2d(self, lines, i, ndim, nelem):
        """reads a 2d SU2 zone"""
        # elements
        tris = []
        quads = []
        for unused_ne in range(nelem):
            # what's the 0th slot?
            #Line           3
            #Triangle       5
            #Quadrilateral  9
            #Tetrahedral    10
            #Hexahedral     12
            #Wedge          13
            #Pyramid        14
            line = lines[i]
            data = line.split()[:-1]
            etype, *nodes = data
            if etype == '9':
                assert len(nodes) == 4, nodes
                quads.append(nodes)
            elif etype == '5':
                assert len(nodes) == 3, nodes
                tris.append(nodes)
            else:
                raise NotImplementedError(etype)
            i += 1

        tris = np.asarray(tris, dtype='int32')
        quads = np.asarray(quads, dtype='int32')

        i, points_dict = read_header(lines, i, self.log)
        nnodes = points_dict['npoints']
        nodes = np.zeros((nelem, 2), dtype='float32')
        for inode in range(nnodes):
            line = lines[i]
            sline = line.split()
            assert len(sline) == 3, sline
            x, y, unused_z = sline
            nodes[inode, :] = [x, y]
            i += 1

        # boundary conditions
        i, bcs = read_header(lines, i, self.log)
        nmark = bcs['nmark']
        for imark in range(nmark):
            nelements_mark = bcs['marker_elems']

            lines2d = []
            #np.zeros((nelements_mark, 2), dtype='int32')
            for unused_ne in range(nelements_mark):
                # what are the 3 slots?
                #Line          (2D)     3
                #Triangle      (3D)     5
                #Quadrilateral (3D)     9
                line = lines[i]
                sline = line.split()

                etype, *nodes = sline
                if etype == '3':
                    assert len(nodes) == 2, nodes
                    lines2d.append(nodes)
                else:
                    raise NotImplementedError(etype)
                i += 1

            lines2d = np.array(lines2d, dtype='int32')
            if imark != nmark - 1:
                i, bcs = read_header(lines, i, self.log)

        assert len(tris) > 0 or len(quads) > 0
        elements = {
            5 : tris,
            9 : quads,
        }
        regions = {3 : lines2d}
        return i, nodes, elements, regions

    def read_3d(self, lines, i, unused_ndim, nelem):
        """reads a 3d SU2 zone"""
        tets = []
        hexs = []
        wedges = []
        #pents = []
        pyramids = []
        for unused_ne in range(nelem):
            line = lines[i]
            etype, *nodes = line.split()[:-1]
            if etype == '10':
                tets.append(nodes)
            elif etype == '12':
                hexs.append(nodes)
            elif etype == '13':
                wedges.append(nodes)
            elif etype == '14':
                pyramids.append(nodes)
            else:
                raise NotImplementedError(f'etype={etype} data={data}')
            i += 1
        tets = np.array(tets, dtype='int32')
        #pents = np.array(pents, dtype='int32')
        wedges = np.array(wedges, dtype='int32')
        pyramids = np.array(pyramids, dtype='int32')
        elements = {
            10 : tets,
            12 : hexs,
            13 : wedges,
            14 : pyramids,
        }

        i, points_dict = read_header(lines, i, self.log)
        nnodes = points_dict['npoints']
        nodes = np.zeros((nelem, 3), dtype='float32')
        for inode in range(nnodes):
            line = lines[i]
            x, y, z = line.split()[:-1]
            nodes[inode, :] = [x, y, z]
            i += 1

        i, bcs = read_header(lines, i, self.log)
        nmark = bcs['nmark']
        for imark in range(nmark):
            nelements_mark = bcs['marker_elems']
            tris = []
            quads = []
            for unused_ne in range(nelements_mark):
                etype, *nodes = line.split()[:-1]
                if etype == '9':
                    assert nodes == 4, nodes
                    quads.append(nodes)
                elif etype == '5':
                    assert nodes == 3, nodes
                    tris.append(nodes)
                i += 1
            if imark != nmark - 1:
                i, bcs = read_header(lines, i, self.log)

            tris = np.array(tris, dtype='int32')
            quads = np.array(quads, dtype='int32')

        regions = {
            5 : tris,
            9 : quads,
        }
        return i, nodes, elements, regions

    def read_su2(self, su2_filename):
        """reads a 2d/3d SU2 single zone model"""
        self.su2_filename = su2_filename

        with open(su2_filename, 'r') as su2_file:
            lines = su2_file.readlines()

        i = 0
        i, zone = read_header(lines, i, self.log)
        nzones = zone['nzones']
        for izone in range(nzones):
            ndim = zone['ndim']
            ndim = zone['ndim']
            nelem = zone['nelem']

            if ndim == 2:
                i, nodes, elements, regions = self.read_2d(lines, i, ndim, nelem)
            elif ndim == 3:
                i, nodes, elements, regions = self.read_3d(lines, i, ndim, nelem)
            else:
                raise RuntimeError(ndim)
            break
        return ndim, nodes, elements, regions

    def write_su2(self, su2_filename, nodes, elements, unused_regions):
        nnodes, ndim = nodes.shape
        nnodes, ndim = nodes.shape
        with open(su2_filename, 'wb') as su2_file:
            su2_file.write('NDIM = %i\n' % ndim)
            if ndim == 2:
                for etype, elementsi in sorted(elements.items()):
                    element_num = self.etype_nnodes_map[etype]
                    fmt = '%%s' + ' %%s' * (element_num-1) + '\n'
                    for element in elementsi:
                        su2_file.write(fmt % element)

                su2_file.write('NPOINTS = %i\n' % nnodes)
                for inode, node in enumerate(nodes):
                    su2_file.write('%i %i %i\n' % (node[0], node[1], inode))
            elif ndim == 3:
                su2_file.write('NPOINTS = %i\n' % nnodes)
                for inode, node in enumerate(nodes):
                    su2_file.write('%i %i %i %i\n' % (node[0], node[1], node[2], inode))

def read_header(lines, i, log):
    """reads the header and stores it as a dictionary"""
    i0 = i
    headers = {}
    line = lines[i]

    while '=' in line:
        word, value = line.split('=')
        if word == 'NDIME':
            headers['ndim'] = int(value)
        elif word == 'NELEM':
            headers['nelem'] = int(value)
        elif word == 'NZONE':
            headers['nzones'] = int(value)
        elif word == 'IZONE':
            headers['izone'] = int(value)

        elif word == 'NMARK':
            headers['nmark'] = int(value)
        elif word == 'MARKER_TAG':
            headers['marker_tag'] = value.strip()
        elif word == 'MARKER_ELEMS':
            headers['marker_elems'] = int(value)

        elif word == 'NPOIN':
            headers['npoints'] = int(value)
        else:
            raise NotImplementedError('word=%r' % word)
        i += 1
        line = lines[i]

    assert len(headers) > 0, headers
    if i0 == 0 and 'nzones' not in headers and ('izone' in headers or'ndim' in headers):
        headers['nzones'] = 1

    line = lines[i]
    log.debug(str(headers))
    return i, headers
