import numpy as np


def read_su2(su2_filename, log=None, debug=False):
    """reads a 2d/3d SU2 multi-zone model"""
    model = SU2Reader(log=log, debug=debug)
    unused_ndim, zones = model.read_su2(su2_filename)
    #model.to_cart3d()
    return model, zones

class SU2Reader:
    """SU2 reader/writer"""
    etype_nnodes_map = {
        3 : 2,  # line
        5 : 3,  # triangle
        9 : 4,  # quad
        10 : 4,  # tetra
        12 : 8,  # hexa
        13 : 6,  # wedge / penta
        14 : 5,  # pyramid
    }
    def __init__(self, log=None, debug=False):
        """initializes the SU2Reader object"""
        self.log = log
        self.debug = debug
        self.su2_filename = None
        self.ndim = None
        self.zones = {}

    def _read_2d(self, lines, i, ndim, nelem):
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
            etype, *nodesi = data
            if etype == '9':
                assert len(nodesi) == 4, nodesi
                quads.append(nodesi)
            elif etype == '5':
                assert len(nodesi) == 3, nodesi
                tris.append(nodesi)
            else:
                raise NotImplementedError(etype)
            i += 1

        tris = np.asarray(tris, dtype='int32')
        quads = np.asarray(quads, dtype='int32')

        i, points_dict = read_header(lines, i, self.log)
        nnodes = points_dict['npoints']
        nodes = np.zeros((nnodes, 2), dtype='float32')
        for inode in range(nnodes):
            line = lines[i]
            sline = line.split()
            assert len(sline) == 3, sline
            x, y, unused_z = sline
            nodes[inode, :] = [x, y]
            i += 1

        # boundary conditions
        i, bcs = read_header(lines, i, self.log)
        i, regions = self._read_2d_bcs(lines, i, bcs)

        assert len(tris) > 0 or len(quads) > 0
        elements = {}
        if tris.shape[0] > 0:
            elements[5] = tris
        if quads.shape[0] > 0:
            elements[9] = quads
        assert len(elements), 'no 2d elements found'

        return i, nodes, elements, regions

    def _read_2d_bcs(self, lines, i, bcs):
        """reads the 1d boundary conditions for a 2d SU2 zone"""
        regions = {}
        nmark = bcs['nmark']
        line = lines[i]
        name_old = None
        for imark in range(nmark):
            nelements_mark = bcs['marker_elems']
            name = bcs['marker_tag']
            assert name != name_old, 'name={name!r} name_old={name_old!r}'

            lines2d = []
            #np.zeros((nelements_mark, 2), dtype='int32')
            for unused_ne in range(nelements_mark):
                # what are the 3 slots?
                #Line          (2D)     3
                #Triangle      (3D)     5
                #Quadrilateral (3D)     9
                line = lines[i]
                sline = line.split()

                etype, *nodesi = sline
                if etype == '3':
                    assert len(nodesi) == 2, f'nnodes={len(nodesi)} nodes={nodesi}'
                    lines2d.append(nodesi)
                else:  # pragma: no cover
                    raise NotImplementedError(etype)
                i += 1

            lines2d = np.array(lines2d, dtype='int32')
            if imark != nmark - 1:
                i, bcs = read_header(lines, i, self.log)
            regionsi = {3 : lines2d}
            assert len(lines2d), 'no 3d bcs found'
            regions[name] = regionsi
            name_old = name
        return i, regions

    def _read_3d(self, lines, i, unused_ndim, nelem):
        """reads a 3d SU2 zone"""
        tets = []
        hexs = []
        wedges = []
        #pents = []
        pyramids = []
        for unused_ne in range(nelem):
            line = lines[i]
            etype, *nodesi = line.split()[:-1]

            #Tetrahedral    10
            #Hexahedral     12
            #Wedge          13
            #Pyramid        14
            if etype == '10':  # tetra
                assert len(nodesi) == 4, f'nnodes={len(nodesi)} nodes={nodesi}'
                tets.append(nodesi)
            elif etype == '12':  # hexa
                assert len(nodesi) == 8, f'nnodes={len(nodesi)} nodes={nodesi}'
                hexs.append(nodesi)
            elif etype == '13':  # wedge
                assert len(nodesi) == 6, f'nnodes={len(nodesi)} nodes={nodesi}'
                wedges.append(nodesi)
            elif etype == '14':  # pyram
                assert len(nodesi) == 5, f'nnodes={len(nodesi)} nodes={nodesi}'
                pyramids.append(nodesi)
            else:  # pragma: no cover
                raise NotImplementedError(f'etype={etype} data={data}')
            i += 1
        tets = np.array(tets, dtype='int32')
        hexs = np.array(hexs, dtype='int32')
        wedges = np.array(wedges, dtype='int32')  # penta element
        pyramids = np.array(pyramids, dtype='int32')

        elements = {}
        if tets.shape[0] > 0:
            elements[10] = tets
        if hexs.shape[0] > 0:
            elements[12] = hexs
        if wedges.shape[0] > 0:
            elements[13] = wedges
        if pyramids.shape[0] > 0:
            elements[14] = pyramids
        assert len(elements), 'no 3d elements found'

        i, points_dict = read_header(lines, i, self.log)
        nnodes = points_dict['npoints']
        nodes = np.zeros((nnodes, 3), dtype='float32')
        for inode in range(nnodes):
            line = lines[i]
            x, y, z = line.split()[:-1]
            nodes[inode, :] = [x, y, z]
            i += 1

        i, bcs = read_header(lines, i, self.log)
        i, regions = self._read_3d_bcs(lines, i, bcs)

        return i, nodes, elements, regions


    def _read_3d_bcs(self, lines, i, bcs):
        """reads the 2d boundary conditions for a 3d SU2 zone"""
        regions = {}
        nmark = bcs['nmark']
        name_old = None
        for imark in range(nmark):
            #print(f'imark={imark} nmark={nmark} look_ahead={imark != nmark - 1}')
            nelements_mark = bcs['marker_elems']
            name = bcs['marker_tag']
            assert name != name_old, 'name=%r name_old=%r' % (name, name_old)
            tris = []
            quads = []
            for unused_ne in range(nelements_mark):
                line = lines[i]
                etype, *nodesi = line.split()
                if etype == '9':
                    assert len(nodesi) == 4, f'nnodes={len(nodesi)} nodes={nodesi}'
                    quads.append(nodesi)
                elif etype == '5':
                    assert len(nodesi) == 3, f'nnodes={len(nodesi)} nodes={nodesi}'
                    tris.append(nodesi)
                else:  # pragma: no cover
                    raise RuntimeError(line.strip())
                i += 1
            if imark != nmark - 1:
                i, bcs = read_header(lines, i, self.log)

            tris = np.array(tris, dtype='int32')
            quads = np.array(quads, dtype='int32')

            regionsi = {}
            if tris.shape[0] > 0:
                regionsi[5] = tris
            if quads.shape[0] > 0:
                regionsi[9] = quads
            assert len(regionsi), 'no 3d bcs found'
            regions[name] = regionsi
            name_old = name
        return i, regions

    def read_su2(self, su2_filename):
        """reads a 2d/3d SU2 single zone model"""
        self.su2_filename = su2_filename

        with open(su2_filename, 'r') as su2_file:
            lines = su2_file.readlines()

        i = 0
        zones = {}
        i, zone_dict = read_header(lines, i, self.log)
        nzones = zone_dict['nzones']
        for izone in range(nzones):
            ndim = zone_dict['ndim']
            nelem = zone_dict['nelem']

            if ndim == 2:
                i, nodes, elements, regions = self._read_2d(lines, i, ndim, nelem)
            elif ndim == 3:
                i, nodes, elements, regions = self._read_3d(lines, i, ndim, nelem)
            else:  # pragma: no cover
                raise RuntimeError(ndim)

            if izone != nzones - 1:
                i, zone_dict = read_header(lines, i, self.log)
            zones[izone] = (nodes, elements, regions)

        self.ndim = ndim
        self.zones = zones
        return ndim, zones

    def write_su2(self, su2_filename, zones=None):
        """writes a 2d/3d multizone SU2 file"""
        if zones is None:
            zones = self.zones

        nzones = len(zones)
        with open(su2_filename, 'w') as su2_file:
            su2_file.write(f'NZONE = {nzones:d}\n')
            for izone, zone in zones.items():
                nodes, elements, regions = zone
                nnodes, ndim = nodes.shape
                su2_file.write(f'IZONE = {nzones:d}\n')
                su2_file.write(f'NDIM = {ndim:d}\n')

                nelements = 0
                for etype, elementsi in elements.items():
                    nelements += elementsi.shape[0]
                su2_file.write(f'NELEM = {nelements:d}\n')

                for etype, elementsi in elements.items():
                    nnodesi = self.etype_nnodes_map[etype]
                    fmt = str(etype) + ' %s' * nnodesi + '\n'
                    for element in elementsi:
                        su2_file.write(fmt % tuple(element))

                su2_file.write('NPOINTS = %i\n' % nnodes)
                if ndim == 2:
                    for inode, node in enumerate(nodes):
                        su2_file.write('%f %f %i\n' % (node[0], node[1], inode))
                elif ndim == 3:
                    su2_file.write('NPOINTS = %i\n' % nnodes)
                    for inode, node in enumerate(nodes):
                        su2_file.write('%f %f %f %i\n' % (node[0], node[1], node[2], inode))
                else:  # pragma: no cover
                    raise RuntimeError(ndim)

                nregions = len(regions)
                su2_file.write(f'NMARK = {nregions:d}\n')
                for name, region in regions.items():
                    nmarker_elementsi = 0
                    for etype, elementsi in region.items():
                        nmarker_elementsi += elementsi.shape[0]

                    su2_file.write(f'MARKER_ELEMS = {nmarker_elementsi:d}\n')
                    su2_file.write(f'MARKER_TAG = {name}\n')
                    for etype, elementsi in region.items():
                        nnodesi = self.etype_nnodes_map[etype]
                        fmt = str(etype) + ' %s' * (nnodesi) + '\n'
                        for element in elementsi:
                            su2_file.write(fmt % tuple(element))


def read_header(lines, i, log):
    """reads the header and stores it as a dictionary"""
    i0 = i
    headers = {}
    line = lines[i]
    while '=' not in line:
        line = line.split('%')[0].strip()
        if line == '':
            i += 1
            line = lines[i]
            continue
        #line_upper = line.upper()
        #if not line.startswith(('N', 'I', 'M')):
        raise RuntimeError('expected a header; found=%r' % line)

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
        else:  # pragma: no cover
            raise NotImplementedError('word=%r' % word)
        i += 1
        line = lines[i]

    assert len(headers) > 0, headers

    # if this is the first zone and nzones isn't defined, add it
    if i0 == 0 and 'nzones' not in headers and ('izone' in headers or'ndim' in headers):
        headers['nzones'] = 1

    line = lines[i]
    log.debug(str(headers))
    return i, headers
