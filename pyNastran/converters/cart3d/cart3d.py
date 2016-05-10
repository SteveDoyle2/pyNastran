from __future__ import print_function, unicode_literals
from six import iteritems, b, PY2
from six.moves import zip, range
import sys
from struct import Struct, pack, unpack
from math import ceil
from collections import defaultdict
from codecs import open as codec_open

from numpy import zeros, where, savetxt, sqrt, amax, amin
from numpy import arange, ravel, cross
from numpy.linalg import norm

from pyNastran.utils import is_binary_file, _filename
from pyNastran.utils.log import get_logger


def comp2tri(self, in_filenames, out_filename,
             is_binary=False, float_fmt='%6.7f'):
    """
    Combines multiple Cart3d files (binary or ascii) into a single file.

    Parameters
    ----------
    in_filenames : List[str]
        list of filenames
    out_filename : str
        output filename
    is_binary : bool; default=False
        is the output filename binary
    float_fmt : str; default='%6.7f'
        the format string to use for ascii writing

    .. note:: assumes loads is None
    """
    points = []
    elements = []
    regions = []

    #ne = 0
    npi = 0
    nri = 0
    model = Cart3D()
    for infilename in in_filenames:
        model.read_cart3d(infilename)
        np = point.shape[0]
        nr = len(unique(region))
        element += npi - 1
        region += nri

        points.append(model.nodes)
        elements.append(model.elements)
        regions.append(model.regions)
        npi += np
        nri += nr

    points = vstack(points)
    elements = vstack(elements)
    regions = vstack(regions)
    model.write_cart3d(out_filename,
                      points, elements, regions,
                      loads=None, is_binary=False, float_fmt=float_fmt)


class Cart3dIO(object):
    """
    Cart3d IO class
    """
    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')
        self._endian = b''
        self._encoding = 'latin1'
        self.n = 0
        self.infile = None
        # self.readHalf = False
        # self.nPoints = None
        # self.nElements = None
        self.infilename = None
        self.points = None
        self.elements = None
        self.regions = None
        self.loads = {}

    def _write_header(self, outfile, points, elements, is_loads, is_binary=False):
        npoints = points.shape[0]
        nelements = elements.shape[0]

        if is_binary:
            if is_loads:
                msg = pack('>iiiii', 3*4, npoints, nelements, 6, 4)
            else:
                msg = pack('>iiii', 2*4, npoints, nelements, 4)

            int_fmt = None
        else:
            if is_loads:
                msg = b"%i %i 6\n" % (npoints, nelements)
            else:
                msg = b"%i %i\n" % (npoints, nelements)

            # take the max value, string it, and length it
            # so 123,456 is length 6
            int_fmt = '%%%si' % len(str(nelements))
        outfile.write(msg)
        return int_fmt

    def _write_points(self, outfile, points, is_binary, float_fmt='%6.6f'):
        if is_binary:
            four = pack(b'>i', 4)
            outfile.write(four)

            npoints = points.shape[0]
            fmt = b'>%if' % (npoints * 3)
            floats = pack(fmt, *ravel(points))

            outfile.write(floats)
            outfile.write(four)
        else:
            if isinstance(float_fmt, str):
                fmt = float_fmt
            else:
                fmt = float_fmt.encode('latin1')
            savetxt(outfile, points, fmt)

    def _write_elements(self, outfile, elements, is_binary, int_fmt='%6i'):
        min_e = elements.min()
        assert min_e == 1, 'min(elements)=%s' % min_e
        if is_binary:
            four = pack(b'>i', 4)
            outfile.write(four)
            nelements = elements.shape[0]
            fmt = b'>%ii' % (nelements * 3)
            ints = pack(fmt, *ravel(elements))

            outfile.write(ints)
            outfile.write(four)
        else:
            if isinstance(int_fmt, str):
                fmt = int_fmt
            else:
                fmt = int_fmt.encode('latin1')
            savetxt(outfile, elements, fmt)

    def _write_regions(self, outfile, regions, is_binary):
        if is_binary:
            four = pack('>i', 4)
            outfile.write(four)

            nregions = len(regions)
            fmt = b'>%ii' % nregions
            ints = pack(fmt, *regions)
            outfile.write(ints)

            outfile.write(four)
        else:
            savetxt(outfile, regions, b'%i')

    def _write_loads(self, outfile, loads, is_binary, float_fmt='%6.6f'):
        if is_binary:
            raise NotImplementedError('is_binary=%s' % is_binary)
        else:
            Cp = loads['Cp']
            rho = loads['rho']
            rhoU = loads['rhoU']
            rhoV = loads['rhoV']
            rhoW = loads['rhoW']
            E = loads['E']

            nrows = len(Cp)
            fmt = '%s\n%s %s %s %s %s\n' % (float_fmt, float_fmt, float_fmt, float_fmt, float_fmt, float_fmt)
            for (cpi, rhoi, rhou, rhov, rhoe, e) in zip(Cp, rho, rhoU, rhoV, rhoW, E):
                outfile.write(fmt % (cpi, rhoi, rhou, rhov, rhoe, e))

    def _read_header_ascii(self):
        line = self.infile.readline()
        sline = line.strip().split()
        if len(sline) == 2:
            npoints, nelements = int(sline[0]), int(sline[1])
            nresults = 0
        elif len(sline) == 3:
            npoints = int(sline[0])
            nelements = int(sline[1])
            nresults = int(sline[2])
        else:
            raise ValueError('invalid result type')
        return npoints, nelements, nresults

    @property
    def nresults(self):
        if isinstance(self.loads, dict):
            return len(self.loads)
        return 0

    @property
    def nnodes(self):
        return self.npoints

    @property
    def npoints(self):
        return self.points.shape[0]

    @property
    def nodes(self):
        return self.points

    @nodes.setter
    def nodes(self, points):
        self.points = points

    @property
    def nelements(self):
        return self.elements.shape[0]

    def _read_points_ascii(self, npoints):
        """
        A point is defined by x,y,z and the ID is the location in points.
        """
        p = 0
        data = []
        assert npoints > 0, 'npoints=%s' % npoints
        points = zeros((npoints, 3), dtype='float32')
        while p < npoints:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                x = data.pop(0)
                y = data.pop(0)
                z = data.pop(0)
                points[p] = [x, y, z]
                p += 1

        #maxX = self.get_max(points, 0)
        #maxY = self.get_max(points, 1)
        #maxZ = self.get_max(points, 2)

        #minX = self.get_min(points, 0)
        #minY = self.get_min(points, 1)
        #minZ = self.get_min(points, 2)

        #self.log.debug("X  max=%g min=%g" % (maxX, minX))
        #self.log.debug("Y  max=%g min=%g" % (maxY, minY))
        #self.log.debug("Z  max=%g min=%g" % (maxZ, minZ))
        return points

    def _read_elements_ascii(self, nelements):
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.
        """
        assert nelements > 0, 'npoints=%s nelements=%s' % (self.npoints, nelements)
        elements = zeros((nelements, 3), dtype='int32')

        e = 0
        data = []
        while e < nelements:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                n1 = int(data.pop(0))
                n2 = int(data.pop(0))
                n3 = int(data.pop(0))
                elements[e] = [n1, n2, n3]
                e += 1
        return elements

    def _read_regions_ascii(self, nelements):
        regions = zeros(nelements, dtype='int32')
        r = 0
        data = []
        while r < nelements:
            data = self.infile.readline().strip().split()
            ndata = len(data)
            regions[r : r + ndata] = data
            r += ndata
        return regions

    def _read_header_binary(self):
        data = self.infile.read(4)
        size_little, = unpack(b'<i', data)
        size_big, = unpack(b'>i', data)
        if size_big in [12, 8]:
            self._endian = b'>'
            size = size_big
        elif size_little in [8, 12]:
            self._endian = b'<'
            size = size_little
        else:
            self.rewind()
            self.show(100)
            raise RuntimeError('unknown endian')

        self.n += 4
        data = self.infile.read(size)
        self.n += size

        so4 = size // 4  # size over 4
        if so4 == 3:
            (npoints, nelements, nresults) = unpack(self._endian + b'iii', data)
            self.log.info("npoints=%s nelements=%s nresults=%s" % (npoints, nelements, nresults))
        elif so4 == 2:
            (npoints, nelements) = unpack(self._endian + b'ii', data)
            nresults = 0
            self.log.info("npoints=%s nelements=%s" % (npoints, nelements))
        else:
            self.rewind()
            self.show(100)
            raise RuntimeError('in the wrong spot...endian...size/4=%s' % so4)
        self.infile.read(8)  # end of first block, start of second block
        return (npoints, nelements, nresults)

    def _read_points_binary(self, npoints):
        size = npoints * 12  # 12=3*4 all the points

        n = 0
        points = zeros(npoints * 3, dtype='float32')
        s = Struct(self._endian + b'3000f') # 3000 floats; 1000 points
        while size > 12000:  # 12k = 4 bytes/float*3 floats/point*1000 points
            data = self.infile.read(4 * 3000)

            node_xyzs = s.unpack(data)
            points[n:n+3000] = node_xyzs
            n += 3000
            size -= 4 * 3000

        assert size >= 0, 'size=%s' % size

        if size > 0:
            data = self.infile.read(size)
            bin_format = self._endian + b'%if' % (size // 4)

            node_xyzs = unpack(bin_format, data)
            points[n:] = node_xyzs

        points = points.reshape((npoints, 3))

        self.infile.read(8)  # end of second block, start of third block
        return points

    def _read_elements_binary(self, nelements):
        size = nelements * 12  # 12=3*4 all the elements

        elements = zeros(nelements * 3, dtype='int32')

        n = 0
        s = Struct(self._endian + b'3000i')
        while size > 12000:  # 4k is 1000 elements
            data = self.infile.read(4 * 3000)
            nodes = s.unpack(data)
            elements[n : n + 3000] = nodes
            size -= 4 * 3000
            n += 3000

        assert size >= 0, 'size=%s' % size
        if size > 0:
            data = self.infile.read(size)
            bin_format = self._endian + b'%ii' % (size // 4)

            nodes = unpack(bin_format, data)
            elements[n:] = nodes

        elements2 = elements.reshape((nelements, 3))
        self.infile.read(8)  # end of third (element) block, start of regions (fourth) block
        return elements2

    def _read_regions_binary(self, nelements):
        size = nelements * 4  # 12=3*4 all the elements
        s = Struct(self._endian + b'3000i')

        regions = zeros(nelements, dtype='int32')

        nr = 0
        while size > 12000:  # 12k is 3000 elements
            data = self.infile.read(4 * 3000)
            try:
                region_data = s.unpack(data)
            except:
                print("len =", len(data))
                raise

            delta = len(region_data)
            regions[nr:nr+delta] = region_data
            nr += delta
            size -= 4 * 3000

        assert size >= 0, 'size=%s' % size
        if size > 0:
            data = self.infile.read(size)
            bin_format = self._endian + b'%ii' % (size // 4)
            try:
                region_data = unpack(bin_format, data)
            except:
                print("len =", len(data))
                raise

            r = nelements
            if len(regions[nr:]) != len(region_data):
                msg = 'len(regions[nr:]=%s len(region_data)=%s' % (len(regions[nr:]), len(region_data))
                raise RuntimeError(msg)
            regions[nr:] = region_data
            size = 0

        self.infile.read(4)  # end of regions (fourth) block
        return regions

    def _read_results_binary(self, i, infile, result_names=None):
        pass

    def rewind(self):
        self.n = 0
        self.infile.seek(self.n)

    def show(self, n, types='ifs', endian=None):
        assert self.n == self.infile.tell(), 'n=%s tell=%s' % (self.n, self.infile.tell())
        nints = n // 4
        data = self.infile.read(4 * n)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        self.infile.seek(self.n)
        return strings, ints, floats

    def show_data(self, data, types='ifs', endian=None):
        return self.write_data(sys.stdout, data, types=types, endian=endian)

    def write_data(self, outfile, data, types='ifs', endian=None):
        """
        Useful function for seeing what's going on locally when debugging.
        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._endian

        if 's' in types:
            strings = unpack(b'%s%is' % (endian, n), data)
            outfile.write("strings = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack(b'%s%ii' % (endian, nints), data)
            outfile.write("ints    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack(b'%s%if' % (endian, nints), data)
            outfile.write("floats  = %s\n" % str(floats))

        if 'l' in types:
            longs = unpack(b'%s%il' % (endian, nints), data)
            outfile.write("long  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack(b'%s%iI' % (endian, nints), data)
            outfile.write("unsigned int = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack(b'%s%iL' % (endian, nints), data)
            outfile.write("unsigned long = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack(b'%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            outfile.write("long long = %s\n" % str(longs))
        return strings, ints, floats

    def show_ndata(self, n, types='ifs'):
        return self.write_ndata(sys.stdout, n, types=types)

    def write_ndata(self, outfile, n, types='ifs'):
        """
        Useful function for seeing what's going on locally when debugging.
        """
        nold = self.n
        data = self.infile.read(n)
        self.n = nold
        self.infile.seek(self.n)
        return self.write_data(outfile, data, types=types)


class Cart3D(Cart3dIO):
    """
    Cart3d interface class
    """
    model_type = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        Cart3dIO.__init__(self, log=log, debug=debug)
        self.loads = {}
        self.points = None
        self.elements = None

    def make_mirror_model(self, nodes, elements, regions, loads, axis='y', tol=0.000001):
        """
        Makes a full cart3d model about the given axis.

        Parameters
        ----------
        nodes : (nnodes, 3) ndarray
            the nodes
        elements : (nelements, 3) ndarray
            the elmements
        regions :  (nelements) ndarray
            the regions
        loads : dict[str] = (nnodes) ndarray
            not supported
        axis : str; {"x", "y", "z", "-x", "-y", "-z"}
            a string of the axis
        tol : float; default=0.000001
            the tolerance for the centerline points
        """
        raise NotImplementedError()
        self.log.info('---starting make_mirror_model---')
        assert tol >= 0, 'tol=%r' % tol #  prevents hacks to the axis

        nnodes = nodes.shape[0]
        assert nnodes > 0, 'nnodes=%s' % nnodes

        nelements = elements.shape[0]
        assert nelements > 0, 'nelements=%s' % nelements

        ax = self._get_ax(axis)
        if ax in [0, 1, 2]:  # positive x, y, z values; mirror to -side
            iy0 = where(nodes[:, ax] > tol)[0]
            ax2 = ax
        elif ax in [3, 4, 5]:  # negative x, y, z values; mirror to +side
            iy0 = where(nodes[:, ax-3] < -tol)[0]
            ax2 = ax - 3 # we create ax2 in order to generalize the data updating
        else:
            raise NotImplementedError(axis)

        # the nodes to be duplicated are the nodes that aren't below the tolerance
        nodes_upper = nodes[iy0]
        nodes_upper[:, ax2] *= -1.0  # flip the nodes about the axis

        nodes2 = vstack([nodes, nodes_upper])
        nnodes2 = nodes2.shape[0]
        assert nnodes2 > nnodes, 'nnodes2=%s nnodes=%s' % (nnodes2, nnodes)

        nnodes_upper = nodes_upper.shape[0]
        elements_upper = elements.copy()
        nelements = elements.shape[0]

        # remap the mirrored nodes with the new node ids
        for eid in range(nelements):
            element = elements_upper[eid, :]
            for i, eidi in enumerate(element):
                if eidi in iy0:
                    elements_upper[eid][i] = nnodes_upper + eidi

                # we need to reverse the element in order to get
                # the proper normal vector
                elements_upper[eid] = elements_upper[eid, ::-1]

        elements2 = vstack([elements, elements_upper])
        nelements2 = elements2.shape[0]
        assert nelements2 > nelements, 'nelements2=%s nelements=%s' % (nelements2, nelements)

        nregions = len(unique(regions))
        regions_upper = regions.copy() + nregions
        regions2 = hstack([regions, regions_upper])

        loads2 = {}
        for key, data in iteritems(loads):

            # flip the sign on the flipping axis terms
            if((key in ['U', 'rhoU'] and ax2 == 0) or
               (key in ['V', 'rhoV'] and ax2 == 1) or
               (key in ['W', 'rhoW'] and ax2 == 2)):
                data_upper = -data[iy0]
            else:
                data_upper = data[iy0]
            loads2[key] = hstack([data, data_upper])

        self.log.info('---finished make_mirror_model---')
        return (nodes2, elements2, regions2, loads2)

    def _get_ax(self, axis):
        axis = axis.lower().strip()
        if axis in ['+x', 'x', 0]:
            ax = 0
        elif axis in ['+y', 'y', 1]:
            ax = 1
        elif axis in ['+z', 'z', 2]:
            ax = 2

        elif axis in ['-x', 3]:
            ax = 3
        elif axis == ['-y', 4]:
            ax = 4
        elif axis == ['-z', 5]:
            ax = 5
        else:
            raise NotImplementedError('axis=%r' % axis)
        self.log.info("axis=%r ax=%s" % (axis, ax))
        return ax

    def make_half_model(self, axis='y', remap_nodes=True):
        """
        Makes a half model from a full model

        ... note:: Cp is really loads['Cp'] and was meant for loads analysis only
        """
        nodes = self.nodes
        elements = self.elements
        loads = self.loads
        if loads is None:
            loads = {}

        nnodes = nodes.shape[0]
        assert nnodes > 0, 'nnodes=%s'  % nnodes

        nelements = elements.shape[0]
        assert nelements > 0, 'nelements=%s'  % nelements

        inodes_remove = set([])
        self.log.info('---starting make_half_model---')
        ax = self._get_ax(axis)

        if ax in [0, 1, 2]:  # remove values > 0
            inodes_save = where(nodes[:, ax] >= 0.0)[0]
        elif ax in [3, 4, 5]:  # remove values < 0
            inodes_save = where(nodes[:, ax-3] <= 0.0)[0]
        else:
            raise NotImplementedError('axis=%r ax=%s' % (axis, ax))
        inodes_save.sort()

        inodes_map = arange(len(inodes_save))
        assert 0 < len(inodes_save) < nnodes, 'len(inodes_save)=%s nnodes=%s'  % (len(inodes_save), nnodes)

        nodes2 = nodes[inodes_save, :]
        nnodes2 = nodes2.shape[0]
        assert 0 < nnodes2 < nnodes, 'nnodes=%s nnodes2=%s'  % (nnodes, nnodes2)

        inodes_save += 1  # +1 is so we don't have to shift inode
        # .. todo:: still need to handle element's node id renumbering
        ielements_save = set([])
        for ielement in range(nelements):
            save_element = True
            element = elements[ielement, :]

            # could be faster...
            for inode in element:
                if inode not in inodes_save:
                    save_element = False
                    break

            if save_element:
                ielements_save.add(ielement)

        ielements_save = list(ielements_save)
        ielements_save.sort()

        elements2 = elements[ielements_save]
        regions2 = regions[ielements_save]

        # renumbers mesh
        nelements2 = elements2.shape[0]
        assert 0 < nelements2 < nelements, 'nelements=%s nelements2=%s'  % (nelements, nelements2)

        remap_nodes = False
        if amax(elements2) > len(inodes_save):
            # build a dictionary of old node ids to new node ids
            nodes_map = {}
            for i in range(1, len(inodes_save) + 1):
                nid = inodes_save[i - 1]
                nodes_map[nid] = i

            # update the node ids
            for ielement in range(nelements2):
                element = elements2[ielement, :]
                elements[ielement, :] = [nodes_map[nid] for nid in element]

        loads2 = {} # 'Cp', 'Mach', 'U', etc.
        for key, load in iteritems(loads):
            loads2[key] = load[inodes_save]

        self.log.info('---finished make_half_model---')
        return (nodes2, elements2, regions2, loads2)

    def get_free_edges(self, elements):
        edge_to_eid_map = defaultdict(list)
        for i, element in enumerate(elements):
            edge1 = tuple(sorted([element[0], element[1]]))
            edge2 = tuple(sorted([element[1], element[2]]))
            edge3 = tuple(sorted([element[2], element[0]]))
            edge_to_eid_map[edge1].append(i)
            edge_to_eid_map[edge2].append(i)
            edge_to_eid_map[edge3].append(i)

        free_edges = []
        for edge, eids in sorted(iteritems(edge_to_eid_map)):
            if len(eids) != 2:
                free_edges.append(edge)
                # print(edge)
        return free_edges

    def read_cart3d(self, infilename, result_names=None):
        """extracts the points, elements, and Cp"""
        self.infilename = infilename
        self.log.info("---starting reading cart3d file...%r---" % self.infilename)

        self.infilename = infilename
        if is_binary_file(infilename):
            with open(infilename, 'rb') as self.infile:
                npoints, nelements, nresults = self._read_header_binary()
                self.points = self._read_points_binary(npoints)
                self.elements = self._read_elements_binary(nelements)
                self.regions = self._read_regions_binary(nelements)
                # TODO: loads
        else:
            with codec_open(_filename(infilename), 'r', encoding=self._encoding) as self.infile:
                npoints, nelements, nresults = self._read_header_ascii()
                self.points = self._read_points_ascii(npoints)
                self.elements = self._read_elements_ascii(nelements)
                self.regions = self._read_regions_ascii(nelements)
                self._read_results_ascii(0, self.infile, nresults, result_names=result_names)

        self.log.debug("npoints=%s nelements=%s" % (self.npoints, self.nelements))
        self.log.info("---finished reading cart3d file...%r---" % self.infilename)
        assert self.npoints > 0, 'npoints=%s' % self.npoints
        assert self.nelements > 0, 'nelements=%s' % self.nelements

    def write_cart3d(self, outfilename, is_binary=False, float_fmt='%6.7f'):
        assert len(self.points) > 0, 'len(self.points)=%s' % len(self.points)

        if self.loads is None or self.loads == {}:
            loads = {}
            is_loads = False
            # print("no loads")
        else:
            is_loads = True

        self.log.info("---writing cart3d file...%r---" % outfilename)
        if is_binary:
            form = 'wb'
        else:
            if PY2:
                form = 'w'
            else:
                form = 'wb'
        with codec_open(outfilename, form) as outfile:
            int_fmt = self._write_header(outfile, self.points, self.elements, is_loads, is_binary)
            self._write_points(outfile, self.points, is_binary, float_fmt)
            self._write_elements(outfile, self.elements, is_binary, int_fmt)
            self._write_regions(outfile, self.regions, is_binary)

            if is_loads:
                assert is_binary is False, 'is_binary=%r is not supported for loads' % is_binary
                self._write_loads(outfile, self.loads, is_binary, float_fmt)
        outfile.close()


    def get_min(self, points, i):
        return amin(points[:, i])

    def get_max(self, points, i):
        return amax(points[:, i])

    def _read_results_ascii(self, i, infile, nresults, result_names=None):
        """
        Reads the Cp results.
        Results are read on a nodal basis from the following table:
          Cp
          rho,rhoU,rhoV,rhoW,rhoE

        With the following definitions:
          Cp = (p - 1/gamma) / (0.5*M_inf*M_inf)
          rhoVel^2 = rhoU^2+rhoV^2+rhoW^2
          M^2 = rhoVel^2/rho^2

        Thus:
          p = (gamma-1)*(e- (rhoU**2+rhoV**2+rhoW**2)/(2.*rho))
          p_dimensional = qInf * Cp + pInf

        # ???
        rho,rhoU,rhoV,rhoW,rhoE

        Parameters
        ----------
        result_names : List[str]; default=None (All)
            result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                            'Mach', 'U', 'V', 'W', 'E']
        """
        if nresults == 0:
            return
        loads = {}
        if result_names is None:
            result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                            'Mach', 'U', 'V', 'W', 'E', 'a', 'T', 'Pressure', 'q']
        self.log.info('---starting read_results---')

        results = zeros((self.npoints, 6), dtype='float32')

        nresult_lines = int(ceil(nresults / 5.)) - 1
        for ipoint in range(self.npoints):
            # rho rhoU,rhoV,rhoW,pressure/rhoE/E
            sline = infile.readline().strip().split()
            i += 1
            for n in range(nresult_lines):
                sline += infile.readline().strip().split()  # Cp
                i += 1
                #gamma = 1.4
                #else:
                #    p=0.
            sline = _get_list(sline)

            # Cp
            # rho       rhoU      rhoV      rhoW      E
            # 0.416594
            # 1.095611  0.435676  0.003920  0.011579  0.856058
            results[ipoint, :] = sline

            #p=0
            #cp = sline[0]
            #rho = float(sline[1])
            #if(rho > abs(0.000001)):
                #rhoU = float(sline[2])
                #rhoV = float(sline[3])
                #rhoW = float(sline[4])
                #rhoE = float(sline[5])
                #mach2 = (rhoU) ** 2 + (rhoV) ** 2 + (rhoW) ** 2 / rho ** 2
                #mach = sqrt(mach2)
                #if mach > 10:
                    #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (pointNum, cp, mach, rho, rhoU, rhoV, rhoW))
            #print("pt=%s i=%s Cp=%s p=%s" %(pointNum,i,sline[0],p))
        del sline
        self.loads = self._calculate_results(result_names, results)

    def _calculate_results(self, result_names, results, loads=None):
        """
        Takes the Cart3d variables and calculates additional variables

        Parameters
        ----------
        result_names : List[str]
            the variables to calculate
        results : (n,6) ndarray
            the non-dimensional prmitive flow variables
        """
        if loads is None:
            loads = {}
        Cp = results[:, 0]
        rho = results[:, 1]
        rhoU = results[:, 2]
        rhoV = results[:, 3]
        rhoW = results[:, 4]
        E = results[:, 5]

        ibad = where(rho <= 0.000001)[0]
        if len(ibad) > 0:

            if 'Mach' in result_names:
                Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2)# / rho
                Mach[ibad] = 0.0
            if 'U' in result_names:
                U = rhoU / rho
                U[ibad] = 0.0
            if 'U' in result_names:
                V = rhoV / rho
                V[ibad] = 0.0
            if 'W' in result_names:
                W = rhoW / rho
                W[ibad] = 0.0
            #if 'rhoE' in result_names:
                #rhoE = rhoE / rho
                #e[ibad] = 0.0

            is_bad = True
            n = 0
            #for i in ibad:
                #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (i, Cp[i], Mach[i], rho[i],
                #                                                                  rhoU[i], rhoV[i], rhoW[i]))
                #Mach[i] = 0.0
                #n += 1
                #if n > 10:
                #    break
        else:
            is_bad = False


        #loc = locals()
        if 'Cp' in result_names:
            loads['Cp'] = Cp
        if 'rhoU' in result_names:
            loads['rhoU'] = rhoU
        if 'rhoV' in result_names:
            loads['rhoV'] = rhoV
        if 'rhoW' in result_names:
            loads['rhoW'] = rhoW
        #if 'rhoE' in result_names:
            #loads['rhoE'] = rhoE

        if 'rho' in result_names:
            loads['rho'] = rho

        if 'Mach' in result_names:
            if not is_bad:
                #Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2) / rho
                Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2)
            loads['Mach'] = Mach

        if 'U' in result_names:
            if not is_bad:
                U = rhoU / rho
            loads['U'] = U
        if 'V' in result_names:
            if not is_bad:
                V = rhoV / rho
            loads['V'] = V
        if 'W' in result_names:
            if not is_bad:
                W = rhoW / rho
            loads['W'] = W
        if 'E' in result_names:
            #if not is_bad:
                #E = rhoE / rho
            loads['E'] = E

        gamma = 1.4
        qinf = 1.0
        pinf = 1. / gamma
        Tinf = 1.0
        #Cp = (p - pinf) / qinf
        p = Cp * qinf + pinf

        T = (Tinf * gamma) * p / rho
        q = 0.5 * rho * Mach ** 2

        if 'a' in result_names:
            loads['a'] = sqrt(T)
        if 'T' in result_names:
            loads['T'] = T

        if 'Pressure' in result_names:
            loads['Pressure'] = p
        if 'q' in result_names:
            loads['q'] = q
        # dynamic pressure
        # speed of sound
        # total pressure = p0/rhoi*ainf**2
        # total density
        # entropy
        # kinetic energy
        # enthalpy
        # energy, E
        # total energy
        # total enthalpy

        #i = where(Mach == max(Mach))[0][0]
        #self.log.info("i=%s Cp=%s rho=%s rhoU=%s rhoV=%s rhoW=%s Mach=%s" % (i, Cp[i], rho[i], rhoU[i],
        #              rhoV[i], rhoW[i], Mach[i]))
        self.log.info('---finished read_results---')
        return loads

    def get_normals(self, shift_nodes=True):
        """
        Gets the centroidal normals

        Parameters
        ----------
        shift_nodes : boolean; default=True
            shifts element IDs such that the
              - node IDs start at 0 instead of 1
                  True : nodes start at 1
                  False : nodes start at 0

        Returns
        -------
        cnormals : (n, 3) ndarray
            normalized centroidal normal vectors
        """
        elements = self.elements
        nodes = self.nodes
        if shift_nodes:
            p1 = nodes[elements[:, 0] - 1, :]
            p2 = nodes[elements[:, 1] - 1, :]
            p3 = nodes[elements[:, 2] - 1, :]
        else:
            p1 = nodes[elements[:, 0], :]
            p2 = nodes[elements[:, 1], :]
            p3 = nodes[elements[:, 2], :]

        ne = elements.shape[0]
        a = p2 - p1
        b = p3 - p1
        n = cross(a, b)
        assert len(n) == ne, 'len(n)=%s ne=%s' % (len(n), ne)

        ni = norm(n, axis=1)
        assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)

        assert ni.min() > 0.0, ni[where(ni <= 0.0)[0]]
        n /= ni[:, None]  # normal vector
        return n

    def get_normals_at_nodes(self, cnormals, shift_nodes=True):
        """
        Gets the nodal normals

        Parameters
        ----------
        cnormals : (n, 3) ndarray
            normalized centroidal normal vectors
        shift_nodes : bool; default=True
            shifts element IDs such that the node IDs start at 0
            instead of 1
              True : nodes start at 1
              False : nodes start at 0

        Returns
        -------
        nnormals : (n, 3) ndarray
            normalized nodal normal vectors
        """
        elements = self.elements
        nodes = self.nodes
        nnodes = self.nnodes
        nid_to_eids = defaultdict(list)

        if shift_nodes:
            for eid, element in enumerate(elements - 1):
                n1, n2, n3 = element
                nid_to_eids[n1].append(eid)
                nid_to_eids[n2].append(eid)
                nid_to_eids[n3].append(eid)
        else:
            # find the elements to consider for each node
            for eid, element in enumerate(elements):
                n1, n2, n3 = element
                nid_to_eids[n1].append(eid)
                nid_to_eids[n2].append(eid)
                nid_to_eids[n3].append(eid)

        nnormals = zeros((nnodes, 3), dtype='float64')
        for nid in range(nnodes):
            eids = nid_to_eids[nid]
            if len(eids) == 0:
                raise RuntimeError('nid=%s is not used' % nid)
            ni_avg = cnormals[eids, :]
            nnormals[nid] = cnormals[eids, :].sum(axis=0)
        ni = norm(nnormals, axis=1)
        assert ni.min() > 0, ni
        nnormals /= ni[:, None]  # normal vector
        return nnormals

def convert_to_float(svalues):
    """Takes a list of strings and converts them to floats."""
    values = []
    for value in svalues:
        values.append(float(value))
    return values

def _get_list(sline):
    """Takes a list of strings and converts them to floats."""
    try:
        sline2 = convert_to_float(sline)
    except ValueError:
        print("sline = %s" % sline)
        raise SyntaxError('cannot parse %s' % sline)
    return sline2
