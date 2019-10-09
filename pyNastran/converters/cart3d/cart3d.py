"""
Defines:
  - Cart3D(log=None, debug=False)
     - read_cart3d(self, infilename, result_names=None)
     - write_cart3d(self, outfilename, is_binary=False, float_fmt='%6.7f')

     - flip_model()
     - make_mirror_model(self, nodes, elements, regions, loads, axis='y', tol=0.000001)
     - make_half_model(self, axis='y', remap_nodes=True)
     - get_free_edges(self, elements)
     - get_area(self)
     - get_normals(self)
     - get_normals_at_nodes(self, cnormals)

  - comp2tri(in_filenames, out_filename,
             is_binary=False, float_fmt='%6.7f')

"""
import sys
from struct import pack, unpack
from math import ceil
from collections import defaultdict
from typing import Tuple, Union

import numpy as np
from cpylog import get_logger2

from pyNastran.utils import is_binary_file, _filename, b


def read_cart3d(cart3d_filename, log=None, debug=False, result_names=None):
    """loads a Cart3D file"""
    model = Cart3D(log=log, debug=debug)
    model.read_cart3d(cart3d_filename, result_names)
    return model


def comp2tri(in_filenames, out_filename,
             is_binary=False, float_fmt='%6.7f',
             log=None, debug=False):
    """
    Combines multiple Cart3d files (binary or ascii) into a single file.

    Parameters
    ----------
    in_filenames : List[str]
        list of filenames
    out_filename : str
        output filename
    is_binary : bool; default=False
        is the output file binary
    float_fmt : str; default='%6.7f'
        the format string to use for ascii writing

    Notes
    -----
    assumes loads is None

    """
    points = []
    elements = []
    regions = []

    npoints = 0
    nregions = 0
    model = Cart3D(log=log, debug=debug)
    for infilename in in_filenames:
        model.read_cart3d(infilename)
        npointsi = model.nodes.shape[0]
        nregionsi = len(np.unique(model.regions))

        points.append(model.nodes)
        elements.append(model.elements + npoints)
        regions.append(model.regions + nregions)
        npoints += npointsi
        nregions += nregionsi

    points = np.vstack(points)
    elements = np.vstack(elements)
    regions = np.vstack(regions)
    model.points = points
    model.elements = elements
    model.regions = regions

    if out_filename:
        model.write_cart3d(out_filename, is_binary=is_binary, float_fmt=float_fmt)
    return model

class Cart3dIO:
    """
    Cart3d IO class
    """
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)
        self._endian = b''
        self._encoding = 'latin1'
        self.n = 0
        self.infile = None
        self.infilename = None
        self.points = None
        self.elements = None
        self.regions = None
        self.loads = {}

    def _write_header(self, outfile, points, elements, is_loads, is_binary=False):
        """
        writes the cart3d header

        Without results
        ---------------
        npoints nelements

        With results
        ------------
        npoints nelements nresults

        """
        npoints = points.shape[0]
        nelements = elements.shape[0]

        if is_binary:
            if is_loads:
                fmt = self._endian + b('iiiii')
                msg = pack(fmt, 3*4, npoints, nelements, 6, 4)
            else:
                fmt = self._endian + b('iiii')
                msg = pack(fmt, 2*4, npoints, nelements, 4)

            int_fmt = None
        else:
            # this is ASCII data
            if is_loads:
                msg = b("%i %i 6\n" % (npoints, nelements))
            else:
                msg = b("%i %i\n" % (npoints, nelements))

            # take the max value, string it, and length it
            # so 123,456 is length 6
            int_fmt = b('%%%si' % len(str(nelements)))
        outfile.write(msg)
        return int_fmt

    def _write_points(self, outfile, points, is_binary, float_fmt='%6.6f'):
        """writes the points"""
        if is_binary:
            four = pack(self._endian + b('i'), 4)
            outfile.write(four)

            npoints = points.shape[0]
            fmt = self._endian + b('%if' % (npoints * 3))
            floats = pack(fmt, *np.ravel(points))

            outfile.write(floats)
            outfile.write(four)
        else:
            if isinstance(float_fmt, bytes):
                fmt_ascii = float_fmt
            else:
                fmt_ascii = float_fmt.encode('latin1')
            np.savetxt(outfile, points, fmt_ascii)

    def _write_elements(self, outfile, elements, is_binary, int_fmt='%6i'):
        """writes the triangles"""
        min_e = elements.min()
        assert min_e == 0, 'min(elements)=%s' % min_e
        if is_binary:
            fmt = self._endian + b('i')
            four = pack(fmt, 4)
            outfile.write(four)
            nelements = elements.shape[0]
            fmt = self._endian + b('%ii' % (nelements * 3))
            ints = pack(fmt, *np.ravel(elements+1))

            outfile.write(ints)
            outfile.write(four)
        else:
            if isinstance(int_fmt, bytes):
                fmt_ascii = int_fmt
            else:
                fmt_ascii = int_fmt.encode('latin1')
            np.savetxt(outfile, elements+1, fmt_ascii)

    def _write_regions(self, outfile, regions, is_binary):
        """writes the regions"""
        if is_binary:
            fmt = self._endian + b('i')
            four = pack(fmt, 4)
            outfile.write(four)

            nregions = len(regions)
            fmt = self._endian + b('%ii' % nregions)
            ints = pack(fmt, *regions)
            outfile.write(ints)

            outfile.write(four)
        else:
            fmt = b'%i'
            np.savetxt(outfile, regions, fmt)

    def _write_loads(self, outfile, loads, is_binary, float_fmt='%6.6f'):
        """writes the *.triq loads"""
        if is_binary:
            raise NotImplementedError('is_binary=%s' % is_binary)
        else:
            Cp = loads['Cp']
            rho = loads['rho']
            rhoU = loads['rhoU']
            rhoV = loads['rhoV']
            rhoW = loads['rhoW']
            E = loads['E']
            npoints = self.points.shape[0]
            assert len(Cp) == npoints, 'len(Cp)=%s npoints=%s' % (len(Cp), npoints)
            #nrows = len(Cp)
            fmt = '%s\n%s %s %s %s %s\n' % (float_fmt, float_fmt, float_fmt,
                                            float_fmt, float_fmt, float_fmt)
            for (cpi, rhoi, rhou, rhov, rhoe, e) in zip(Cp, rho, rhoU, rhoV, rhoW, E):
                outfile.write(fmt % (cpi, rhoi, rhou, rhov, rhoe, e))

    def _read_header_ascii(self):
        """
        Reads the header::

          npoints nelements          # geometry
          npoints nelements nresults # results

        """
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
        """get the number of results"""
        if isinstance(self.loads, dict):
            return len(self.loads)
        return 0

    @property
    def nnodes(self):
        """alternate way to access number of points"""
        return self.npoints

    @property
    def npoints(self):
        """get the number of points"""
        return self.points.shape[0]

    @property
    def nodes(self):
        """alternate way to access the points"""
        return self.points

    @nodes.setter
    def nodes(self, points):
        """alternate way to access the points"""
        self.points = points

    @property
    def nelements(self):
        """get the number of elements"""
        return self.elements.shape[0]

    def _read_points_ascii(self, npoints):
        """
        A point is defined by x,y,z and the ID is the location in points.

        """
        p = 0
        data = []
        assert npoints > 0, 'npoints=%s' % npoints
        points = np.zeros((npoints, 3), dtype='float32')
        while p < npoints:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                x = data.pop(0)
                y = data.pop(0)
                z = data.pop(0)
                points[p] = [x, y, z]
                p += 1
        return points

    def _read_elements_ascii(self, nelements):
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.

        """
        assert nelements > 0, 'npoints=%s nelements=%s' % (self.npoints, nelements)
        elements = np.zeros((nelements, 3), dtype='int32')

        ieid = 0
        data = []
        while ieid < nelements:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                n1 = int(data.pop(0))
                n2 = int(data.pop(0))
                n3 = int(data.pop(0))
                elements[ieid] = [n1, n2, n3]
                ieid += 1

        nid_min = elements.min()
        if nid_min != 1:
            nid_max = elements.max()
            nnodes = self.nodes.shape[0]
            if nid_max == nnodes:
                msg = (
                    'Possible Cart3d error due to unused nodes\n'
                    'min(nids)=%s; expected 1; nid_max=%s nnodes=%s' % (
                        nid_min, nid_max, nnodes))
                self.log.warning(msg)
            else:
                msg = 'elements:\n%s\nmin(nids)=%s; expected 1; nid_max=%s nnodes=%s' % (
                    elements, nid_min, nid_max, nnodes, )
                raise RuntimeError(msg)
            #assert elements.min() == 1, elements.min()
        return elements - 1

    def _read_regions_ascii(self, nelements):
        """reads the region section"""
        regions = np.zeros(nelements, dtype='int32')
        iregion = 0
        data = []
        while iregion < nelements:
            data = self.infile.readline().strip().split()
            ndata = len(data)
            regions[iregion : iregion + ndata] = data
            iregion += ndata
        return regions

    def _read_header_binary(self):
        """
        Reads the header::

          npoints nelements          # geometry
          npoints nelements nresults # results

        """
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
            self._rewind()
            self.show(100)
            raise RuntimeError('unknown endian')

        self.n += 4
        data = self.infile.read(size)
        self.n += size

        so4 = size // 4  # size over 4
        if so4 == 3:
            (npoints, nelements, nresults) = unpack(self._endian + b('iii'), data)
            self.log.info("npoints=%s nelements=%s nresults=%s" % (npoints, nelements, nresults))
        elif so4 == 2:
            (npoints, nelements) = unpack(self._endian + b('ii'), data)
            nresults = 0
            self.log.info("npoints=%s nelements=%s" % (npoints, nelements))
        else:
            self._rewind()
            self.show(100)
            raise RuntimeError('in the wrong spot...endian...size/4=%s' % so4)
        self.infile.read(8)  # end of first block, start of second block
        return (npoints, nelements, nresults)

    def _read_points_binary(self, npoints):
        """reads the xyz points"""
        size = npoints * 12  # 12=3*4 all the points
        data = self.infile.read(size)

        dtype = np.dtype(self._endian + b('f4'))
        points = np.frombuffer(data, dtype=dtype).reshape((npoints, 3)).copy()

        self.infile.read(8)  # end of second block, start of third block
        return points

    def _read_elements_binary(self, nelements):
        """reads the triangles"""
        size = nelements * 12  # 12=3*4 all the elements
        data = self.infile.read(size)

        dtype = np.dtype(self._endian + b('i4'))
        elements = np.frombuffer(data, dtype=dtype).reshape((nelements, 3)).copy()

        self.infile.read(8)  # end of third (element) block, start of regions (fourth) block
        assert elements.min() == 1, elements.min()
        return elements - 1

    def _read_regions_binary(self, nelements):
        """reads the regions"""
        size = nelements * 4  # 12=3*4 all the elements
        data = self.infile.read(size)

        regions = np.zeros(nelements, dtype='int32')
        dtype = self._endian + b'i'
        regions = np.frombuffer(data, dtype=dtype).copy()

        self.infile.read(4)  # end of regions (fourth) block
        return regions

    def _read_results_binary(self, i, infile, result_names=None):
        """binary results are not supported"""
        pass

    def _rewind(self):  # pragma: no cover
        """go back to the beginning of the file"""
        self.n = 0
        self.infile.seek(self.n)

    def show(self, n, types='ifs', endian=None):  # pragma: no cover
        assert self.n == self.infile.tell(), 'n=%s tell=%s' % (self.n, self.infile.tell())
        #nints = n // 4
        data = self.infile.read(4 * n)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        self.infile.seek(self.n)
        return strings, ints, floats

    def show_data(self, data, types='ifs', endian=None):  # pragma: no cover
        return self._write_data(sys.stdout, data, types=types, endian=endian)

    def _write_data(self, outfile, data, types='ifs', endian=None):  # pragma: no cover
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

    def show_ndata(self, n, types='ifs'):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types)

    def _write_ndata(self, outfile, n, types='ifs'):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.

        """
        nold = self.n
        data = self.infile.read(n)
        self.n = nold
        self.infile.seek(self.n)
        return self._write_data(outfile, data, types=types)


class Cart3D(Cart3dIO):
    """Cart3d interface class"""
    model_type = 'cart3d'
    is_structured = False
    is_outward_normals = True

    def __init__(self, log=None, debug=False):
        Cart3dIO.__init__(self, log=log, debug=debug)
        self.loads = {}
        self.points = None
        self.elements = None

    def flip_model(self):
        """flip the model about the y-axis"""
        self.points[:, 1] *= -1.
        self.elements = np.hstack([
            self.elements[:, 0:1],
            self.elements[:, 2:3],
            self.elements[:, 1:2],
        ])
        print(self.elements.shape)

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
        #self.log.info('---starting make_mirror_model---')
        #assert tol >= 0, 'tol=%r' % tol #  prevents hacks to the axis

        #nnodes = nodes.shape[0]
        #assert nnodes > 0, 'nnodes=%s' % nnodes

        #nelements = elements.shape[0]
        #assert nelements > 0, 'nelements=%s' % nelements

        #ax, ax0 = self._get_ax(axis)
        #if ax in [0, 1, 2]:  # positive x, y, z values; mirror to -side
            #iy0 = np.where(nodes[:, ax] > tol)[0]
            #ax2 = ax
        #elif ax in [3, 4, 5]:  # negative x, y, z values; mirror to +side
            #iy0 = np.where(nodes[:, ax-3] < -tol)[0]
            #ax2 = ax - 3 # we create ax2 in order to generalize the data updating
        #else:
            #raise NotImplementedError(axis)

        ## the nodes to be duplicated are the nodes that aren't below the tolerance
        #nodes_upper = nodes[iy0]
        #nodes_upper[:, ax2] *= -1.0  # flip the nodes about the axis

        #nodes2 = np.vstack([nodes, nodes_upper])
        #nnodes2 = nodes2.shape[0]
        #assert nnodes2 > nnodes, 'nnodes2=%s nnodes=%s' % (nnodes2, nnodes)

        #nnodes_upper = nodes_upper.shape[0]
        #elements_upper = elements.copy()
        #nelements = elements.shape[0]

        ## remap the mirrored nodes with the new node ids
        #for eid in range(nelements):
            #element = elements_upper[eid, :]
            #for i, eidi in enumerate(element):
                #if eidi in iy0:
                    #elements_upper[eid][i] = nnodes_upper + eidi

                ## we need to reverse the element in order to get
                ## the proper normal vector
                #elements_upper[eid] = elements_upper[eid, ::-1]

        #elements2 = np.vstack([elements, elements_upper])
        #nelements2 = elements2.shape[0]
        #assert nelements2 > nelements, 'nelements2=%s nelements=%s' % (nelements2, nelements)

        #nregions = len(np.unique(regions))
        #regions_upper = regions.copy() + nregions
        #regions2 = np.hstack([regions, regions_upper])

        #loads2 = {}
        #for key, data in loads.items():

            ## flip the sign on the flipping axis terms
            #if((key in ['U', 'rhoU'] and ax2 == 0) or
               #(key in ['V', 'rhoV'] and ax2 == 1) or
               #(key in ['W', 'rhoW'] and ax2 == 2)):
                #data_upper = -data[iy0]
            #else:
                #data_upper = data[iy0]
            #loads2[key] = np.hstack([data, data_upper])

        #self.log.info('---finished make_mirror_model---')
        #return (nodes2, elements2, regions2, loads2)

    def _get_ax(self, axis: Union[str, int]) -> Tuple[int, int]:
        """helper method to convert an axis_string into an integer"""
        if isinstance(axis, str):
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
        else:  # pragma: no cover
            raise NotImplementedError('axis=%r' % axis)
        self.log.debug("axis=%r ax=%s" % (axis, ax))

        # shift ax to the actual column index in the nodes array
        ax0 = ax if ax in [0, 1, 2] else ax - 3
        return ax, ax0

    def make_half_model(self, axis='y', remap_nodes=True):
        """
        Makes a half model from a full model

        Notes
        -----
        Cp is really loads['Cp'] and was meant for loads analysis only

        """
        nodes = self.nodes
        elements = self.elements
        regions = self.regions
        loads = self.loads
        if loads is None:
            loads = {}

        nnodes = nodes.shape[0]
        assert nnodes > 0, 'nnodes=%s'  % nnodes

        nelements = elements.shape[0]
        assert nelements > 0, 'nelements=%s'  % nelements

        self.log.info('---starting make_half_model---')
        ax, ax0 = self._get_ax(axis)

        self.log.debug('find elements to remove')
        #print(f'axis={axis} ax={ax}')

        ynode = nodes[:, ax0]
        y1 = ynode[elements[:, 0]]
        y2 = ynode[elements[:, 1]]
        y3 = ynode[elements[:, 2]]
        if ax in [0, 1, 2]:  # keep values > 0
            ielements_save = np.where((y1 >= 0.0) & (y2 >= 0.0) & (y3 >= 0.0))[0]
        elif ax in [3, 4, 5]:  # keep values < 0
            ielements_save = np.where((y1 <= 0.0) & (y2 <= 0.0) & (y3 <= 0.0))[0]
        else:
            raise NotImplementedError('axis=%r ax=%s' % (axis, ax))

        elements2 = elements[ielements_save, :]
        regions2 = regions[ielements_save]

        nelements2 = elements2.shape[0]
        assert 0 < nelements2 < nelements, 'nelements=%s nelements2=%s'  % (nelements, nelements2)

        # location of nodes to keep
        inodes_save = np.unique(elements2.ravel())

        #------------------------------------------
        #checks

        #inodes_map = np.arange(len(inodes_save))
        is_nodes = 0 < len(inodes_save) < nnodes
        if not is_nodes:
            msg = 'There are no nodes in the model; len(inodes_save)=%s nnodes=%s'  % (
                len(inodes_save), nnodes)
            raise RuntimeError(msg)

        nodes2 = nodes[inodes_save, :]
        nnodes2 = nodes2.shape[0]
        assert 0 < nnodes2 < nnodes, 'no nodes were removed; nnodes=%s nnodes2=%s'  % (nnodes, nnodes2)
        #------------------------------------------
        # renumbers mesh

        # node id renumbering
        emin0 = elements2.min()
        emax0 = elements2.max()

        # only renumber if we need to
        if emin0 > 0 or emax0 != nelements2:
            # we're making a 0-based node map of old -> new node_id
            nodes_map = {}
            for i, nid in enumerate(inodes_save):
                nodes_map[nid] = i

            # update the node ids
            for ielement in range(nelements2):
                element2 = elements2[ielement, :]
                elements2[ielement, :] = [nodes_map[nid] for nid in element2]
        assert len(elements2) > 0, elements2.shape

        min_e = elements2.min()
        assert min_e == 0, 'min(elements)=%s' % min_e
        #------------------------------------------

        loads2 = {} # 'Cp', 'Mach', 'U', etc.
        for key, load in loads.items():
            loads2[key] = load[inodes_save]

        self.log.info('---finished make_half_model---')
        self.nodes = nodes2
        self.elements = elements2
        self.regions = regions2
        self.loads = loads2
        return (nodes2, elements2, regions2, loads2)

    def get_free_edges(self, elements):
        """
        Cart3d must be a closed model with each edge shared by 2 elements
        The free edges indicate the problematic areas.

        Returns
        -------
        free edges : (nedges, 2) int ndarray
            the free edge node ids

        """
        edge_to_eid_map = defaultdict(list)
        for i, element in enumerate(elements):
            edge1 = tuple(sorted([element[0], element[1]]))
            edge2 = tuple(sorted([element[1], element[2]]))
            edge3 = tuple(sorted([element[2], element[0]]))
            edge_to_eid_map[edge1].append(i)
            edge_to_eid_map[edge2].append(i)
            edge_to_eid_map[edge3].append(i)

        free_edges = []
        for edge, eids in sorted(edge_to_eid_map.items()):
            if len(eids) != 2:
                free_edges.append(edge)
        return np.array(free_edges, dtype='int32')

    def read_cart3d(self, infilename, result_names=None):
        """extracts the points, elements, and Cp"""
        self.infilename = infilename
        self.log.info("---reading cart3d...%r---" % self.infilename)

        self.infilename = infilename
        if is_binary_file(infilename):
            with open(infilename, 'rb') as self.infile:
                try:
                    npoints, nelements, nresults = self._read_header_binary()
                    self.points = self._read_points_binary(npoints)
                    self.elements = self._read_elements_binary(nelements)
                    self.regions = self._read_regions_binary(nelements)
                    # TODO: loads
                except:
                    msg = 'failed reading %r' % infilename
                    self.log.error(msg)
                    raise

        else:
            with open(_filename(infilename), 'r', encoding=self._encoding) as self.infile:
                try:
                    npoints, nelements, nresults = self._read_header_ascii()
                    self.points = self._read_points_ascii(npoints)
                    self.elements = self._read_elements_ascii(nelements)
                    self.regions = self._read_regions_ascii(nelements)
                    self._read_results_ascii(0, self.infile, nresults, result_names=result_names)
                except:
                    msg = 'failed reading %r' % infilename
                    self.log.error(msg)
                    raise

        self.log.debug("npoints=%s nelements=%s" % (self.npoints, self.nelements))
        assert self.npoints > 0, 'npoints=%s' % self.npoints
        assert self.nelements > 0, 'nelements=%s' % self.nelements

    def write_cart3d(self, outfilename, is_binary=False, float_fmt='%6.7f'):
        """
        writes a cart3d file

        """
        assert len(self.points) > 0, 'len(self.points)=%s' % len(self.points)

        if self.loads is None or self.loads == {}:
            #loads = {}
            is_loads = False
        else:
            is_loads = True

        self.log.info("---writing cart3d...%r---" % outfilename)
        with open(outfilename, 'wb') as outfile:
            int_fmt = self._write_header(outfile, self.points, self.elements, is_loads, is_binary)
            self._write_points(outfile, self.points, is_binary, float_fmt)
            self._write_elements(outfile, self.elements, is_binary, int_fmt)
            self._write_regions(outfile, self.regions, is_binary)

            if is_loads:
                assert is_binary is False, 'is_binary=%r is not supported for loads' % is_binary
                self._write_loads(outfile, self.loads, is_binary, float_fmt)

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
        #loads = {}
        if result_names is None:
            result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                            'Mach', 'U', 'V', 'W', 'E', 'a', 'T', 'Pressure', 'q']
        self.log.debug('---starting read_results---')

        results = np.zeros((self.npoints, 6), dtype='float32')

        nresult_lines = int(ceil(nresults / 5.)) - 1
        for ipoint in range(self.npoints):
            # rho rhoU,rhoV,rhoW,pressure/rhoE/E
            sline = infile.readline().strip().split()
            i += 1
            for unused_n in range(nresult_lines):
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
                    #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (
                        #pointNum, cp, mach, rho, rhoU, rhoV, rhoW))
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
        loads : dict; default=None -> {}
            key : ???
            value : ???

        """
        if loads is None:
            loads = {}
        Cp = results[:, 0]
        rho = results[:, 1]
        rho_u = results[:, 2]
        rho_v = results[:, 3]
        rho_w = results[:, 4]
        E = results[:, 5]

        ibad = np.where(rho <= 0.000001)[0]
        if len(ibad) > 0:

            if 'Mach' in result_names:
                Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2)# / rho
                Mach[ibad] = 0.0
            if 'U' in result_names:
                U = rho_u / rho
                U[ibad] = 0.0
            if 'U' in result_names:
                V = rho_v / rho
                V[ibad] = 0.0
            if 'W' in result_names:
                W = rho_w / rho
                W[ibad] = 0.0
            #if 'rhoE' in result_names:
                #rho_e = rhoE / rho
                #e[ibad] = 0.0

            is_bad = True
            #n = 0
            #for i in ibad:
                #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (
                    #i, Cp[i], Mach[i], rho[i], rho_u[i], rho_v[i], rho_w[i]))
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
            loads['rhoU'] = rho_u
        if 'rhoV' in result_names:
            loads['rhoV'] = rho_v
        if 'rhoW' in result_names:
            loads['rhoW'] = rho_w
        #if 'rhoE' in result_names:
            #loads['rhoE'] = rho_e

        if 'rho' in result_names:
            loads['rho'] = rho

        if 'Mach' in result_names:
            if not is_bad:
                #Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2) / rho
                Mach = np.sqrt(rho_u**2 + rho_v**2 + rho_w**2)
            loads['Mach'] = Mach

        if 'U' in result_names:
            if not is_bad:
                U = rho_u / rho
            loads['U'] = U
        if 'V' in result_names:
            if not is_bad:
                V = rho_v / rho
            loads['V'] = V
        if 'W' in result_names:
            if not is_bad:
                W = rho_w / rho
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
            #print('T: min=%s max=%s' % (T.min(), T.max()))
            loads['a'] = np.sqrt(T)
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
        #self.log.info("i=%s Cp=%s rho=%s rho_u=%s rho_v=%s rho_w=%s Mach=%s" % (
            #i, Cp[i], rho[i], rho_u[i], rho_v[i], rho_w[i], Mach[i]))
        self.log.debug('---finished read_results---')
        return loads

    def _get_area_vector(self):
        """
        Gets the area vector (unnormalized normal vector)
        Returns
        -------
        normals : (n, 3) ndarray
            unnormalized centroidal normal vectors

        """
        elements = self.elements
        nodes = self.nodes
        p1 = nodes[elements[:, 0], :]
        p2 = nodes[elements[:, 1], :]
        p3 = nodes[elements[:, 2], :]

        ne = elements.shape[0]
        avec = p2 - p1
        bvec = p3 - p1
        n = np.cross(avec, bvec)
        assert len(n) == ne, 'len(n)=%s ne=%s' % (len(n), ne)

        return n

    def get_area(self):
        """
        Gets the element area

        Returns
        -------
        area : (n, 3) ndarray
            the element areas
        """

        ne = self.elements.shape[0]
        n = self._get_area_vector()
        ni = np.linalg.norm(n, axis=1)
        assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)
        return 0.5 * ni

    def get_normals(self):
        """
        Gets the centroidal normals

        Returns
        -------
        cnormals : (n, 3) ndarray
            normalized centroidal normal vectors

        """
        ne = self.elements.shape[0]
        n = self._get_area_vector()
        ni = np.linalg.norm(n, axis=1)
        assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)

        assert ni.min() > 0.0, ni[np.where(ni <= 0.0)[0]]
        n /= ni[:, None]  # normal vector
        return n

    def get_normals_at_nodes(self, cnormals):
        """
        Gets the nodal normals

        Parameters
        ----------
        cnormals : (n, 3) ndarray
            normalized centroidal normal vectors

        Returns
        -------
        nnormals : (n, 3) ndarray
            normalized nodal normal vectors

        """
        elements = self.elements
        #nodes = self.nodes
        nnodes = self.nnodes
        nid_to_eids = defaultdict(list)

        # find the elements to consider for each node
        for eid, element in enumerate(elements):
            n1, n2, n3 = element
            nid_to_eids[n1].append(eid)
            nid_to_eids[n2].append(eid)
            nid_to_eids[n3].append(eid)

        nnormals = np.zeros((nnodes, 3), dtype='float64')
        for nid in range(nnodes):
            eids = nid_to_eids[nid]
            if len(eids) == 0:
                raise RuntimeError('nid=%s is not used' % nid)
            #ni_avg = cnormals[eids, :]
            nnormals[nid] = cnormals[eids, :].sum(axis=0)
        ni = np.linalg.norm(nnormals, axis=1)
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
