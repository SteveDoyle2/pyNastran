"""
models from:
    http://people.sc.fsu.edu/~jburkardt/data/tec/tec.html
"""
import sys
import os
from struct import unpack
import itertools
from typing import List

import numpy as np
from cpylog import get_logger2

from pyNastran.utils import is_binary_file
from pyNastran.converters.tecplot.zone import Zone, CaseInsensitiveDict, is_3d

def read_tecplot(tecplot_filename: str, use_cols=None, dtype=None, log=None, debug=False):
    """loads a tecplot file"""
    tecplot = Tecplot(log=log, debug=debug)
    if use_cols:
        tecplot.use_cols = use_cols
        tecplot.dtype = dtype
    tecplot.read_tecplot(tecplot_filename)
    return tecplot


class Tecplot:
    """
    Parses a hexa binary/ASCII Tecplot 360 file.
    Writes an ASCII Tecplot 10 file.

    Supports:
     - title
     - single zone only?
     - unstructured:
       - nodal results
       - single element type; ZONETYPE = [FETRIANGLE, FEQUADRILATERAL, FETETRAHEDRON, FEBRICK]
       - DATAPACKING = [POINT, ELEMENT] writing
       - 2D/3D
       - full support for writing
     - structured:
       - nodal results
       - F=POINT
       - 2d/3d (POINT, I/J/K)
       - full support for writing
       - no reshaping of xyz to make slicing easier!

    Doesn't support:
     - text
     - geometry
     - transient writing
     - centroidal results
     - non-sequential node ids
     - data lists (100*0.0)
    """
    def __repr__(self):
        msg = ''
        for zone in self.zones:
            msg += str(zone)
        return msg

    def __init__(self, log=None, debug=False):
        # defines binary file specific features
        self._endian = b'<'
        self._n = 0

        self.tecplot_filename = ''
        self.log = get_logger2(log, debug=debug)
        self.debug = debug

        # mesh = None : model hasn't been read
        self.is_mesh = None

        self.title = 'tecplot geometry and solution file'
        self.variables = None

        self.zones = []
        # mesh = True : this is a structured/unstructured grid

        # mesh = False : this is a plot file
        self.use_cols = None
        self.dtype = None

        self._uendian = ''
        self.n = 0

    @property
    def nzones(self):
        """gets the number of zones"""
        return len(self.zones)

    @property
    def headers_dict(self):
        raise RuntimeError('this data member has been removed')
    @property
    def xy(self):
        raise RuntimeError('this data member has been removed')
    @property
    def xyz(self):
        raise RuntimeError('this data member has been removed')
    @property
    def tri_elements(self):
        raise RuntimeError('this data member has been removed')
    @property
    def quad_elements(self):
        raise RuntimeError('this data member has been removed')
    @property
    def tet_elements(self):
        raise RuntimeError('this data member has been removed')
    @property
    def hexa_elements(self):
        raise RuntimeError('this data member has been removed')

    @headers_dict.setter
    def headers_dict(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @xy.setter
    def xy(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @xyz.setter
    def xyz(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @tri_elements.setter
    def tri_elements(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @quad_elements.setter
    def quad_elements(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @tet_elements.setter
    def tet_elements(self, unused_x):
        raise RuntimeError('this data member has been removed')
    @hexa_elements.setter
    def hexa_elements(self, unused_x):
        raise RuntimeError('this data member has been removed')

    def read_tecplot(self, tecplot_filename):
        """
        Reads an ASCII/binary Tecplot file.

        The binary file reader must have ONLY CHEXAs and be Tecplot 360 with
        `rho`, `u`, `v`, `w`, and `p`.
        The ASCII file reader has only been tested with Tecplot 10, but will
        probably work on Tecplot360.  It **should** work with any set of
        variables.
        """
        if is_binary_file(tecplot_filename):
            return self.read_tecplot_binary(tecplot_filename)
        return self.read_tecplot_ascii(tecplot_filename)

    def read_tecplot_ascii(self, tecplot_filename, nnodes=None, nelements=None):
        """
        Reads a Tecplot ASCII file.

        Supports:
         - CTRIA3
         - CQUAD4
         - CTETRA
         - CHEXA

        .. note :: assumes single typed results
        .. warning:: BLOCK option doesn't work if line length isn't the same...
        """
        self.tecplot_filename = tecplot_filename
        assert os.path.exists(tecplot_filename), tecplot_filename

        iline = 0
        nnodes = -1
        nelements = -1
        with open(tecplot_filename, 'r') as tecplot_file:
            lines = tecplot_file.readlines()
        del tecplot_file

        quads_list = []
        hexas_list = []
        tris_list = []
        tets_list = []
        xyz_list = []
        results_list = []

        line = lines[iline].strip()
        iline += 1
        iblock = 0
        while 1:
            #print('start...')
            iline, title_line, header_lines, line = _read_header_lines(
                lines, iline, line, self.log)
            #print('header_lines', header_lines)
            headers_dict = _header_lines_to_header_dict(title_line, header_lines, self.variables)
            if headers_dict is None:
                break
            zone = Zone(self.log)
            zone.headers_dict = headers_dict
            self.variables = headers_dict['VARIABLES']
            #print('self.variables', self.variables)

            #print(headers_dict.keys())
            if 'ZONETYPE' in headers_dict:
                zone_type = headers_dict['ZONETYPE'].upper() # FEBrick
                data_packing = headers_dict['DATAPACKING'].upper() # block
                iline = self._read_zonetype(
                    zone, zone_type, lines, iline, iblock, headers_dict, line,
                    nnodes, nelements,
                    xyz_list, hexas_list, tets_list, quads_list, tris_list,
                    results_list,
                    data_packing=data_packing)
            elif 'F' in headers_dict:
                fe = headers_dict['F'] # FEPoint
                assert isinstance(fe, str), headers_dict
                zone_type = fe.upper() # FEPoint
                self.log.debug('zone_type = %r' % zone_type[0])
                iline = self._read_zonetype(
                    zone, zone_type, lines, iline, iblock, headers_dict, line,
                    nnodes, nelements,
                    xyz_list, hexas_list, tets_list, quads_list, tris_list,
                    results_list,
                    fe=fe)
                iline -= 1
            elif (('ZONE' in headers_dict) and
                  (headers_dict['ZONE'] is None) and
                  ('T' in headers_dict)):
                A, line = self.read_table(lines, iline, iblock, headers_dict, line)
                self.A = A
                #print('read_table...')
                return
            else:
                msg = 'headers=%s\n' % str(headers_dict)
                msg += 'line = %r' % line.strip()
                raise NotImplementedError(msg)

            self.zones.append(zone)

            #sline = line.split()
            #print('stack...')
            _stack(zone, xyz_list, quads_list, tris_list, tets_list, hexas_list, results_list, self.log)
            #print(zone)
            if line is None:
                return
            try:
                line = lines[iline].strip()
            except IndexError:
                break
            quads_list = []
            hexas_list = []
            tris_list = []
            tets_list = []
            xyz_list = []
            results_list = []

    def read_table(self, tecplot_file, unused_iblock, headers_dict, line):
        """
        reads a space-separated tabular data block
        """
        variables = [var.strip('" ') for var in headers_dict['VARIABLES']]
        #print('variables = %s' % variables)
        #self.dtype[]
        use_cols = [variables.index(var) for var in self.use_cols]

        # add on the preceding line to the line "list"
        # that's not a hack at all...
        lines = itertools.chain((line, ), iter(tecplot_file))
        A = np.loadtxt(lines, dtype=self.dtype, comments='#', delimiter=None,
                       converters=None, skiprows=0,
                       usecols=use_cols, unpack=False, ndmin=0)
        return A, None

    def _read_zonetype(self, zone, zone_type, lines, iline, iblock, headers_dict, line,
                       nnodes, nelements,
                       xyz_list, hexas_list, tets_list, quads_list, tris_list,
                       results_list,
                       data_packing=None, fe=None):
        """
        Parameters
        ----------
        zone_type : str

        fe : str
          - a zone_type.upper() string???
          - FEPOINT

        reads:
          - ZONE E
          - ZONE T

        ZONE is a flag, T is title, E is number of elements

        -------------  ---------  ----------  ----------------------------------------------
        Parameter      Ordered    Finite      Description
                       Data       Element
        -------------  ---------  ----------  ----------------------------------------------
        T="title"      Yes        Yes         Zone title.
        I=imax         Yes        No          Number of points in 1st dimension.
        J=jmax         Yes        No          Number of points in 2nd dimension.
        K=kmax         Yes        No          Number of points in 3rd dimension.
        C=colour       Yes        Yes         Colour from WHITE, BLACK, RED, GREEN,
                                                          BLUE, CYAN, YELLOW, PURPLE,
                                                          CUST1, CUST2,....CUST8.
        F=format       Yes        Yes         POINT or BLOCK for ordered data.
                                              FEPOINT or FEBLOCK for finite element.
        D=(list)       Yes        Yes         A list of variable names to to include
                                              from the last zone.
        DT=(list)      Yes        Yes         A list of datatypes for each variable.
                                              SINGLE, DOUBLE, LONGINT, SHORTINT, BYTE, BIT.
        N=num          No         Yes         Number of nodes.
        E=num          No         Yes         Number of elements.
        ET=type        No         Yes         Element type from TRIANGLE, BRICK,
                                                                QUADRILATERAL, TETRAHEDRON.
        NV=variable    No         Yes         Variable for node value.
        -------------  ---------  ----------  ----------------------------------------------
        http://paulbourke.net/dataformats/tp/
        """
        #print('self.variables', self.variables)
        #ndim = zone.ndim
        #print('iblock =', iblock)
        if iblock == 0:
            variables = headers_dict['VARIABLES']
            zone.variables = [variable.strip(' \r\n\t"\'') for variable in variables]
            self.log.debug('zone.variables = %s' % zone.variables)
            nresults = len(variables) - 3 # x, y, z, rho, u, v, w, p
            self.log.debug('nresults = %s' % nresults)

        self.log.debug(str(headers_dict))
        is_unstructured = False
        is_structured = False
        if zone_type in ['FETRIANGLE', 'FEQUADRILATERAL', 'FETETRAHEDRON', 'FEBRICK']:
            nnodesi = headers_dict['N']
            nelementsi = headers_dict['E']
            is_unstructured = True
        elif zone_type in ['POINT', 'BLOCK']: #  structured
            ni = headers_dict['I']
            if 'J' in headers_dict:
                nj = headers_dict['J']
                if 'K' in headers_dict:
                    # 3d
                    nk = headers_dict['K']
                    nnodesi = ni * nj * nk
                    nelementsi = (ni - 1) * (nj - 1) * (nk - 1)
                else:
                    # 2d
                    nnodesi = ni * nj
                    nelementsi = (ni - 1) * (nj - 1)
            else:
                assert 'K' not in headers_dict, list(headers_dict.keys())
                nnodesi = ni
                nelementsi = (ni - 1)
            assert nelementsi >= 0, nelementsi
            #nelementsi = 0
            elements = None # np.zeros((nelementsi, 8), dtype='int32')
            is_structured = True
        else:
            raise NotImplementedError('zone_type = %r' % zone_type)
        self.log.info(f'zone_type={zone_type} data_packing={data_packing} '
                      f'nnodes={nnodesi} nelements={nelementsi}')

        assert nnodesi > 0, nnodesi
        assert nresults >= 0, 'nresults=%s' % nresults
        xyz = np.zeros((nnodesi, 3), dtype='float32')
        results = np.zeros((nnodesi, nresults), dtype='float32')
        if zone_type == 'FEBRICK':
            # hex
            elements = np.zeros((nelementsi, 8), dtype='int32')
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETETRAHEDRON'):
            # quads / tets
            elements = np.zeros((nelementsi, 4), dtype='int32')
        elif zone_type == 'FETRIANGLE':
            # tris
            elements = np.zeros((nelementsi, 3), dtype='int32')
        #elif zone_type == 'FEBLOCK':
            #pass
        elif  zone_type in ['POINT', 'BLOCK']:
            # already handled
            #print('data')
            pass
        else:
            #if isinstance(zone_type, list):
                #raise NotImplementedError(zone_type[0])
            raise NotImplementedError(zone_type)

        sline = split_line(line.strip())
        if zone_type in ('FEBRICK', 'FETETRAHEDRON'):
            if data_packing == 'POINT':
                for inode in range(nnodesi):
                    if inode == 0:
                        self.log.debug('zone_type=%s sline=%s' %(zone_type, sline))
                    if not len(sline[3:]) == len(results[inode, :]):
                        msg = 'sline[3:]=%s results[inode, :]=%s' % (sline[:3], results[inode, :])
                        raise RuntimeError(msg)

                    try:
                        xyz[inode, :] = sline[:3]
                        results[inode, :] = sline[3:]
                    except ValueError:
                        msg = 'i=%s line=%r\n' % (inode, line)
                        msg += 'sline = %s' % str(sline)
                        print(msg)
                        raise

                    iline, line, sline = get_next_sline(lines, iline)
            elif data_packing == 'BLOCK':
                iline, line, sline = read_zone_block(lines, iline, xyz, results, nresults, zone_type,
                                                     sline, nnodesi, self.log)
                #print('sline =', sline)
            else:
                raise NotImplementedError(data_packing)
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETRIANGLE'):
            sline = split_line(line.strip())
            for inode in range(nnodesi):
                #print(iline, inode, sline)
                xyz[inode, :] = sline[:3]
                #if abs(xyz[inode, 1]) <= 5.0:
                    #msg = 'inode=%s xyz=%s'  % (inode, xyz[inode, :])
                    #raise RuntimeError(msg)

                results[inode, :] = sline[3:]
                iline, line, sline = get_next_sline(lines, iline)

        elif zone_type == 'POINT':
            nvars = len(zone.variables)
            iline, line, sline = read_point(lines, iline, xyz, results, zone_type,
                                            line, sline, nnodesi, nvars, self.log)
        elif zone_type == 'BLOCK':
            nvars = len(zone.variables)
            iline = read_block(lines, iline, xyz, results, zone_type,
                               line, sline, nnodesi, nvars, self.log)
        else:  # pragma: no cover
            raise NotImplementedError(zone_type)

        #print(elements.shape)
        #print('xyz[0 , :]', xyz[0, :])
        #print('xyz[-1, :]', xyz[-1, :])
        #print(sline)
        if is_structured:
            pass
        elif is_unstructured:
            iline, line, sline = read_unstructured_elements(lines, iline, sline, elements, nelementsi)

            #print(f.readline())

            if zone_type == 'FEBRICK':
                hexas_list.append(elements + nnodes)
            elif zone_type == 'FETETRAHEDRON':
                tets_list.append(elements + nnodes)
            elif zone_type in ('FEPOINT', 'FEQUADRILATERAL'):
                # TODO: why are points stuck in the quads?
                quads_list.append(elements + nnodes)
            elif zone_type == 'FETRIANGLE':
                tris_list.append(elements + nnodes)
            else:
                raise NotImplementedError(zone_type)
        else:
            raise RuntimeError()
        xyz_list.append(xyz)
        results_list.append(results)
        nnodes += nnodesi
        nelements += nelementsi
        self.log.debug('nnodes=%s nelements=%s (0-based)' % (nnodes, nelements))
        del headers_dict
        iblock += 1
        if iblock == 10:
            return
        self.log.debug('final sline=%s' % sline)
        return iline

    def read_tecplot_binary(self, tecplot_filename, nnodes=None,
                            nelements=None):
        """
        The binary file reader must have ONLY CHEXAs and be Tecplot 360
        with:
        `rho`, `u`, `v`, `w`, and `p`.
        """
        self.tecplot_filename = tecplot_filename
        assert os.path.exists(tecplot_filename), tecplot_filename
        with open(tecplot_filename, 'rb') as tecplot_file:
            self.f = tecplot_file
            self._uendian = '<'
            self.n = 0
            self.variables = ['rho', 'u', 'v', 'w', 'p']

            data = tecplot_file.read(8)
            self.n += 8
            word, = unpack(b'8s', data)
            self.log.debug('word = %r' % word)

            #self.show(100, endian='<')

            # http://home.ustc.edu.cn/~cbq/360_data_format_guide.pdf
            # page 151
            if 1:
                values = []
                ii = 0
                for ii in range(100):
                    datai = tecplot_file.read(4)
                    vali, = unpack(b'i', datai)
                    valf, = unpack(b'f', datai)
                    self.n += 4
                    values.append((vali, valf))
                    if vali == 9999:
                        #print('breaking...')
                        break
                #for j, vals in enumerate(values):
                    #print('  ', j, vals)
                assert ii < 100, ii

                nbytes = 3 * 4
                data = tecplot_file.read(nbytes)
                self.n += nbytes
                self.show_data(data, types='if', endian='<')

            nbytes = 1 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes
            zone_type, = unpack(b'i', data)
            self.log.debug('zone_type = %s' % zone_type)
            self.show(100, types='if', endian='<')

            nbytes = 11 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes
            #self.show_data(data, types='if', endian='<') # 'if'?
            s = unpack('2f 9i', data)
            self.log.debug(str(s))
            #assert self.n == 360, self.n
            #print('----------')

            nbytes = 2 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes
            nnodes2, nelements2 = unpack('2i', data)
            assert nnodes2 > 0, nnodes2
            assert nelements2 > 0, nelements2
            #self.show_data(data, types='if', endian='<') # 'if'?
            if nnodes and nelements:
                self.log.debug('nnodes=%s nelements=%s' % (nnodes, nelements))
                self.log.debug('nnodes2=%s nelements2=%s' % (nnodes2, nelements2))
            else:
                nnodes = nnodes2
                nelements = nelements2

            self.log.info('nnodes=%s nelements=%s' % (nnodes, nelements))
            assert nnodes == nnodes2
            assert nelements == nelements2
            #assert nnodes2 < 10000, nnodes
            #assert nelements2 < 10000, nelements

            nbytes = 35 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes
            #self.show_data(data, types='if', endian='<')
            #print('----------')

            nbytes = 30 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes

            #print('----------------------')
            #self.show_data(data, types='if', endian='<')
            #print('----------------------')

            # 0 - ORDERED (meaning?)
            # 1 - FELINESEG (meaning?)
            # 2 - FETRIANGLE
            # 3 - FEQUADRILATERAL
            # 4 - FETETRAHEDRON
            # 5 - FEBRICK
            assert zone_type in [0, 1, 2, 3, 4, 5], zone_type

            # p.93
            # zone_title
            # zone_type
            #   0 = ORDERED
            #   1 = FELINESEG
            #   2 = FETRIANGLE
            #   3 = FEQUADRILATERAL
            #   4 = FETETRAHEDRON
            #   5 = FEBRICK
            # i_max_or_num_points
            # j_max_or_num_elements
            # k_max
            # i_cell_max
            # j_cell_max
            # k_cell_max
            # solution_time
            # strand_id
            # parent_zone
            # is_block (0=POINT, 1=BLOCK)
            # num_face_connections
            # face_neighbor_mode
            # passive_var_list
            # value_location (0=cell-centered; 1=node-centered)
            # share_var_from_zone
            # share_connectivity_from_zone

            # http://www.hgs.k12.va.us/tecplot/documentation/tp_data_format_guide.pdf

            # 0=POINT
            # 1=BLOCK
            is_block = False

            # value_location:
            #   0=cell-centered
            #   1=node-centered
            #value_location = None
            if is_block:
                raise NotImplementedError('is_block=%s' % is_block)
            else:
                # is_point
                #print('----------')
                # the variables: [x, y, z]
                nvars = 3
                #nnodes = 3807

                ni = nnodes * nvars
                nbytes = ni * 4
                #print('nbytes =', nbytes)
                data = tecplot_file.read(nbytes)
                self.n += nbytes
                xyzvals = unpack('%if' % ni, data)
                xyz = np.array(xyzvals, dtype='float32').reshape(3, nnodes).T

                # the variables: [rho, u, v, w, p]
                nvars = 5
                dunno = 0    # what's with this...
                ni = nnodes * nvars + dunno
                nbytes = ni * 4
                data = tecplot_file.read(nbytes)
                self.n += nbytes
                resvals = unpack('%if' % ni, data)
                nodal_results = np.array(resvals, dtype='float32').reshape(nvars, nnodes).T


                # 7443 elements
                if zone_type == 5:
                    # CHEXA
                    nnodes_per_element = 8 # 8 nodes/elements
                    nvals = nnodes_per_element * nelements
                #elif zone_type == 1:
                    #asdf
                elif zone_type == 0:
                    # CQUAD4
                    nnodes_per_element = 4
                    nvals = nnodes_per_element * nelements
                    self.log.debug('nvals = %s' % nvals)
                else:
                    raise NotImplementedError('zone_type=%s' % zone_type)

                nbytes = nvals * 4
                node_ids = unpack(b'%ii' % nvals, tecplot_file.read(nbytes))
                self.n += nbytes

                elements = np.array(node_ids).reshape(nelements, nnodes_per_element)
                #print(elements)

                #self.show_data(data, types='ifs', endian='<')
                #print(vals)

                #self.show(100, endian='<')
                zone = Zone(self.log)
                if zone_type == 5:
                    zone.hexa_elements = elements
                elif zone_type == 0:
                    zone.quad_elements = elements
                else:
                    raise NotImplementedError(zone_type)
            del self.f

        zone.xyz = xyz
        zone.nodal_results = nodal_results
        self.zones.append(zone)
        #self.log.debug('done...')

    def show(self, n, types='ifs', endian=None):  # pragma: no cover
        assert self.n == self.f.tell()
        nints = n // 4
        data = self.f.read(4 * nints)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        self.f.seek(self.n)
        return strings, ints, floats

    def show_data(self, data, types='ifs', endian=None):  # pragma: no cover
        """
        Shows a data block as various types

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        .. warning:: 's' is apparently not Python 3 friendly

        """
        return self._write_data(sys.stdout, data, types=types, endian=endian)

    def _write_data(self, f, data, types='ifs', endian=None):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._uendian
            assert endian is not None, endian

        f.write('\nndata = %s:\n' % n)
        for typei in types:
            assert typei in 'sifdq lIL', 'type=%r is invalid' % typei

        if 's' in types:
            strings = unpack('%s%is' % (endian, n), data)
            f.write("  strings = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack('%s%ii' % (endian, nints), data)
            f.write("  ints    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack('%s%if' % (endian, nints), data)
            f.write("  floats  = %s\n" % str(floats))
        if 'd' in types:
            doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
            f.write("  doubles (float64) = %s\n" % str(doubles))

        if 'l' in types:
            longs = unpack('%s%il' % (endian, nints), data)
            f.write("  long  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack('%s%iI' % (endian, nints), data)
            f.write("  unsigned int = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack('%s%iL' % (endian, nints), data)
            f.write("  unsigned long = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write("  long long (int64) = %s\n" % str(longs))
        f.write('\n')
        return strings, ints, floats

    def show_ndata(self, n, types='ifs'):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types)

    def _write_ndata(self, f, n, types='ifs'):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.
        """
        nold = self.n
        data = self.f.read(n)
        self.n = nold
        self.f.seek(self.n)
        return self._write_data(f, data, types=types)

    def slice_x(self, xslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        zone = self.zones[0]
        x = zone.xyz[:, 0]
        self._slice_plane(zone, x, xslice)

    def slice_y(self, yslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        zone = self.zones[0]
        y = zone.xyz[:, 1]
        self._slice_plane(zone, y, yslice)

    def slice_z(self, zslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        zone = self.zones[0]
        z = zone.xyz[:, 2]
        self._slice_plane(zone, z, zslice)

    def slice_xyz(self, xslice, yslice, zslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        zone = self.zones[0]
        x = zone.xyz[:, 0]
        y = zone.xyz[:, 1]
        z = zone.xyz[:, 2]

        inodes = []
        if xslice is not None:
            xslice = float(xslice)
            inodes.append(np.where(x < xslice)[0])
        if yslice is not None:
            yslice = float(yslice)
            inodes.append(np.where(y < yslice)[0])
        if zslice is not None:
            zslice = float(zslice)
            inodes.append(np.where(z < zslice)[0])

        nodes = None
        if len(inodes) == 1:
            nodes = inodes[0]
        elif len(inodes) == 2:
            nodes = np.intersect1d(inodes[0], inodes[1], assume_unique=True)
        elif len(inodes) == 3:
            nodes = np.intersect1d(
                np.intersect1d(inodes[0], inodes[1], assume_unique=True),
                inodes[2], assume_unique=True)
            #inodes = arange(self.nodes.shape[0])
            # nodes = unique(hstack(inodes))
        if nodes is not None:
            zone._slice_plane_inodes(nodes)

    def _slice_plane(self, zone, y, slice_value):
        """
        - Only works for CHEXA
        - Doesn't remove unused nodes/renumber elements
        """
        slice_value = float(slice_value)
        inodes = np.where(y < slice_value)[0]
        zone._slice_plane_inodes(inodes)

    def _get_write_header(self, res_types):
        """gets the tecplot header"""
        is_y = True
        is_z = True
        #is_results = False
        assert self.nzones >= 1, self.nzones
        for zone in self.zones:
            variables = zone.variables
            is_y = 'Y' in zone.headers_dict['VARIABLES']
            is_z = 'Z' in zone.headers_dict['VARIABLES']
            is_results = bool(len(zone.nodal_results))
            if res_types is None:
                res_types = zone.variables
            elif isinstance(res_types, str):
                res_types = [res_types]
            break

        #"tecplot geometry and solution file"
        title = self.title
        if '"' in title or "'" in title:
            msg = 'TITLE = %s\n' % self.title
        else:
            msg = 'TITLE = "%s"\n' % self.title
        msg += 'VARIABLES = "X"\n'
        if is_y:
            msg += '"Y"\n'
        if is_z:
            msg += '"Z"\n'

        result_indices_to_write = []
        if is_results:
            #msg += '"rho"\n'
            #msg += '"u"\n'
            #msg += '"v"\n'
            #msg += '"w"\n'
            #msg += '"p"\n'
            #msg += 'ZONE T="%s"\n' % r'\"processor 1\"'
            #print('res_types =', res_types)
            #print('vars =', variables)
            for ivar, var in enumerate(res_types):
                if var not in variables:
                    raise RuntimeError('var=%r not in variables=%s' % (var, variables))
                #print('adding %s' % var)
                result_indices_to_write.append(variables.index(var))
            ivars = np.unique(result_indices_to_write)
            ivars.sort()
            for ivar in ivars:
                var = variables[ivar]
                msg += '"%s"\n' % var
            #print('ivars =', ivars)
        else:
            #if res_types is None:
            assert len(res_types) == 0, len(res_types)
            ivars = []
        return msg, ivars

    def write_tecplot(self, tecplot_filename, res_types=None, adjust_nids=True):
        """
        Only handles single type writing

        Parameters
        ----------
        tecplot_filename : str
            the path to the output file
        res_types : str; List[str, str, ...]; default=None -> all
            the results that will be written (must be consistent with
            self.variables)
        adjust_nids : bool; default=True
            element_ids are 0-based in binary and must be switched to
            1-based in ASCII

        """
        self.log.info('writing tecplot %s' % tecplot_filename)
        msg, ivars = self._get_write_header(res_types)

        with open(tecplot_filename, 'w') as tecplot_file:
            tecplot_file.write(msg)
            for zone in self.zones:
                nnodes = zone.nnodes
                nelements = zone.nelements
                (is_structured, is_unstructured, is_points, zone_type,
                 is_tris, is_quads, is_tets, is_hexas) = zone.determine_element_type()
                #print(is_structured, is_unstructured, is_points, zone_type)
                #print(is_tris, is_quads, is_tets, is_hexas)

                if is_unstructured:
                    zone.write_unstructured_zone(tecplot_file, ivars, is_points, nnodes, nelements, zone_type, self.log,
                                                 is_tris, is_quads, is_tets, is_hexas, adjust_nids=adjust_nids)
                elif is_structured:
                    zone.write_structured_zone(tecplot_file, ivars, self.log, zone.headers_dict, adjust_nids=adjust_nids)
                else:  # pragma: no cover
                    raise RuntimeError('only structured/unstructured')

    #def skin_elements(self):
        #sss
        #return tris, quads

    #def get_free_faces(self):
        #"""get the free faces for hexa elements"""
        #sss
        #return free_faces

    def extract_y_slice(self, y0, tol=0.01, slice_filename=None):
        """
        doesn't work...
        """
        self.log.info('slicing...')
        zone = model.zones[0]
        y = self.xyz[:, 1]
        nodes = self.xyz
        assert tol > 0.0, tol
        elements = zone.hexa_elements
        results = zone.nodal_results

        iy = np.where((y0 - tol <= y) & (y <= y0 + tol))[0]

        self.log.debug(y[iy])
        self.log.debug(nodes[iy, 1].min(), nodes[iy, 1].max())
        #iy = np.where(y <= y0 + tol)[0]
        assert len(iy) > 0, iy
        #inode = iy + 1


        # find all elements that have iy within tolerance
        #slots = np.where(elements == iy)
        #slots = np.where(element for element in elements
                          #if any(iy in element))
        #slots = where(iy == elements.ravel())[0]
        ielements = np.unique([ie for ie, unused_elem in enumerate(elements)
                               for i in range(8)
                               if i in iy])
        #print(slots)
        #ri, ci = slots
        #ri = unique(hstack([where(element == iy)[0] for element in elements]))
        #ri = [ie for ie, element in enumerate(elements)
              #if [n for n in element
                  #if n in iy]]
        #ri = [np.where(element == iy)[0] for element in elements if np.where(element == iy)[0]]
        #print(ri)
        #ielements = np.unique(ri)
        self.log.debug(ielements)
        assert len(ielements) > 0, ielements

        # find nodes
        elements2 = elements[ielements, :]
        inodes = np.unique(elements2)
        assert len(inodes) > 0, inodes

        # renumber the nodes
        nidmap = {}
        for inode, nid in enumerate(inodes):
            nidmap[nid] = inode
        elements3 = np.array(
            [[nidmap[nid] for nid in element]
             for element in elements2],
            dtype='int32')

        self.log.debug(inodes)
        nodes2 = nodes[inodes, :]
        nodal_results2 = results[inodes, :]
        model = Tecplot()
        zone = Zone(self.log)
        zone.xyz = nodes2
        zone.nodal_results = nodal_results2
        zone.hexa_elements = elements3
        model.zones = [zone]

        if slice_filename:
            model.write_tecplot(slice_filename)
        return model


def split_headers(header_in):
    #allowed_keys = ['TITLE', 'VARIABLES', 'T', 'ZONETYPE', 'DATAPACKING',
                    #'N', 'E', 'F', 'DT', 'SOLUTIONTIME', 'STRANDID',
                    #'I', 'J', 'K'
                    #]
    #print(f'header1 = {header_in}')
    header = header_in.replace('""', '","')
    #print(f'header2 = {header}')
    cheaders = header.split(',')

    #print(header)
    #print(cheaders)
    #header = cheaders[0]
    #headers = [header]
    #i = 1
    #while i < len(cheaders):
        #headeri = cheaders[i]
        #uheaderi = headeri.upper().replace(' ', '')
        #is_key = [uheaderi.startswith(key+'=') for key in allowed_keys]
        #if any(is_key):
            #print('key!', headeri)
            #header = headeri
            #headers.append(header.lstrip())
        #else:
            #headers[-1] += ',' + headeri
        #i += 1
    #print('headers')
    #for headeri in headers:
        #print('  ', headeri)

    #print(headers)
    #print(header.replace('""', '","'))
    #if ''
    #headers = header.replace('""', '","').split(',')
    return cheaders

def _header_lines_to_header_dict(title_line: str, header_lines: List[str], variables: List[str]):
    """parses the parsed header lines"""
    #print('header_lines', header_lines)
    #headers_dict = {}
    headers_dict = CaseInsensitiveDict()
    if title_line:
        title_sline = title_line.split('=', 1)
        title = title_sline[1]
    else:
        title = 'tecplot geometry and solution file'
    headers_dict['TITLE'] = title

    if len(header_lines) == 0:
        #raise RuntimeError(header_lines)
        return None
    header = ','.join(header_lines)

    # this is so overly complicataed and probably not even enough...
    # what about the following 'quote' style?
    headers = split_headers(header)
    #headers = header.replace('""', '","').split(',')

    #TITLE = "Weights=1/6,6,1"
    #Variables = "x","y","z","psi"
    #Zone N = 125, E = 64, DATAPACKING = POINT, ZONETYPE = FEBRICK

    nheaders = len(headers) - 1
    for iheader, header in enumerate(headers):
        header = header.strip()
        #print(f'{iheader} {header!r}')

    for iheader, header in enumerate(headers):
        header = header.strip()
        #print('%2i %s' % (iheader, header))
        #print('iheader=%s header=%r' % (iheader, header))
        if '=' in header:
            sline = header.split('=', 1)
            parse = False
            #print('iheader=%s nheaders=%s' % (iheader, nheaders))
            if iheader == nheaders:
                parse = True
            elif '=' in headers[iheader + 1]:
                parse = True
        elif header.upper() == 'ZONE':
            # apparently the null key is also a thing...
            # we'll use 'ZONE' because...
            headers_dict['ZONE'] = None
            parse = True
            #continue
        elif '"' in header:
            sline += [header]
            parse = False
            if iheader == nheaders:
                parse = True
            elif '=' in headers[iheader + 1]:
                parse = True
        else:
            raise NotImplementedError('header=%r headers=%r' % (header, headers))

        if parse:
            #print('  parsing')
            #print('sline =', sline)
            key = sline[0].strip().upper()
            if key.startswith('ZONE '):
                # the key is not "ZONE T" or "ZONE E"
                # ZONE is a flag, T is title, E is number of elements
                key = key[5:].strip()

            value = [val.strip() for val in sline[1:]]
            if len(value) == 1:
                value = value[0].strip()
            #assert not isinstance(value, list), value
            headers_dict[key] = value
            #print('  ', value)
            #value = value.strip()

            # 'T', 'ZONE T',  ???
            #   'DT', 'SOLUTIONTIME', 'STRANDID', # tecplot 360 specific things not supported
            allowed_keys = ['VARIABLES', 'T', 'ZONETYPE', 'DATAPACKING', # 'TITLE',
                            'N', 'E', 'F', 'DT', 'SOLUTIONTIME', 'STRANDID',
                            'I', 'J', 'K']
            assert key in allowed_keys, 'key=%r; allowed=[%s]' % (key, ', '.join(allowed_keys))
            parse = False
    #print('headers_dict', headers_dict)
    #print(headers_dict.keys())

    _simplify_header(headers_dict, variables)
    assert len(headers_dict) > 0, headers_dict
    return headers_dict

def _simplify_header(headers_dict, variables: List[str]) -> None:
    """cast the integer headers adn sets the variables"""
    # unstructured
    if 'N' in headers_dict: # nnodes
        headers_dict['N'] = int(headers_dict['N'])
    if 'E' in headers_dict: # nelements
        headers_dict['E'] = int(headers_dict['E'])

    # structured
    if 'I' in headers_dict:
        headers_dict['I'] = int(headers_dict['I'])
    if 'J' in headers_dict:
        headers_dict['J'] = int(headers_dict['J'])
    if 'K' in headers_dict:
        headers_dict['K'] = int(headers_dict['K'])

    #print('_simplify_header', variables, headers_dict)
    if 'TITLE' not  in headers_dict:
        headers_dict['TITLE'] = 'tecplot geometry and solution file'

    if 'VARIABLES' in headers_dict and variables is None:
        #print('VARIABLES' in headers_dict, variables is None)
        _simplify_variables(headers_dict)
    elif 'VARIABLES' in headers_dict:
        _simplify_variables(headers_dict)
    elif variables is not None:
        headers_dict['VARIABLES'] = variables
    else:
        raise RuntimeError('no variables...')

def _simplify_variables(headers_dict) -> None:
    variables = headers_dict['VARIABLES']
    headers_dict['VARIABLES'] = [var.strip('"') for var in variables]

def _stack(zone, xyz_list, quads_list, tris_list, tets_list, hexas_list, results_list, log):
    """
    elements are read as a list of lines, so we need to stack them
    and cast them while we're at it.
    """
    log.debug('stacking elements')
    if len(hexas_list):
        zone.hexa_elements = np.vstack(hexas_list)
    if len(tets_list):
        zone.tet_elements = np.vstack(tets_list)
    if len(quads_list):
        zone.quad_elements = np.vstack(quads_list)
    if len(tris_list):
        zone.tri_elements = np.vstack(tris_list)

    log.debug('stacking nodes')
    if len(xyz_list) == 1:
        xyz = xyz_list[0]
    else:
        xyz = np.vstack(xyz_list)

    #self.elements = elements - 1
    #print(self.elements)

    if is_3d(zone.headers_dict):
        zone.xyz = xyz

        nresults = len(results_list)
        if nresults == 1:
            results = results_list[0]
        else:
            results = np.vstack(results_list)
        zone.nodal_results = results
    else:
        zone.xy = xyz[:, :2]
        nresults = len(results_list) + 1
        nnodes_temp = xyz.shape[0]
        if nresults == 1:
            zone.nodal_results = xyz[:, 2].reshape(nnodes_temp, 1)
        else:
            inputs = [xyz[:, 2].reshape(nnodes_temp, 1), *results_list]
            zone.nodal_results = np.hstack(inputs)
        del nnodes_temp

    zone.variables = [var for var in zone.variables if var not in ['X', 'Y', 'Z']]

def read_zone_block(lines, iline, xyz, results, nresults, zone_type,
                    sline, nnodes, log):
    """a zone can be structured or unstructred"""
    #print('***', iline, sline)

    # read all data
    #result = sline
    #iresult = len(sline)
    #nresult = len(sline)

    result = []
    iresult = 0
    nresult = 0
    nnodes_max = (3 + nresults) * nnodes
    #print('nnodes_max =', nnodes_max)
    while nresult < nnodes_max: #  changed from iresult to nresult
        #print('zb', iline, sline, len(sline))
        result += sline
        nresult += len(sline)
        if iresult >= nnodes_max:
            log.debug('breaking...')
            #break
        iline, line, sline = get_next_sline(lines, iline)
        if iresult == 0:
            log.debug('zone_type=%s sline=%s' % (zone_type, sline))
        iresult += len(sline)
        #print('len', iresult, nresult, len(result))
    #print(result, len(result))
    for i, value in enumerate(result):
        assert '.' in value, 'i=%i value=%s' % (i, value)
    assert len(result) == nnodes_max, 'len(result)=%s expected=%s' % (len(result), nnodes_max)
    #-----------------

    # pack data
    for ires in range(3 + nresults):
        i0 = ires * nnodes
        i1 = (ires + 1) * nnodes #+ 1
        if len(result[i0:i1]) != nnodes:
            msg = 'ires=%s len=%s nnodes=%s' % (
                ires, len(result[i0:i1]), nnodes)
            raise RuntimeError(msg)
        if ires in [0, 1, 2]:
            log.debug('ires=%s nnodes=%s len(result)=%s' % (ires, nnodes, len(result)))
            xyz[:, ires] = result[i0:i1]
        else:
            results[:, ires - 3] = result[i0:i1]

    # setup
    #iline, line, sline = get_next_sline(lines, iline)
    return iline, line, sline

def read_unstructured_elements(lines, iline, sline, elements, nelements):
    assert '.' not in sline[0], sline

    i = 0
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(iline, i, sline)
        try:
            elements[i, :] = sline
        except IndexError:
            raise RuntimeError('i=%s sline=%s' % (i, str(sline)))
        except ValueError:
            raise RuntimeError('i=%s sline=%s' % (i, str(sline)))

        iline, line, sline = get_next_sline(lines, iline)
        #line = lines.readline()
        #iline += 1
        #sline = line.strip().split()
    return iline, line, sline

def read_point(lines, iline, xyz, results, zone_type, line, sline, nnodes, nvars, log):
    """a POINT grid is a structured grid"""
    log.debug(f'start of POINT (structured); nnodes={nnodes} nvars={nvars} zone_type={zone_type}')
    for inode in range(nnodes):
        iline, sline = get_next_nsline(lines, iline, sline, nvars)
        #print(iline, inode, sline)

        #if inode == 0:
            #log.debug('zone_type=%s sline=%s' %(zone_type, sline))


        if not len(sline[3:]) == len(results[inode, :]):
            msg = 'sline[3:]=%s results[inode, :]=%s' % (sline[:3], results[inode, :])
            raise RuntimeError(msg)

        try:
            xyz[inode, :] = sline[:3]
            results[inode, :] = sline[3:]
        except ValueError:
            msg = 'i=%s line=%r\n' % (inode, line)
            msg += 'sline = %s' % str(sline)
            print(msg)
            raise
        iline, line, sline = get_next_sline(lines, iline)
        #log.debug(sline)
    log.debug('end of POINT')
    return iline, line, sline

def read_block(lines, iline, xyz, results, zone_type, line, sline, nnodes, nvars, log):
    """
    BLOCK format is similar to PLOT3D in that you read all the X values before the Ys,
    Zs, and results.  The alternative format is POINT, which reads them on a per node
    basis.
    """
    log.debug('start of BLOCK')
    #print('nnodes =', nnodes)
    #print('nvars =', nvars)
    ndata = nnodes * nvars
    #print('ndata =', ndata)
    results = []

    while len(results) < ndata:
        sline = split_line(line)
        results += sline
        #print('block:', iline, sline, len(results))
        if len(sline) == 0:
            raise
        iline, line, sline = get_next_sline(lines, iline)
        #log.debug(sline)
    #print(len(results))
    assert len(results) == ndata, 'len(results)=%s expected=%s' % (len(results), ndata)
    log.debug('end of BLOCK')

    #TODO: save results
    raise RuntimeError('not done...save results')
    return iline

def get_next_line(lines, iline):
    """Read the next line from the file.  Handles comments."""
    try:
        line = lines[iline].strip()
    except IndexError:
        line = None
        return iline, line

    iline += 1
    igap = 0
    ngap_max = 10
    while len(line) == 0 or line[0] == '#':
        try:
            line = lines[iline].strip()
        except IndexError:
            line = None
            return iline, line
        iline += 1
        if igap > ngap_max:
            break
        igap += 1
    return iline, line

def split_line(line):
    """splits a comma or space separated line"""
    if ',' in line:
        line2 = line.replace(',', ' ')
        sline = line2.split()
    else:
        sline = line.split()
    return sline

def get_next_sline(lines, iline):
    """Read the next split line from the file.  Handles comments."""
    iline, line = get_next_line(lines, iline)
    if line is None:
        return iline, None, None
    sline = split_line(line)
    return iline, line, sline


def get_next_nsline(lines, iline, sline, nvars):
    #print(iline, sline)
    while len(sline) != nvars:  # long line was split
        #print(sline, nvars)
        iline, line, slinei = get_next_sline(lines, iline)
        #print(iline, line, slinei, nvars)
        assert len(slinei) > 0, slinei
        sline += slinei
        #print(sline, '\n')
        #iline += 1
    assert len(sline) == nvars, 'iline=%i sline=%s nvars=%s' % (iline, sline, nvars)
    return iline, sline

def _read_header_lines(lines, iline, line, log):
    """
    reads a tecplot header

    Examples
    --------
    **Example 1**

    TITLE     = "tecplot geometry and solution file"
    VARIABLES = "x"
    "y"
    "z"
    "rho"
    "u"
    "v"
    "w"
    "p"
    ZONE T="\"processor 1\""
    n=522437, e=1000503, ZONETYPE=FEBrick
    DATAPACKING=BLOCK

    **Example 2**

    title="Force and Momment Data for forces"
    variables="Iteration"
    "C_L","C_D","C_M_x","C_M_y","C_M_z""C_x","C_y","C_z","C_Lp","C_Dp", "C_Lv", "C_Dv""C_M_xp"
    "C_M_yp","C_M_zp","C_M_xv","C_M_yv""C_M_zv","C_xp","C_yp","C_zp","C_xv","C_yv""C_zv
    "Mass flow","<greek>r</greek>","u"
    "p/p<sub>0</sub>","T","p<sub>t</sub>/p<sub>0</sub>"
    "T<sub>t</sub>","Mach"
    "Simulation Time"
    zone,t="forces"
    """
    i = 0
    title_line = ''
    #variables_line = ''
    active_key = None

    vars_found = []
    header_lines = []
    #print('-----------------------------')
    #for iii, linei in enumerate(lines):
        #if iii > 10:
            #break
        #print(linei)
    #print('-----------------------------')
    while i < 30:
        #print(iline, i, line.strip())
        #self.n = 0
        if len(line) == 0 or line[0] == '#':
            line = lines[iline].strip()
            iline += 1
            i += 1
            continue
        if line[0].isdigit() or line[0] == '-':
            #print(line)
            log.debug('breaking...')
            break

        uline = line.upper()
        uline2 = uline.replace(' ', '')
        if 'TITLE=' in uline2:
            title_line += line
            vars_found.append('TITLE')
            active_key = 'TITLE'
        elif 'VARIABLES' in uline2:
            vars_found.append('VARIABLES')
            #variables_line += line
            active_key = 'VARIABLES'
        else:
            #if 'ZONE T' in line:
                #vars_found.append('ZONE T')
            if 'ZONE' in uline2:
                vars_found.append('ZONE')
                active_key = 'ZONE'
            #if 'ZONE N' in uline:
                #vars_found.append('N')
            if 'ZONETYPE' in uline2:
                vars_found.append('ZONETYPE')
                active_key = 'ZONE'
            if 'DATAPACKING' in uline2:
                vars_found.append('DATAPACKING')
                active_key = 'ZONE'

        #print(active_key, line)
        if active_key in ['ZONE', 'VARIABLES']:
            header_lines.append(line.strip())
        #if len(vars_found) == 5:
            #break

        #if active_key
        i += 1
        line = lines[iline].strip()
        iline += 1

    log.debug('vars_found = %s' % vars_found)
    #print('header_lines', header_lines)
    #print("title = %r" % title_line)
    #print("variables_line = %r" % variables_line)

    return iline, title_line, header_lines, line

def main():  # pragma: no cover
    #plt = Tecplot()
    fnames = os.listdir(r'Z:\Temporary_Transfers\steve\output\time20000')

    #datai = [
        #'n=3807, e=7443',
        #'n=3633, e=7106',
        #'n=3847, e=7332',
        #'n=3873, e=6947',
        #'n=4594, e=8131',
        #'n=4341, e=7160',
        #'n=4116, e=8061',
        #'n=4441, e=8105',
        #'n=4141, e=8126',
        #'n=4085, e=8053',
        #'n=4047, e=8215',
        #'n=4143, e=8123',
        #'n=4242, e=7758',
        #'n=3830, e=7535',
        #'n=3847, e=7936',
        #'n=3981, e=7807',
        #'n=3688, e=7415',
        #'n=4222, e=8073',
        #'n=4164, e=7327',
        #'n=3845, e=8354',
        #'n=4037, e=6786',
        #'n=3941, e=8942',
        #'n=4069, e=7345',
        #'n=4443, e=8001',
        #'n=3895, e=7459',
        #'n=4145, e=7754',
        #'n=4224, e=8152',
        #'n=4172, e=7878',
        #'n=4138, e=8864',
        #'n=3801, e=7431',
        #'n=3984, e=6992',
        #'n=4195, e=7967',
        #'n=4132, e=7992',
        #'n=4259, e=7396',
        #'n=4118, e=7520',
        #'n=4176, e=7933',
        #'n=4047, e=8098',
        #'n=4064, e=8540',
        #'n=4144, e=8402',
        #'n=4144, e=7979',
        #'n=3991, e=6984',
        #'n=4080, e=8465',
        #'n=3900, e=7981',
        #'n=3709, e=8838',
        #'n=4693, e=8055',
        #'n=4022, e=7240',
        #'n=4028, e=8227',
        #'n=3780, e=7551',
        #'n=3993, e=8671',
        #'n=4241, e=7277',
        #'n=4084, e=6495',
        #'n=4103, e=8165',
        #'n=4496, e=5967',
        #'n=3548, e=8561',
        #'n=4143, e=7749',
        #'n=4136, e=8358',
        #'n=4096, e=7319',
        #'n=4209, e=8036',
        #'n=3885, e=7814',
        #'n=3800, e=8232',
        #'n=3841, e=7837',
        #'n=3874, e=7571',
        #'n=3887, e=8079',
        #'n=3980, e=7834',
        #'n=3763, e=7039',
        #'n=4287, e=7130',
        #'n=4110, e=8336',
        #'n=3958, e=7195',
        #'n=4730, e=7628',
        #'n=4087, e=8149',
        #'n=4045, e=8561',
        #'n=3960, e=7320',
        #'n=3901, e=8286',
        #'n=4065, e=7013',
        #'n=4160, e=7906',
        #'n=3628, e=7140',
        #'n=4256, e=8168',
        #'n=3972, e=8296',
        #'n=3661, e=7879',
        #'n=3922, e=8093',
        #'n=3972, e=6997',
        #'n=3884, e=7603',
        #'n=3609, e=6856',
        #'n=4168, e=7147',
        #'n=4206, e=8232',
        #'n=4631, e=8222',
        #'n=3970, e=7569',
        #'n=3998, e=7617',
        #'n=3855, e=7971',
        #'n=4092, e=7486',
        #'n=4407, e=7847',
        #'n=3976, e=7627',
        #'n=3911, e=8483',
        #'n=4144, e=7919',
        #'n=4033, e=8129',
        #'n=3976, e=7495',
        #'n=3912, e=7739',
        #'n=4278, e=8522',
        #'n=4703, e=8186',
        #'n=4230, e=7811',
        #'n=3971, e=7699',
        #'n=4081, e=8242',
        #'n=4045, e=7524',
        #'n=4532, e=5728',
        #'n=4299, e=8560',
        #'n=3885, e=7531',
        #'n=4452, e=8405',
        #'n=4090, e=7661',
        #'n=3937, e=7739',
        #'n=4336, e=7612',
        #'n=4101, e=7461',
        #'n=3980, e=8632',
        #'n=4523, e=7761',
        #'n=4237, e=8463',
        #'n=4013, e=7856',
        #'n=4219, e=8013',
        #'n=4248, e=8328',
        #'n=4529, e=8757',
        #'n=4109, e=7496',
        #'n=3969, e=8026',
        #'n=4093, e=8506',
        #'n=3635, e=7965',
        #'n=4347, e=8123',
        #'n=4703, e=7752',
        #'n=3867, e=8124',
        #'n=3930, e=7919',
        #'n=4247, e=7154',
        #'n=4065, e=8125',
    #]
    fnames = [os.path.join(r'Z:\Temporary_Transfers\steve\output\time20000', fname)
              for fname in fnames]
    tecplot_filename_out = None
    #tecplot_filename_out = 'tecplot_joined.plt'
    from pyNastran.converters.tecplot.utils import merge_tecplot_files
    model = merge_tecplot_files(fnames, tecplot_filename_out)

    y0 = 0.0
    model.extract_y_slice(y0, tol=0.014, slice_filename='slice.plt')

    return
    #for iprocessor, fname in enumerate(fnames):
        #nnodes, nelements = datai[iprocessor].split(',')
        #nnodes = int(nnodes.split('=')[1])
        #nelements = int(nelements.split('=')[1])

        #ip = iprocessor + 1
        #tecplot_filename = 'model_final_meters_part%i_tec_volume_timestep20000.plt' % ip
        #print(tecplot_filename)
        #try:
            #plt.read_tecplot_binary(tecplot_filename, nnodes=nnodes, nelements=nelements)
            #plt.write_tecplot('processor%i.plt' % ip)
        #except:
            #raise
        ##break

def main2():  # pragma: no cover
    """tests slicing"""
    plt = Tecplot()
    #fnames = os.listdir(r'Z:\Temporary_Transfers\steve\output\time20000')
    #fnames = [os.path.join(r'Z:\Temporary_Transfers\steve\output\time20000', fname)
    #          for fname in fnames]
    fnames = ['slice.plt']
    # tecplot_filename_out = None
    #tecplot_filename_out = 'tecplot_joined.plt'
    #model = merge_tecplot_files(fnames, tecplot_filename_out)

    for iprocessor, tecplot_filename in enumerate(fnames):
        plt.read_tecplot(tecplot_filename)
        plt.write_tecplot('processor_%i.plt' % iprocessor)

if __name__ == '__main__':   # pragma: no cover
    main()
