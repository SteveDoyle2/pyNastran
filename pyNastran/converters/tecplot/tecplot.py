# coding: utf-8
"""
models from:
    http://people.sc.fsu.edu/~jburkardt/data/tec/tec.html
"""
from __future__ import print_function
import sys
import os
from struct import unpack
from collections import defaultdict
import itertools

from six import string_types
import numpy as np
from numpy import (
    array, vstack, hstack, where, unique, zeros, loadtxt, savetxt, intersect1d, in1d)
#import numpy as np
from cpylog import get_logger2

from pyNastran.utils import is_binary_file


def read_tecplot(tecplot_filename, use_cols=None, dtype=None, log=None, debug=False):
    """loads a tecplot file"""
    tecplot = Tecplot(log=log, debug=debug)
    if use_cols:
        tecplot.use_cols = use_cols
        tecplot.dtype = dtype
    tecplot.read_tecplot(tecplot_filename)
    return tecplot


class Tecplot(object):
    """
    Parses a hexa binary/ASCII Tecplot 360 file.
    Writes an ASCII Tecplot 10 file (no transient support).
    """
    def __init__(self, log=None, debug=False):
        # defines binary file specific features
        self._endian = b'<'
        self._n = 0

        self.tecplot_filename = ''
        self.log = get_logger2(log, debug=debug)
        self.debug = debug

        # mesh = None : model hasn't been read
        self.is_mesh = None

        # mesh = True : this is a structured/unstructured grid
        self.xyz = array([], dtype='float32')
        self.tet_elements = array([], dtype='int32')
        self.hexa_elements = array([], dtype='int32')
        self.quad_elements = array([], dtype='int32')
        self.tri_elements = array([], dtype='int32')
        self.results = array([], dtype='float32')
        self.variables = []

        # mesh = False : this is a plot file
        self.use_cols = None
        self.dtype = None
        self.A = None

        self._uendian = ''
        self.n = 0


    @property
    def result_names(self):
        """gets the variables"""
        return self.variables

    @result_names.setter
    def result_names(self, vals):
        """sets the variables"""
        self.variables = vals

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

    def read_header_lines(self, tecplot_file, line):
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
        vars_found = []
        header_lines = []
        while i < 30:
            #print(i, line.strip())
            #self.n = 0
            if len(line) == 0 or line[0] == '#':
                line = tecplot_file.readline().strip()
                i += 1
                continue
            if line[0].isdigit() or line[0] == '-':
                self.log.debug('breaking...')
                break
            if 'TITLE' in line:
                vars_found.append('TITLE')
            if 'VARIABLES' in line:
                vars_found.append('VARIABLES')
            if 'ZONE T' in line:
                vars_found.append('ZONE T')
            if 'ZONE' in line:
                vars_found.append('ZONE')
            #if 'ZONE N' in line:
                #vars_found.append('N')
            if 'ZONETYPE' in line:
                vars_found.append('ZONETYPE')
            if 'DATAPACKING' in line:
                vars_found.append('DATAPACKING')
            header_lines.append(line.strip())
            #if len(vars_found) == 5:
                #break
            i += 1
            line = tecplot_file.readline().strip()

        self.log.debug('vars_found = %s' % vars_found)

        return header_lines, i, line

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

        nnodes = -1
        nelements = -1
        with open(tecplot_filename, 'r') as tecplot_file:
            quads_list = []
            hexas_list = []
            tris_list = []
            tets_list = []

            xyz_list = []
            results_list = []

            line = tecplot_file.readline().strip()
            iblock = 0
            while 1:
                header_lines, unused_i, line = self.read_header_lines(tecplot_file, line)
                #print(header_lines)
                headers_dict = _header_lines_to_header_dict(header_lines)
                if headers_dict is None:
                    break

                #print(headers_dict.keys())
                if 'ZONETYPE' in headers_dict:
                    zone_type = headers_dict['ZONETYPE'].upper() # FEBrick
                    data_packing = headers_dict['DATAPACKING'].upper() # block
                    self._read_zonetype(zone_type, tecplot_file, iblock, headers_dict, line,
                                        nnodes, nelements,
                                        xyz_list, hexas_list, tets_list, quads_list, tris_list,
                                        results_list,
                                        data_packing=data_packing)
                elif 'F' in headers_dict:
                    fe = headers_dict['F'] # FEPoint
                    assert isinstance(fe, str), headers_dict
                    zone_type = fe.upper() # FEPoint
                    assert zone_type == 'FEPOINT', zone_type
                    self.log.debug('zone_type = %r' % zone_type[0])
                    self._read_zonetype(zone_type, tecplot_file, iblock, headers_dict, line,
                                        nnodes, nelements,
                                        xyz_list, hexas_list, tets_list, quads_list, tris_list,
                                        results_list,
                                        fe=fe)
                elif (('ZONE' in headers_dict) and
                      (headers_dict['ZONE'] is None) and
                      ('T' in headers_dict)):
                    A, line = self.read_table(tecplot_file, iblock, headers_dict, line)
                    self.A = A
                    return

                else:
                    msg = 'headers=%s\n' % str(headers_dict)
                    msg += 'line = %r' % line.strip()
                    raise NotImplementedError(msg)

        self.log.debug('stacking elements')
        if len(hexas_list):
            self.hexa_elements = vstack(hexas_list)
        if len(tets_list):
            self.tet_elements = vstack(tets_list)
        if len(quads_list):
            self.quad_elements = vstack(quads_list)
        if len(tris_list):
            self.tri_elements = vstack(tris_list)

        self.log.debug('stacking nodes')
        if len(xyz_list) == 1:
            xyz = xyz_list[0]
        else:
            xyz = vstack(xyz_list)
        if len(results_list) == 1:
            results = results_list[0]
        else:
            results = vstack(results_list)
        #self.elements = elements - 1
        #print(self.elements)

        self.xyz = xyz
        self.results = results

    def read_table(self, tecplot_file, unused_iblock, headers_dict, line):
        """
        reads a space-separated tabular data block
        """
        variables = [var.strip('" ') for var in headers_dict['VARIABLES']]
        print('variables = %s' % variables)
        #self.dtype[]
        use_cols = [variables.index(var) for var in self.use_cols]

        # add on the preceding line to the line "list"
        # that's not a hack at all...
        lines = itertools.chain((line, ), iter(tecplot_file))
        A = loadtxt(lines, dtype=self.dtype, comments='#', delimiter=None,
                    converters=None, skiprows=0,
                    usecols=use_cols, unpack=False, ndmin=0)
        return A, None

    def _read_zonetype(self, zone_type, tecplot_file, iblock, headers_dict, line,
                       nnodes, nelements,
                       xyz_list, hexas_list, tets_list, quads_list, tris_list,
                       results_list,
                       data_packing=None, fe=None):
        """
        reads:
          - ZONE E
          - ZONE T

        ZONE is a flag, T is title, E is number of elements
        """
        if iblock == 0:
            variables = headers_dict['VARIABLES']
            self.variables = [variable.strip(' \r\n\t"\'') for variable in variables]
            self.log.debug('self.variables = %s' % self.variables)
            nresults = len(variables) - 3 # x, y, z, rho, u, v, w, p
            self.log.debug('nresults = %s' % nresults)

        self.log.debug(headers_dict)
        nnodesi = int(headers_dict['N'])
        nelementsi = int(headers_dict['E'])
        self.log.info('nnodes=%s nelements=%s' % (nnodesi, nelementsi))
        xyz = zeros((nnodesi, 3), dtype='float32')
        results = zeros((nnodesi, nresults), dtype='float32')
        if zone_type == 'FEBRICK':
            # hex
            elements = zeros((nelementsi, 8), dtype='int32')
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETETRAHEDRON'):
            # quads / tets
            elements = zeros((nelementsi, 4), dtype='int32')
        elif zone_type == 'FETRIANGLE':
            # tris
            elements = zeros((nelementsi, 3), dtype='int32')
        #elif zone_type == 'FEBLOCK':
            #pass
        else:
            #if isinstance(zone_type, list):
                #raise NotImplementedError(zone_type[0])
            raise NotImplementedError(zone_type)

        sline = line.strip().split()
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

                    line = tecplot_file.readline().strip()
                    while len(line) == 0 or line[0] == '#':
                        line = tecplot_file.readline().strip()
                    sline = line.split()
            elif data_packing == 'BLOCK':
                result = []
                iresult = 0
                nresult = 0
                nnodes_max = (3 + nresults) * nnodesi
                while iresult < nnodes_max:
                    result += sline
                    nresult += len(sline)
                    if iresult >= nnodes_max:
                        self.log.debug('breaking...')
                        #break
                    line = tecplot_file.readline().strip()
                    while line[0] == '#':
                        line = tecplot_file.readline().strip()
                    sline = line.split()
                    if iresult == 0:
                        self.log.debug('zone_type=%s sline=%s' % (zone_type, sline))
                    iresult += len(sline)

                for ires in range(3 + nresults):
                    i0 = ires * nnodesi
                    i1 = (ires + 1) * nnodesi #+ 1
                    if len(result[i0:i1]) != nnodesi:
                        msg = 'ires=%s len=%s nnodesi=%s' % (
                            ires, len(result[i0:i1]), nnodesi)
                        raise RuntimeError(msg)
                    if ires in [0, 1, 2]:
                        self.log.debug('ires=%s nnodesi=%s len(result)=%s' %
                                       (ires, nnodesi, len(result)))
                        xyz[:, ires] = result[i0:i1]
                    else:
                        results[:, ires - 3] = result[i0:i1]
            else:
                raise NotImplementedError(data_packing)
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETRIANGLE'):
            sline = line.strip().split()
            for inode in range(nnodesi):
                #print('inode ', inode, sline)
                xyz[inode, :] = sline[:3]
                #if abs(xyz[inode, 1]) <= 5.0:
                    #msg = 'inode=%s xyz=%s'  % (inode, xyz[inode, :])
                    #raise RuntimeError(msg)

                results[inode, :] = sline[3:]
                line = tecplot_file.readline().strip()
                while line[0] == '#':
                    line = tecplot_file.readline().strip()
                sline = line.split()
        else:
            raise NotImplementedError(zone_type)

        i = 0
        #print(elements.shape)
        #print('xyz[0 , :]', xyz[0, :])
        #print('xyz[-1, :]', xyz[-1, :])
        #print(sline)
        assert '.' not in sline[0], sline

        i = 0
        for i in range(nelementsi):
            try:
                elements[i, :] = sline
            except IndexError:
                raise RuntimeError('i=%s sline=%s' % (i, str(sline)))
            except ValueError:
                raise RuntimeError('i=%s sline=%s' % (i, str(sline)))
            line = tecplot_file.readline()
            sline = line.strip().split()
        #print(f.readline())

        if zone_type == 'FEBRICK':
            hexas_list.append(elements + nnodes)
        elif zone_type == 'FETETRAHEDRON':
            tets_list.append(elements + nnodes)
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL'):
            quads_list.append(elements + nnodes)
        elif zone_type == 'FETRIANGLE':
            tris_list.append(elements + nnodes)
        else:
            raise NotImplementedError(zone_type)
        xyz_list.append(xyz)
        results_list.append(results)
        nnodes += nnodesi
        nelements += nelementsi
        self.log.debug('nnodes=%s nelements=%s' % (nnodes, nelements))
        del headers_dict
        iblock += 1
        if iblock == 10:
            return
        self.log.debug('final sline=%s' % sline)

    @property
    def nnodes(self):
        """gets the number of nodes in the model"""
        return self.xyz.shape[0]

    @property
    def nelements(self):
        """gets the number of elements in the model"""
        return (self.hexa_elements.shape[0] + self.tet_elements.shape[0] +
                self.quad_elements.shape[0] + self.tri_elements.shape[0])

    def read_word(self, tecplot_file):
        n0 = self.n
        n = 0
        while 1:
            data = tecplot_file.read(4)
            value, = unpack(b'f', data)
            n += 4
            if value == 299.0:  # 299.0 is the zone_marker
                break
        tecplot_file.seek(n0)
        data = tecplot_file.read(n)
        self.show_data(data)
        self.n += n

        nbytes_title = n - 4
        n_title_letters = (n - 4) // 4
        int_letters = unpack('i'*n_title_letters, data[:nbytes_title])
        word = ''
        words = []
        for i, letter in enumerate(int_letters):
            if letter == 0:
                words.append(word)
                word = ''
                break

            char = chr(letter)
            word += char
            # letters = [chr(letter) for letter in int_letters]
        self.n -= nbytes_title
        #self.n += len(words[0])
        tecplot_file.seek(self.n)
        return words[0]

    def read_var(self, tecplot_file):
        var = ''
        cont = True
        while cont:
            data = tecplot_file.read(4)
            self.n += 4
            chari, = unpack(b'i', data)
            if chari != 0:
                var += chr(data[0])
            else:
                cont = False
        return var

    def read_tecplot_binary(self, tecplot_filename, nnodes=None,
                            nelements=None):
        """
        The binary file reader must have ONLY CHEXAs and be Tecplot 360
        with:
        `rho`, `u`, `v`, `w`, and `p`.

        """
        self.tecplot_filename = tecplot_filename
        assert os.path.exists(tecplot_filename), tecplot_filename

        zone_type_map = {
            0 : 'ORDERED',
            1 : 'FELINESEG',
            2 : 'FETRIANGLE',
            3 : 'FEQUADRILATERAL',
            4 : 'FETETRAHEDRON',
            5 : 'FEBRICK',
            6 : 'FEPOLYGON',
            7 : 'FEPOLYHEDRON',
        }
        data_packing_map = {
            0 : 'BLOCK',
            1 : 'POINT',
        }
        var_location_map = {
            0 : 'NODE',
            1 : 'CELL-CENTERED',
        }

        with open(tecplot_filename, 'rb') as tecplot_file:
            self.f = tecplot_file
            self._uendian = '<'
            self.n = 0
            self.variables = ['rho', 'u', 'v', 'w', 'p']

            # http://home.ustc.edu.cn/~cbq/360_data_format_guide.pdf
            # ------------------------------------------
            # header section

            #  i. Magic number, Version number
            data = tecplot_file.read(8)
            #self.show_data(data)
            self.n += 8
            word, = unpack(b'8s', data)
            word = word.decode('utf8')
            try:
                magic_number, version = word.split('V')
            except:
                print('word = %r' %  word)
                raise
            self.log.debug('word=%r magic_number=%r version=%r' % (word, magic_number, version))

            # ii. Integer value of 1.
            data = tecplot_file.read(4)
            self.n += 4
            byte_order = unpack(b'i', data)
            self.log.debug('byte_order=%r' % (byte_order))

            # iii. Title and variable names.
            # FileType:
            #    0 = FULL
            #    1 = GRID
            #    2 = SOLUTION
            data = tecplot_file.read(4)
            self.n += 4
            file_type = unpack(b'i', data)
            self.log.debug('file_type=%r (0=full, 1=grid, 2=solution)' % (file_type))

            # Title (INT32*N)
            title = self.read_var(tecplot_file)
            self.log.debug('TITLE = "%s"' % title)

            # Number of variables (INT32)
            data = tecplot_file.read(4)
            self.n += 4
            num_vars, = unpack(b'i', data)
            self.log.debug('Number of variables = %d' % num_vars)

            # Variable names
            variable_names = []
            for i in range(num_vars):
                var = self.read_var(tecplot_file)
                variable_names.append(var.lower())
            self.log.debug('Variables: %s' % variable_names)
            #self.show(164)

            #header = ''.join(letters)
            #print("header='%s'" % header, len(header))
            #self.show(200, types='i', endian='<')

            #  iv. Zones (p.152)
            zones = []
            zone_marker = 0.0
            izone = 0
            while zone_marker != 357.0:
                print('------------------------------------------------------------')
                # Zone marker (FLOAT32) (should be 299.0)
                data = tecplot_file.read(4)
                self.n += 4

                zone_marker, = unpack(b'f', data)
                self.log.debug('zone_marker = %s' %  zone_marker)
                if zone_marker == 357.0:
                    break
                #assert zone_marker == 0.0, zone_marker

                #--------------------------------------------------------------------
                # you need to loop over everything until you get 357.0 (end of header)
                # If you get more that don't work passed this point, you may have to
                # add:
                #       399.0   (Geometry marker)
                #       499.0   (Text)
                #       599.0   (CustomLabel)
                #       699.0   (UserRec)
                #       899.0   (Variable Auxiliary data)
                #--------------------------------------------------------------------
                if zone_marker == 299.0:
                    zone = {}
                    # Zone name (INT32*N)
                    cont = True
                    zone_name = ''
                    zone_name = self.read_var(tecplot_file)
                    self.log.debug('Zone name = "%s"' % zone_name)
                    zone['NAME'] = zone_name

                    if version == '102':
                        strand_id = None
                        solution_time = None
                        data = tecplot_file.read(16)
                        self.n += 16
                        # i            i            i                 i
                        parent_zone, zone_type_int, data_packing_int, var_location_int = unpack(b'<4i', data)
                        self.log.debug('parent_zone=%s, zone_type=%s' % (
                            parent_zone, zone_type_int))
                    else:
                        data = tecplot_file.read(32)
                        self.n += 32
                        # i             i            d            i        i            i                 i
                        parent_zone, strand_id, solution_time, not_used, zone_type_int, data_packing_int, var_location_int = unpack(b'<2i d iiii', data)
                        self.log.debug('parent_zone=%s,  strand_id=%s, solution_time=%s, not_used=%s, zone_type=%s' % (
                            parent_zone,  strand_id, solution_time, not_used, zone_type_int))

                    # ParentZone: Zero-based zone number within this
                    #             datafile to which this zone is
                    #             a child.

                    # StrandID: -2 = pending strand ID for assignment by Tecplot
                    #           -1 = static strand ID
                    #            0 <= N < 32700 valid strand ID

                    zone_type = zone_type_map[zone_type_int]
                    self.log.debug('zone_type=%s' % zone_type)

                    # Data Packing (INT32)
                    data_packing = data_packing_map[data_packing_int]
                    self.log.debug('Data packing = %s' % data_packing)

                    # Specify Var Location (INT32)
                    #specify_var_location = 'UNKNOWN'
                    if var_location_int == 0:
                        specify_var_location = "Don't specify"
                    elif var_location_int == 1:
                        specify_var_location = 'Specify'
                        # Variable Locations (INT32*NV)
                        # 0 = Node
                        # 1 = Cell Centered
                        nbytes = num_vars * 4
                        data = tecplot_file.read(nbytes)
                        self.n += 4
                        assert nbytes == 4, nbytes
                        var_location_int, = unpack(b'i', data)
                        var_location = var_location_map[var_location_int]
                        #self.log.debug('var_location_int = %d' % var_location_int)
                        self.log.debug('var location = %s' % var_location)
                    else:
                        raise NotImplementedError(var_location_int)
                    self.log.debug('Specify var location = %s' % specify_var_location)

                    # Are raw local 1-to-1 face neighbors supplied?
                    # (0=FALSE 1=TRUE).
                    #
                    # These raw values are a compact form of the local
                    # 1-to-1 face neighbors.  If supplied, Tecplot
                    # assumes that the face neighbors are fully specified.
                    # As such, it will not perform auto face neighbor
                    # assignment.

                    # This improves Tecplot’s time to first plot.
                    # See the data section below for format details.
                    # ORDERED and FELINESEG zones must specify 0 for
                    # this value because raw face neighbors are not
                    # defined for these zone types. FEPOLYGON and
                    # FEPOLYHEDRON zones must specify 0 for this value
                    # since face neighbors are defined in the face map
                    # for these zone types.
                    #
                    # Raw local 1-to-1 face neighbors (INT32)
                    data = tecplot_file.read(4)
                    self.n += 4
                    are_raw_local_one_to_one_face_neighbors, = unpack(b'i', data)
                    self.log.debug('are_raw_local_one_to_one_face_neighbors=%s' % (are_raw_local_one_to_one_face_neighbors))

                    num_misc_user_defned_face_neighbor_connections_int = None
                    if are_raw_local_one_to_one_face_neighbors == 1:
                        # Number of misc user-defined face neighbor connection (INT32)
                        data = tecplot_file.read(4)
                        self.n += 4
                        num_misc_user_defned_face_neighbor_connections_int, = unpack(b'i', data)

                        num_misc_user_defned_face_neighbor_connections_map = {
                            0 : 'Local 1-to-1',
                            1 : 'Local 1-to-many',
                            2 : 'Global 1-to-1',
                            3 : 'Global 1-to-many',
                        }
                        misc_user_defned_face_neighbor_connections =  num_misc_user_defned_face_neighbor_connections_map[
                            num_misc_user_defned_face_neighbor_connections_int]
                        self.log.debug('misc_user_defned_face_neighbor_connections=%s' % (
                            misc_user_defned_face_neighbor_connections))
                        if num_misc_user_defned_face_neighbor_connections_int != 0:
                            # User-defined face neighbor mode (INT32)
                            data = tecplot_file.read(4)
                            self.n += 4
                            user_defined_face_neighbor_mode, = unpack(b'i', data)
                            self.log.debug('user_defined_face_neighbor_mode=%s' % (
                                user_defined_face_neighbor_mode))
                            if  'FE' in zone_type:
                                # If FE face neighbors are specified by misc face neighbors given (INT32)
                                data = tecplot_file.read(4)
                                fe_face_neighbors_specified_by_misc_face_neighbors_given, = unpack(b'i', data)
                                self.log.debug('fe_face_neighbors_specified_by_misc_face_neighbors_given=%s' % (
                                    fe_face_neighbors_specified_by_misc_face_neighbors_given))

                    if zone_type == 'ORDERED':
                        # IMax, JMax, KMax (INT32*3)
                        data = tecplot_file.read(3 * 4)
                        self.n += 12
                        imax, jmax, kmax, = unpack(b'3i', data)
                        self.log.debug('IMax, JMax, KMax = %d, %d, %d' % (imax, jmax, kmax))
                        aasdf

                    elif 'FE' in zone_type:  # FE-ZONE
                        # NumPts (INT32)
                        data = tecplot_file.read(4)
                        self.n += 4
                        num_pts, = unpack(b'i', data)
                        self.log.debug('Num points = %d' % (num_pts))
                        zone['NPNTS'] = num_pts

                        #if FE Zone:
                          #+-------+--------+
                          #| INT32 | NumPts |
                          #+-------+--------+
                          #if ZoneType is FEPOLYGON or FEPOLYHEDRON:

                        if zone_type == 'FEPOLYGON' or zone_type == 'FEPOLYHEDRON':
                            # NumFaces (INT32)
                            data = tecplot_file.read(4)
                            self.n += 4
                            num_faces, = unpack(b'i', data)
                            self.log.debug('Num faces = %d' % (num_faces))

                            # Total number of face nodes (INT32)
                            data = tecplot_file.read(4)
                            self.n += 4
                            num_face_nodes, = unpack(b'i', data)
                            self.log.debug('Num face nodes = %d' % (num_face_nodes))

                            # Total number of boundary faces (INT32)
                            data = tecplot_file.read(4)
                            self.n += 4
                            num_bnd_faces, = unpack(b'i', data)
                            self.log.debug('Num boundary faces = %d' % (num_bnd_faces))

                            # Total number of boundary connections (INT32)
                            data = tecplot_file.read(4)
                            self.n += 4
                            num_bnd_con, = unpack(b'i', data)
                            self.log.debug('Num boundary connections = %d' % (num_bnd_con))

                        #For all zone types (repeat for each Auxiliary data name/value pair)
                        #+-----------+
                        #|   INT32   | 1=Auxiliary name/value pair to follow
                        #+-----------+ 0=No more Auxiliary name/value pairs.
                        #
                        #If the above is 1, then supply the following:
                        #+-----------+
                        #|  INT32*N  | name string (See note 1.)
                        #+-----------+
                        #+-----------+
                        #|   INT32   | Auxiliary Value Format
                        #+-----------+ (Currently only allow 0=AuxDataType_String)
                        #
                        #+-----------+
                        #|  INT32*N  | Value string (See note 1.)
                        #+-----------+
                        #raise NotImplementedError(zone_type)
                        #
                        # NumElements (INT32)
                        data = tecplot_file.read(4)
                        self.n += 4
                        num_elems, = unpack(b'i', data)
                        self.log.debug('Num elements = %d' % (num_elems))
                        zone['NELEMS'] = num_elems

                        # ICellDim, JCellDim, KCellDim (INT32*3)
                        data = tecplot_file.read(3 * 4)
                        self.n += 12
                        icelldim, jcelldim, kcelldim, = unpack(b'3i', data)
                        self.log.debug('ICellDim, JCellDim, KCellDim = %d, %d, %d' % (icelldim, jcelldim, kcelldim))

                        # Aux continue (INT32)
                        data = tecplot_file.read(4)
                        self.n += 4
                        aux_continue, = unpack(b'i', data)
                        self.log.debug('Aux continue = %d' % (aux_continue))
                        if aux_continue == 1:
                            #TODO
                            # Name string
                            # Aux Value Format
                            # Value string
                            pass
                            raise NotImplementedError(aux_continue)
                        else:
                            # TODO: why???
                            # this has been 36 and 48 (9/12 in words)
                            datai = b''
                            ni = 0
                            self.show(4)
                            while ni < 1000:
                                dataii = tecplot_file.read(4)
                                floati, = unpack('f', dataii)
                                if floati in [299., 357.]:
                                    print('***floati', floati)
                                    break
                                if len(dataii) == 0:
                                    raise RuntimeError('bad data')
                                datai += dataii
                                ni += 4
                            if ni == 1000:
                                raise RuntimeError('didnt find end of table')
                            self.n += ni
                            tecplot_file.seek(self.n)
                            self.show(ni)
                            #n0 = self.n
                            #aux_word = self.read_var(tecplot_file)
                            #print('aux_word =', aux_word)
                            #print('dn =', self.n - n0)
                            #self.show(56)
                            #data = tecplot_file.read(36)
                            #self.show_data(data, types='ifs', endian=None)
                            #self.n += 36
                            #raise NotImplementedError(aux_continue)
                        zones.append(zone)
                        print(zone)
                    else:
                        raise NotImplementedError(zone_type_int)
                    #sss
                elif zone_marker in [399., 499., 599., 699., 799., 899.]:
                    raise NotImplementedError('zone_marker = %r' % zone_marker)
                else:
                    raise NotImplementedError('zone_marker = %r' % zone_marker)

                #self.show(40)
                #aaa


            pshell = 1

            # ------------------------------------------------------------------------------------------------------------
            variable_data_format_map = {
                1 : 'float',
                2 : 'double',
                3 : 'long',
                4 : 'short',
                5 : 'byte',
                6 : 'bit',
            }
            #II. DATA SECTION (don’t forget to separate the header from the data
            #                  with an EOHMARKER). The data section contains all
            #                  of the data associated with the zone definitions
            #                  in the header.
            for i, zone in enumerate(zones):
                zone_name = zone['NAME']
                npoints = zone['NPNTS']
                nelements = zone['NELEMS']
                self.log.info('Processing Zone: %s, npts: %d, nelems: %d' % (
                    zone_name, npoints, nelements))
                variables = {}

                # i. For both ordered and fe zones:
                #    Zone marker (FLOAT32); p.
                data = tecplot_file.read(4)
                self.n += 4
                zone_marker, = unpack(b'f', data)
                self.log.debug('  Zone marker: %f' % zone_marker)

                # Variable data format (INT32*N)
                # Variable data format, N=Total number of vars
                #  1=Float, 2=Double, 3=LongInt,
                #  4=ShortInt, 5=Byte, 6=Bit
                data = tecplot_file.read(4 * num_vars)
                self.n += 4 * num_vars
                variable_data_format_ints = np.frombuffer(data, dtype=np.int32)
                variable_data_format = [variable_data_format_map[variable_data_format_int]
                                        for variable_data_format_int in variable_data_format_ints]
                self.log.debug('  Variable data format: %r' % variable_data_format_ints)
                self.log.debug('  Variable data format: %r' % variable_data_format)

                # Has passive variables (INT32)
                data = tecplot_file.read(4)
                self.n += 4
                npassive_vars, = unpack(b'i', data)
                self.log.debug('  Passive Vars: %d' % npassive_vars)
                if npassive_vars != 0:
                    self.log.debug('    Has passive vars')
                    # Is variable passive (INT32*NV)
                    data = tecplot_file.read(4 * num_vars)
                    self.n += 4 * num_vars
                    is_var_passive = np.frombuffer(data, dtype=np.int32)

                # Has variable sharing (INT32)
                data = tecplot_file.read(4)
                self.n += 4
                var_sharing, = unpack(b'i', data)
                self.log.debug('  Variable sharing: %d' % var_sharing)
                if var_sharing != 0:
                    # Zero based zone number to share var with (INT32*NV)
                    data = tecplot_file.read(4 * num_vars)
                    self.n += 4 * num_vars
                    zero_based_zone_num = np.frombuffer(data, dtype=np.int32)
                    self.log.debug('    Has variable sharing; 0-based_zone_num=%s' % zero_based_zone_num)

                # Zero based zone number to share connectivity list with (INT32)
                data = tecplot_file.read(4)
                self.n += 4
                zone_number_to_share_connectivity_list, = unpack(b'i', data)
                self.log.debug('  Zero based zone number to share connectivity list with: %d' % zone_number_to_share_connectivity_list)

                # List of min/max pairs for each non-shared and
                # non-passive variables (FLOAT64)
                data = tecplot_file.read(8 * 2 * num_vars)
                self.n += 8 * 2 * num_vars
                var_minMax = np.frombuffer(data, dtype=np.float64).reshape(num_vars, 2)
                self.log.debug('  Min/Max of vars:\n%r' % var_minMax)

                #5. Cell centered variable (DATA SECTION)
                #    To make reading of cell centered binary data efficient, Tecplot stores
                #    IMax*JMax*KMax numbers of cell centered values, where IMax, JMax,
                #    and KMax represent the number of points in the I, J, and K directions.
                #    Therefore extra zero values (ghost values) are written to the data file
                #    for the slowest moving indices. For example, if your data’s IJK
                #    dimensions are 2x3x2, a cell-centered variable will have 1x2x1
                #    (i.e. (I-1)x(J-1)x(K-1)) significant values. However, 2x3x2 values must
                #    be written out because it must include the ghost values. Assume that the
                #    two significant cell-centered values are 1.5 and 12.5. The ghost values
                #    will be output with a zero value.
                #    So if the zone was dimensioned 2x3x2 its cell centered variable would be
                #    represented as follows:
                #    1.5 0.0 12.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                #    If the zone was dimensioned 3x2x2 its cell centered variable would be
                #    represented as follows:
                #    1.5 12.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                #    and if the zone was dimensioned 2x2x3 its cell centered variable would be
                #    represented as follows:
                #    1.5 0.0 0.0 0.0 12.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                #    For large variables the wasted space is less significant that it
                #    is for the small example above.
                nrows = 5
                num_vars = 733 *2  #  733
                nbytes = 4 * nrows * num_vars
                data = tecplot_file.read(nbytes)
                self.n += nbytes
                    #8.102624
                vals = np.frombuffer(data, dtype='float32').reshape(nrows, num_vars).T
                for i, val in enumerate(vals):
                    print(i, val.tolist())

                #print(vals)
                self.show(500, types='ifs')
                if zone_type == 'ORDERED':
                    # if 'zone number to share connectivity list with' == -1 and
                    #    'num of misc. user defined face neighbor connections' != 0
                    #    +-----------+
                    #    |  INT32*N  | Face neighbor connections.
                    #    +-----------+ N = (number of miscellaneous user defined
                    #                       face neighbor connections) * P
                    #     (See note 5 below).
                    if zone_number_to_share_connectivity_list == -1:
                        if num_misc_user_defned_face_neighbor_connections_int != 0:
                            #TODO
                            pass
                    raise NotImplementedError(zone_type)
                elif 'FE' in zone_type:
                    if zone_type not in ['FEPOLYGON', 'FEPOLYHEDRON']:
                        #if 'zone number to share connectivity lists with' == -1
                        #+-----------+
                        #|  INT32*N  | Zone Connectivity Data N=L*JMax
                        #+-----------+ (see note 2 below ).

                        # Zone connectivity (INT32*N);  N=L*JMax
                        #2. This represents JMax sets of adjacency zero based indices where each
                        #    set contains L values and L is:
                        #    2 for LINESEGS
                        #    3 for TRIANGLES
                        #    4 for QUADRILATERALS
                        #    4 for TETRAHEDRONS
                        #    8 for BRICKS
                        if 0:
                            self.log.debug('  Reading zone connectivity data')
                            npoints_L_map = {
                                #'LINE' : 2,
                                'FETRIANGLE' : 3,
                                'FEQUADRILATERAL' : 4,
                                #'FETET' : 4,
                                #'FEBRICK' : 8,
                            }
                            npoints_L = npoints_L_map[zone_type]

                            self.log.info('npoints_L=%r nelements=%s' % (npoints_L, nelements))
                            if nelements:
                                data = tecplot_file.read(4 * npoints_L * nelements)
                                connectivity = np.frombuffer(data, np.int32).reshape(nelements,  npoints_L)


                        #-----------------

                        if zone_number_to_share_connectivity_list == -1 and  are_raw_local_one_to_one_face_neighbors == 1:
                            note_3
                            #if 'zone number to share connectivity lists with' == -1 &&
                            #   'raw local 1-to-1 face neighbors are supplied'
                            #
                            #+-----------+ Raw local 1-to-1 face neighbor array.
                            #|  INT32*N  | N = (NumElements * NumFacesPerElement)
                            #+-----------+ (See note 3 below).
                            #
                            #3. The raw face neighbor array is dimensioned by (number of elements for
                            #   the zone) times (the number of faces per element), where each member
                            #   of the array holds the zero-based element neighbor of that face. A
                            #   boundary face is one that has no neighboring element and is
                            #   represented by a -1. Faces should only be neighbors if they logically
                            #   share nodes and they should be reciprocal.

                        #if zone_number_to_share_connectivity_list == -1 and num_misc_user_defned_face_neighbor_connections_int != 0:
                            #if 'zone number to share connectivity lists with' == -1 &&
                            #   'num of misc. user defined face neighbor connections' != 0
                            #
                            #+-----------+ Face neighbor connections.
                            #|  INT32*N  | N = (number of miscellaneous user defined
                            #+-----------+      face neighbor connections) * P
                            #              (See note 4 below).
                            #
                            # 4. FaceNeighbor Mode # values Data
                            #    ---------------------------------------------------------------------
                            #    LocalOneToOne 3 cz,fz,cz
                            #    LocalOneToMany nz+4 cz,fz,oz,nz,cz1,cz2,...,czn
                            #    GlobalOneToOne 4 cz,fz,ZZ,CZ
                            #    GlobalOneToMany 2*nz+4 cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
                            #
                            #    Where:
                            #    cz = cell in current zone (zero based)
                            #    fz = face of cell in current zone (zero based)
                            #    oz = face obscuration flag (only applies to one-to-many):
                            #    0 = face partially obscured
                            #    1 = face entirely obscured
                            #    nz = number of cell or zone/cell associations
                            #    (only applies to one-to-many)
                            #    ZZ = remote Zone (zero based)
                            #    CZ = cell in remote zone (zero based)
                            #
                            #    cz,fz combinations must be unique and multiple entries are
                            #    not allowed. Additionally, Tecplot assumes that with the
                            #    one-to-one face neighbor modes, a supplied cell face is
                            #    entirely obscured by its neighbor. With one-to-many, the
                            #    obscuration flag must be supplied.
                            #
                            #    Face neighbors that are not supplied are run through
                            #    Tecplot’s auto face neighbor generator (FE only).
                            #note4
                            #
                    elif zone_type in ['FEPOLYGON', 'FEPOLYHEDRON']:
                        pass
                        # if 'zone number to share face map data with' == -1
                        #
                        #+-----------+ Face node offsets into the face nodes array
                        #|  INT32*F  | below. Does not exist for FEPOLYGON zones.
                        #+-----------+ F = NumFaces+1.
                        #
                        #
                        #+-----------+ Face nodes array containing the node numbers
                        #|  INT32*FN | for all nodes in all faces.
                        #+-----------+ FN = total number of face nodes.
                        #
                        #+-----------+ Elements on the left side of all faces.
                        #|  INT32*F  | Boundary faces use a negative value which is
                        #+-----------+ the negated offset into the face boundary
                        #              connection offsets array. A value of '-1'
                        #              indicates there is no left element.
                        #              F = NumFaces.
                        #
                        #+-----------+ Elements on the right side of all faces. See
                        #|  INT32*F  | description of left elements above for more
                        #+-----------+ details. F = NumFaces.
                        #
                        raise RuntimeError(zone_type)
                        #if 'total number of boundary faces' != 0:
                            #raise RuntimeError('total number of boundary faces=%s' % total_number_of_boudary_faces)
                        #if 'total number of boundary faces' != 0
                        #+-----------+ Boundary face connection offsets into the
                        #| INT32*NBF | boundary face connecion elements array and
                        #+-----------+ the boundary face connection zones array.
                        #              The number of elements for a face (F) is
                        #              determined by offset[-o] - offset[-o-1]
                        #              where ‘o’ is the negative value from either
                        #              the left or right elements arrays above.
                        #              Offset[0] = 0. Offset[1] = 0 so that -1 as
                        #              the left or right element always indicates
                        #              no neighboring element. If the number of
                        #              elements is 0, then there is no neighboring
                        #              element.
                        #              NBF = total number of boundary faces + 1.
                        #
                        #
                        #+-----------+ Boundary face connection elements. A value of
                        #| INT32*NBI | '-1' indicates there is no element on part of
                        #+-----------+ the face.
                        #              NBI = total number of boundary connections.
                        #
                        #+-----------+ Boundary face connection zones. A value of
                        #| INT32*NBI | '-1' indicates the current zone.
                        #+-----------+ NBI = total number of boundary connections.
                        #


                    for ivar, variable_name in enumerate(variable_names):
                        dtype = variable_data_format[ivar]
                        self.log.debug('  Reading variable %s; npoints=%s' % (variable_name, npoints))
                        #npoints = zone['NPNTS']
                        if dtype == 'double':
                            sdtype = 'd'
                            nwords = 2
                        elif dtype == 'float':
                            nwords = 1
                            sdtype = 'f'
                        else:
                            raise NotImplementedError(dtype)
                        data = tecplot_file.read(4 * nwords * npoints)
                        self.n += 4 * nwords * npoints
                        vals = np.frombuffer(data, dtype=sdtype)
                        #print(npoints, vals)
                        variables[variable_name] = vals
                    print(variables)
                else:
                    raise NotImplementedError(zone_type)

                #if debug:
                    #debug_file.write('%12s  ' % 'Point')
                    #for v in vars:
                        #debug_file.write('%12s  ' % v)
                    #debug_file.write('\n')
                    #for n in range(zone['NPNTS']):
                        #debug_file.write('%12i  ' % (n+1))
                        #for v in vars:
                            #debug_file.write('%12g  ' % variables[v][n])
                        #debug_file.write('\n')

                #end_of_zone
            end_of_header
            nbytes = 11 * 4
            data = tecplot_file.read(nbytes)
            self.n += nbytes
            #self.show_data(data, types='if', endian='<') # 'if'?
            s = unpack('2f 9i', data)
            self.log.debug(s)
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
                xyzvals = unpack(b'%sf' % ni, data)
                xyz = array(xyzvals, dtype='float32').reshape(3, nnodes).T

                # the variables: [rho, u, v, w, p]
                nvars = 5
                dunno = 0    # what's with this...
                ni = nnodes * nvars + dunno
                nbytes = ni * 4
                data = tecplot_file.read(nbytes)
                self.n += nbytes
                resvals = unpack(b'%sf' % ni, data)
                results = array(resvals, dtype='float32').reshape(nvars, nnodes).T


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

                elements = array(node_ids).reshape(nelements, nnodes_per_element)
                #print(elements)

                #self.show_data(data, types='ifs', endian='<')
                #print(vals)

                #self.show(100, endian='<')
                if zone_type == 5:
                    self.hexa_elements = elements
                elif zone_type == 0:
                    self.quad_elements = elements
                else:
                    raise NotImplementedError(zone_type)
            del self.f

        self.xyz = xyz
        self.results = results
        self.log.debug('done...')

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
        x = self.xyz[:, 0]
        self._slice_plane(x, xslice)

    def slice_y(self, yslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        y = self.xyz[:, 1]
        self._slice_plane(y, yslice)

    def slice_z(self, zslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        z = self.xyz[:, 2]
        self._slice_plane(z, zslice)

    def slice_xyz(self, xslice, yslice, zslice):
        """TODO: doesn't remove unused nodes/renumber elements"""
        x = self.xyz[:, 0]
        y = self.xyz[:, 1]
        z = self.xyz[:, 2]

        inodes = []
        if xslice is not None:
            xslice = float(xslice)
            inodes.append(where(x < xslice)[0])
        if yslice is not None:
            yslice = float(yslice)
            inodes.append(where(y < yslice)[0])
        if zslice is not None:
            zslice = float(zslice)
            inodes.append(where(z < zslice)[0])

        nodes = None
        if len(inodes) == 1:
            nodes = inodes[0]
        elif len(inodes) == 2:
            nodes = intersect1d(inodes[0], inodes[1], assume_unique=True)
        elif len(inodes) == 3:
            nodes = intersect1d(
                intersect1d(inodes[0], inodes[1], assume_unique=True),
                inodes[2], assume_unique=True)
            #inodes = arange(self.nodes.shape[0])
            # nodes = unique(hstack(inodes))
        if nodes is not None:
            self._slice_plane_inodes(nodes)

    def _slice_plane(self, y, slice_value):
        """
        - Only works for CHEXA
        - Doesn't remove unused nodes/renumber elements
        """
        slice_value = float(slice_value)
        inodes = where(y < slice_value)[0]
        self._slice_plane_inodes(inodes)

    def _slice_plane_inodes(self, inodes):
        """
        TODO: doesn't remove unused nodes/renumber elements
        """
        # old_num = inodes
        # new_num = arange(self.xyz.shape[0], dtype='int32')
        #print('old =', old_num)
        #print('new =', new_num)

        # lookup_table = dict( zip( old_num, new_num ) ) # create your translation dict
        # vect_lookup = np.vectorize( lookup_table.get ) # create a function to do the translation

        nhexas = self.hexa_elements.shape[0]
        if nhexas:
            #boolean_hexa = self.hexa_elements.ravel() == inodes
            #boolean_hexa = (self.hexa_elements.ravel() == inodes)#.all(axis=1)
            boolean_hexa = in1d(self.hexa_elements.ravel(), inodes).reshape(nhexas, 8)
            # print(boolean_hexa)
            # assert len(boolean_hexa) == self.hexa_elements.shape[0]
            assert True in boolean_hexa
            irow = where(boolean_hexa)[0]
            isave = unique(irow)
            nsave = len(isave)
            self.hexa_elements = self.hexa_elements[isave, :]
            # print(self.hexa_elements)
            #self.hexa_elements =

            # vect_lookup(self.hexa_elements) # Reassign the elements you want to change
            self.hexa_elements.reshape(nsave, 8)

        # print(boolean_hexa)
        # for hexa in hexas:
            # if
        # self.hexa_elements

    def write_tecplot(self, tecplot_filename, res_types=None, is_points=True,
                      adjust_nids=True):
        """
        Only handles single type writing

        Parameters
        ----------
        tecplot_filename : str
            the path to the output file
        res_types : str; List[str, str, ...]; default=None -> all
            the results that will be written (must be consistent with
            self.variables)
        is_points : bool; default=True
            write in POINT format vs. BLOCK format
        adjust_nids : bool; default=True
            element_ids are 0-based in binary and must be switched to
            1-based in ASCII
        """
        self.log.info('writing tecplot %s' % tecplot_filename)
        with open(tecplot_filename, 'w') as tecplot_file:
            is_results = bool(len(self.results))
            msg = 'TITLE     = "tecplot geometry and solution file"\n'
            msg += 'VARIABLES = "x"\n'
            msg += '"y"\n'
            msg += '"z"\n'
            if res_types is None:
                res_types = self.variables
            elif isinstance(res_types, string_types):
                res_types = [res_types]
            result_indices_to_write = []
            if is_results:
                #msg += '"rho"\n'
                #msg += '"u"\n'
                #msg += '"v"\n'
                #msg += '"w"\n'
                #msg += '"p"\n'
                # msg += 'ZONE T="%s"\n' % r'\"processor 1\"'
                # print('res_types =', res_types)
                # print('vars =', self.variables)
                for ivar, var in enumerate(res_types):
                    if var not in self.variables:
                        raise RuntimeError('var=%r not in variables=%s' % (var, self.variables))
                    result_indices_to_write.append(self.variables.index(var))
                ivars = unique(result_indices_to_write)
                ivars.sort()
                for ivar in ivars:
                    var = self.variables[ivar]
                    msg += '"%s"\n' % var
                # print('ivars =', ivars)
            else:
                assert len(res_types) == 0, len(res_types)
                ivars = []
            msg += 'ZONE '

            etype_elements = [
                ('CHEXA', self.hexa_elements),
                ('CTETRA', self.tet_elements),
                ('CTRIA3', self.tri_elements),
                ('CQUAD4', self.quad_elements),
            ]
            is_points = True
            is_tets = False
            is_hexas = False
            is_tris = False
            is_quads = False

            nnodes = self.nnodes
            nelements = self.nelements
            for etype, elements in etype_elements:
                if not len(elements):
                    continue
                if etype == 'CHEXA' and len(elements):
                    #print(etype)
                    # is_points = False
                    is_hexas = True
                    #nnodes_per_element = 8
                    zone_type = 'FEBrick'
                elif etype == 'CTETRA' and len(elements):
                    #print(etype)
                    # is_points = False
                    is_tets = True
                    #nnodes_per_element = 4
                    zone_type = 'FETETRAHEDRON'
                elif etype == 'CTRIA3' and len(elements):
                    #print(etype)
                    # is_points = True
                    is_tris = True
                    #nnodes_per_element = 3
                    zone_type = 'FETRIANGLE'
                elif etype == 'CQUAD4' and len(elements):
                    #print(etype)
                    # is_points = True
                    is_quads = True
                    #nnodes_per_element = 4
                    zone_type = 'FEQUADRILATERAL'
                else:
                    self.log.info('etype=%r' % etype)
                    self.log.info(elements)
                    continue
                break

            self.log.info('is_points = %s' % is_points)
            if is_points:
                msg += ' n=%i, e=%i, ZONETYPE=%s, DATAPACKING=POINT\n' % (
                    nnodes, nelements, zone_type)
            else:
                msg += ' n=%i, e=%i, ZONETYPE=%s, DATAPACKING=BLOCK\n' % (
                    nnodes, nelements, zone_type)
            tecplot_file.write(msg)

            # xyz
            assert self.nnodes > 0, 'nnodes=%s' % self.nnodes
            nresults = len(ivars)
            if is_points:
                if nresults:
                    res = self.results[:, ivars]
                    try:
                        data = hstack([self.xyz, res])
                    except ValueError:
                        msg = 'Cant hstack...\n'
                        msg += 'xyz.shape=%s\n' % str(self.xyz.shape)
                        msg += 'results.shape=%s\n' % str(self.results.shape)
                        raise ValueError(msg)
                    fmt = ' %15.9E' * (3 + nresults)
                else:
                    data = self.xyz
                    fmt = ' %15.9E %15.9E %15.9E'

                # works in numpy 1.15.1
                savetxt(tecplot_file, data, fmt=fmt)
            else:
                #nvalues_per_line = 5
                for ivar in range(3):
                    #tecplot_file.write('# ivar=%i\n' % ivar)
                    vals = self.xyz[:, ivar].ravel()
                    msg = ''
                    for ival, val in enumerate(vals):
                        msg += ' %15.9E' % val
                        if (ival + 1) % 3 == 0:
                            tecplot_file.write(msg)
                            msg = '\n'
                    tecplot_file.write(msg.rstrip() + '\n')

                if nresults:
                    # print('nnodes_per_element =', nnodes_per_element)
                    # for ivar in range(nnodes_per_element):
                    for ivar in ivars:
                        #tecplot_file.write('# ivar=%i\n' % ivar)
                        vals = self.results[:, ivar].ravel()
                        msg = ''
                        for ival, val in enumerate(vals):
                            msg += ' %15.9E' % val
                            if (ival + 1) % 5 == 0:
                                tecplot_file.write(msg)
                                msg = '\n'
                        tecplot_file.write(msg.rstrip() + '\n')

            self.log.info('is_hexas=%s is_tets=%s is_quads=%s is_tris=%s' %
                          (is_hexas, is_tets, is_quads, is_tris))
            if is_hexas:
                efmt = ' %i %i %i %i %i %i %i %i\n'
                elements = self.hexa_elements
            elif is_tets:
                efmt = ' %i %i %i %i\n'
                elements = self.tet_elements
            elif is_quads:
                efmt = ' %i %i %i %i\n'
                elements = self.quad_elements
            elif is_tris:
                efmt = ' %i %i %i\n'
                elements = self.tri_elements
            else:
                raise RuntimeError()

            # we do this before the nid adjustment
            node_min = elements.min()
            node_max = elements.max()
            self.log.info('inode: min=%s max=%s' % (node_min, node_max))
            assert node_min >= 0, node_min

            if node_max > nnodes:
                msg = 'elements.min()=node_min=%s elements.max()=node_max=%s nnodes=%s' % (
                    node_min, node_max, nnodes)
                raise RuntimeError(msg)
            # assert elements.min() == 1, elements.min()
            # assert elements.max() == nnodes, elements.max()

            if adjust_nids:
                elements += 1

            for element in elements:
                tecplot_file.write(efmt % tuple(element))

    def skin_elements(self):
        tris = []
        quads = []
        if len(self.tet_elements):
            faces1 = self.tet_elements[:, :3]
            faces2 = self.tet_elements[:, 1:4]
            faces3 = self.tet_elements[:, [2, 3, 0]]
            tris.append(faces1)
            tris.append(faces2)
            tris.append(faces3)

        if len(self.hexa_elements):
            faces1 = self.hexa_elements[:, :4]
            faces2 = self.hexa_elements[:, 4:7]
            assert faces1.shape[1] == 4, faces1.shape
            assert faces2.shape[1] == 4, faces2.shape
            #faces3 = self.hexa_elements[:, [2, 3, 0]]
            # TODO: others CHEXA faces...
            quads.append(faces1)
            quads.append(faces2)

        # if tris:
            # tris = vstack(tris)
            # tris.sort(axis=0)
            # tris = unique_rows(tris)
        # if quads:
            # quads = vstack(quads)
            # quads.sort(axis=0)
            # quads = unique_rows(tris)
        return tris, quads

    def get_free_faces(self):
        """get the free faces for hexa elements"""
        self.log.info('start get_free_faces')
        sort_face_to_element_map = defaultdict(list)
        sort_face_to_face = {}
        for ie, element in enumerate(self.hexa_elements):
            btm = [element[0], element[1], element[2], element[3]]
            top = [element[4], element[5], element[6], element[7]]
            left = [element[0], element[3], element[7], element[4]]
            right = [element[1], element[2], element[6], element[5]]
            front = [element[0], element[1], element[5], element[4]]
            back = [element[3], element[2], element[6], element[7]]
            for face in [btm, top, left, right, front, back]:
                if len(unique(face)) >= 3:
                    sort_face = tuple(sorted(face))
                    sort_face_to_element_map[sort_face].append(ie)
                    sort_face_to_face[sort_face] = face

        free_faces = []
        for sort_face, eids in sort_face_to_element_map.items():
            if len(eids) == 1:
                free_faces.append(sort_face_to_face[sort_face])
        self.log.info('finished get_free_faces')
        return free_faces

    def extract_y_slice(self, y0, tol=0.01, slice_filename=None):
        """
        doesn't work...
        """
        self.log.info('slicing...')
        y = self.xyz[:, 1]
        nodes = self.xyz
        assert tol > 0.0, tol
        elements = self.hexa_elements
        results = self.results

        iy = where((y0 - tol <= y) & (y <= y0 + tol))[0]

        self.log.debug(y[iy])
        self.log.debug(nodes[iy, 1].min(), nodes[iy, 1].max())
        #iy = where(y <= y0 + tol)[0]
        assert len(iy) > 0, iy
        #inode = iy + 1


        # find all elements that have iy within tolerance
        #slots = where(elements == iy)
        #slots = where(element for element in elements
                      #if any(iy in element))
        #slots = where(iy == elements.ravel())[0]
        ielements = unique([ie for ie, unused_elem in enumerate(elements)
                            for i in range(8)
                            if i in iy])
        #print(slots)
        #ri, ci = slots
        #ri = unique(hstack([where(element == iy)[0] for element in elements]))
        #ri = [ie for ie, element in enumerate(elements)
              #if [n for n in element
                  #if n in iy]]
        #ri = [where(element == iy)[0] for element in elements if where(element == iy)[0]]
        #print(ri)
        #ielements = unique(ri)
        self.log.debug(ielements)
        assert len(ielements) > 0, ielements

        # find nodes
        elements2 = elements[ielements, :]
        inodes = unique(elements2)
        assert len(inodes) > 0, inodes

        # renumber the nodes
        nidmap = {}
        for inode, nid in enumerate(inodes):
            nidmap[nid] = inode
        elements3 = array(
            [[nidmap[nid] for nid in element]
             for element in elements2],
            dtype='int32')

        self.log.debug(inodes)
        nodes2 = nodes[inodes, :]
        results2 = results[inodes, :]
        model = Tecplot()
        model.xyz = nodes2
        model.results = results2
        model.hexa_elements = elements3

        if slice_filename:
            model.write_tecplot(slice_filename)
        return model



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


def _header_lines_to_header_dict(header_lines):
    """parses the parsed header lines"""
    headers_dict = {}
    if len(header_lines) == 0:
        #raise RuntimeError(header_lines)
        return None
    header = ','.join(header_lines)

    # this is so overly complicataed and probably not even enough...
    # what about the following 'quote' style?
    headers = header.replace('""', '","').split(',')

    nheaders = len(headers) - 1
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
            raise NotImplementedError(header)

        if parse:
            #print('  parsing')
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
            allowed_keys = ['TITLE', 'VARIABLES', 'T', 'ZONETYPE', 'DATAPACKING',
                            'N', 'E', 'F', 'DT', 'SOLUTIONTIME', 'STRANDID',
                            'I', 'J']
            assert key in allowed_keys, 'key=%r; allowed=[%s]' % (key, ', '.join(allowed_keys))
            parse = False
    #print(headers_dict.keys())
    assert len(headers_dict) > 0, headers_dict
    return headers_dict


if __name__ == '__main__':
    main()
