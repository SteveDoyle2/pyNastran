# encoding: utf-8
import os
import sys
from struct import Struct, unpack, pack
from collections import namedtuple
from typing import BinaryIO, Optional, Union, Any

import numpy as np
from pyNastran.utils import print_bad_path
from pyNastran.converters.tecplot.zone import Zone, CaseInsensitiveDict
from pyNastran.utils import object_attributes, object_methods, object_stats

from cpylog import get_logger2

ZoneTuple = namedtuple('Zone', ['zone_name',
                                'strand_id', 'solution_time',
                                'zone_type', 'zone_type_str',
                                'data_packing', 'data_packing_str',
                                'celldim',
                                'raw_local',
                                'n_misc_neighbor_connections',
                                'nelement'])

class Base:
    def object_attributes(obj: Any, mode: str='public',
                          keys_to_skip: Optional[list[str]]=None,
                          filter_properties: bool=False) -> list[str]:
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings
        filter_properties: bool: default=False
            filters the @property objects

        Returns
        -------
        attribute_names : list[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong

        """
        return object_attributes(obj,
                                 mode=mode,
                                 keys_to_skip=keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(obj: Any, mode: str='public',
                       keys_to_skip: Optional[list[str]]=None) -> list[str]:
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : list[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong

        """
        return object_methods(obj,
                              mode=mode,
                              keys_to_skip=keys_to_skip)

    def object_stats(obj: Any, mode: str='public',
                 keys_to_skip: Optional[list[str]]=None,
                 filter_properties: bool=False) -> str:
        """Prints out an easy to read summary of the object"""
        return object_stats(obj,
                            mode=mode,
                            keys_to_skip=keys_to_skip,
                            filter_properties=filter_properties)


class TecplotBinary(Base):
    def __init__(self, log=None, debug: bool=False):
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

        self.zones: list[Zone] = []
        # mesh = True : this is a structured/unstructured grid

        # mesh = False : this is a plot file
        self.use_cols = None

        # TODO: what is this for?
        self.dtype = None

        self._uendian = ''
        self.n = 0

    def show(self, n: int, types: str='ifs', endian=None):  # pragma: no cover
        assert self.n == self.f.tell()
        nints = n // 4
        data = self.f.read(4 * nints)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        self.f.seek(self.n)
        return strings, ints, floats

    def show_data(self, data: bytes, types: str='ifs', endian=None):  # pragma: no cover
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

    def _write_data(self, f, data: bytes, types: str='ifs', endian=None):  # pragma: no cover
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
            chrs = [chr(val) if 40 <= val <= 120 else val
                    for val in ints]
            f.write("  ints    = %s\n" % str(ints))
            f.write("  chrs    = %s\n" % str(chrs))
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

    def show_ndata(self, n: int, types: str='ifs'):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types)

    def _write_ndata(self, f, n: int, types: str='ifs'):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.
        """
        nold = self.n
        data = self.f.read(n)
        self.n = nold
        self.f.seek(self.n)
        return self._write_data(f, data, types=types)

    def write_tecplot_binary(self, tecplot_filename: str,
                             version: str='102') -> None:
        assert version == '102', version

        nvars = len(self.variables)
        #sprint('vars =', self.variables)
        assert nvars > 1, self.variables
        with open(tecplot_filename, 'wb') as tecplot_file:
            _write_binary_header(self, tecplot_file, version)
            _write_binary_zone_headers(
                self, tecplot_file, nvars, version)
            _write_binary_results(
                self, tecplot_file, nvars, version)

    def read_tecplot_binary(self, tecplot_filename: str,
                            zones_to_exclude: Optional[list[int]]=None) -> None:
        """
        Supports multiblock, but FEQUADs in BLOCK format only
        This is actually a semi-competent reader, so that's the good news

        http://www.hgs.k12.va.us/tecplot/documentation/tp_data_format_guide.pdf
        """
        set_zones_to_exclude = zones_to_exclude_to_set(zones_to_exclude)
        log = self.log
        assert os.path.exists(tecplot_filename), print_bad_path(tecplot_filename)

        self.n = 0
        with open(tecplot_filename, 'rb') as self.f:
            file_obj = self.f
            out = _read_binary_header(self, file_obj)
            version, header_dict, title, file_type, variables, zone_tuples = out
            del out

            self.is_mesh = True
            self.header_dict = header_dict
            self.variables = variables
            self.title = title

            nvars = len(variables)

            #self.show(1000, types='ifdq')
            for izone, zone_tuple in enumerate(zone_tuples):
                zone_name = zone_tuple.zone_name
                log.info(f'-----izone={izone} {zone_name!r}-----')
                tris, quads, tets, hexas, zone_data = _read_binary_results(
                    self, file_obj, zone_tuple, nvars, version)
                if quads is not None:
                    log.debug(f'quads.min/max = {quads.min()} {quads.max()}')
                    assert quads.min() >= 0, quads.min()

                if izone not in set_zones_to_exclude:
                    zone = Zone.set_zone_from_360(
                        log, header_dict, variables,
                        zone_name,
                        strand_id=zone_tuple.strand_id,
                        data_packing=zone_tuple.data_packing_str,
                        zone_type=zone_tuple.zone_type_str,
                        tris=tris, quads=quads,
                        tets=tets, hexas=hexas,
                        zone_data=zone_data,
                        )
                    self.zones.append(zone)
                #assert quads.max() <= 236_064, (zone_name, quads.max())

            # final check
            data = file_obj.read(1)
            if len(data):
                raise RuntimeError("there is data at the end of the file "
                                   "that wasn't read")
        #print(str(self))
        self.header_dict = zone.headers_dict
        zone.headers_dict['DATAPACKING']
        self.variables = zone.variables
        del self.f


def _write_binary_zone_headers(model: TecplotBinary,
                               tecplot_file: BinaryIO,
                               nvars: int,
                               version: bytes) -> None:
    #end 357.0, 299.0,
    start_of_zone_flag = pack(b'<f', 299.0)
    for zone in model.zones:
        tecplot_file.write(start_of_zone_flag)
        zone_type = zone.zone_type_int

        if version == '102':
            #print(f'writing name={zone.name!r}')
            data = _write_string(zone.name)
            tecplot_file.write(data)

            data = pack('<2i', zone.strand_id, zone_type)
            tecplot_file.write(data)

            data = pack('<3i', 0, 0, 0)
            tecplot_file.write(data)

            data = pack('<2i', zone.nnodes, zone.nelements)
            tecplot_file.write(data)

            data = pack('<4i', 0, 0, 0, 0)
            tecplot_file.write(data)

            #specify_var = 0
            #data_packing = 0  # BLOCK
            raw_local = 0
            solution_time = np.nan
        elif version == '112':
            raise RuntimeError(version)
            assert -2 <= strand_id < 32700, strand_id
            print(f'  parent_zone={parent_zone} strand_id={strand_id} solution_time={solution_time:g}')
            print(f'  zone_type={zone_type} data_packing={data_packing} specify_var={specify_var}')
            data = unpack('<iid 4i', (parent_zone, strand_id, solution_time,
                 unused_a, zone_type, data_packing, specify_var))
            tecplot_file.write(data)

            ordered_zone = is_ordered_zone(zone_type)
            fe_poly = is_fe_poly_zone(zone_type)
            #specify_var = 1
            if specify_var:
                ndatai = 4 * nvars; data = file_obj.read(ndatai); n += ndatai
                specify_varsi = np.frombuffer(data, dtype='int32')

            print(f'  raw_local={raw_local} n_misc_neighbor_connections={n_misc_neighbor_connections}')
            assert raw_local in {0, 1}, raw_local
            assert n_misc_neighbor_connections >= 0, n_misc_neighbor_connections
            data = unpack('<2i', raw_local, n_misc_neighbor_connections)
            tecplot_file.write(data)

            # Are raw local 1-to-1 face neighbors supplied?
            # (0=FALSE 1=TRUE).
            if is_fe_zone(zone_type):
                pass
                #show_ndata(file_obj, 100, types='ifsqd', endian='<')
            if ordered_zone:
                raise RuntimeError('ordered zone; p.153')
            if fe_poly:
                raise RuntimeError(' FEPOLYGON or FEPOLYHEDRON; p.153')

            print(f'  nelement={nelement} celldim={celldim}')
            data = unpack('<4i', nelement, *celldim)
            tecplot_file.write(data)

            data = file_obj.read(4); n += 4
        else:
            raise NotImplementedError(version)
    return

def _write_binary_results(model: TecplotBinary,
                          tecplot_file: BinaryIO,
                          nvars: int, version: str) -> None:
    #elif flag == 357.0:
    # EOHMARKER (end of header marker)

    assert isinstance(version, str), version

    #end 357.0, 299.0,
    start_of_results_flag = pack(b'<f', 357.0)
    start_of_zone_flag = pack(b'<f', 299.0)

    #print('writing 357.0')
    tecplot_file.write(start_of_results_flag)

    for zone in model.zones:
        #print('writing 299.0')
        tecplot_file.write(start_of_zone_flag)
        zone_type = zone.zone_type_int

        # INT32*N
        # Variable data format, N=Total number of vars
        # 1=Float, 2=Double, 3=LongInt,
        # 4=ShortInt, 5=Byte, 6=Bit
        results = []
        dtypes = []
        _prep_results_dtypes(results, dtypes, zone.zone_data)

        data = pack(f'<{nvars}i', *dtypes)
        tecplot_file.write(data)

        has_passive_vars = 0
        data = pack(b'<i', has_passive_vars)
        tecplot_file.write(data)

        connectivity_sharing = -1
        data = pack(b'<i', connectivity_sharing)
        tecplot_file.write(data)

        if version != '102': # definitely right
            # Tecplot 10
            min_max_data = _get_min_max(results)
            tecplot_file.write(min_max_data)

        ## Zone Data. Each variable is in data format as
        ## specified above
        for result in results:
            datai = result.tobytes()
            tecplot_file.write(datai)

        if zone_type == 2: # FETRIA
            elements = zone.tri_elements
        elif zone_type == 3: # FEQUAD
            elements = zone.quad_elements
        elif zone_type == 4: # FETETRAHEDRON
            elements = zone.tet_elements
        elif zone_type == 5:  # FEBRICK
            elements = zone.hexa_elements
        else:
            raise RuntimeError(zone_type)
        tecplot_file.write(elements.ravel().tobytes())
    return


def _prep_results_dtypes(results_list: list[np.ndarray],
                         dtypes: list[int],
                         data: np.ndarray,) -> None:
    nvars = data.shape[1]
    assert nvars >= 1, nvars
    for i in range(nvars):
        results_list.append(data[:, i])
    dtypes += _get_dtypes_from_array(data)
    assert len(dtypes) > 0

def _read_binary_results(
        self: TecplotBinary,
        file_obj: BinaryIO,
        zone: ZoneTuple,
        nvars: int,
        version: bytes) -> tuple[np.ndarray,  # xyz
                                 Optional[np.ndarray], # tris
                                 Optional[np.ndarray], # quads
                                 Optional[np.ndarray], # tets
                                 Optional[np.ndarray], # hexas
                                 ]:
    """
    II. DATA SECTION (don’t forget to separate the header from the data
    with an EOHMARKER). The data section contains all of the data
    associated with the zone definitions in the header.
    """
    log = self.log
    tris = None
    quads = None
    tets = None
    hexas = None
    zone_type = zone.zone_type
    data_packing = zone.data_packing
    celldim = zone.celldim
    raw_local = zone.raw_local
    n_misc_neighbor_connections = zone.n_misc_neighbor_connections
    nelement = zone.nelement
    del celldim

    #if 'has passive variables' != 0
    #
    #INT32*NV
    #Is variable passive: 0 = no, 1 = yes
    #(Omit entirely if 'Has passive variables' is 0).
    #
    #INT32
    #Has variable sharing 0 = no, 1 = yes.
    #if 'has variable sharing' != 0
    #
    #INT32*NV
    #Zero based zone number to share variable with
    #(relative to this datafile). (-1 = no sharing).
    #(Omit entirely if 'Has variable sharing' is 0).
    #
    #INT32
    #Zero based zone number to share connectivity
    #list with (-1 = no sharing). FEPOLYGON and
    #FEPOLYHEDRON zones use this zone number to
    #share face map data

    #0=ORDERED,
    #1=FELINESEG,
    #2=FETRIANGLE,
    #3=FEQUADRILATERAL,
    #4=FETETRAHEDRON,
    #5=FEBRICK,
    #6=FEPOLYGON,
    #7=FEPOLYHEDRON
    #
    # NOTE 2. This represents JMax sets of adjacency zero based indices where each
    #    set contains L values and L is
    #    2 for LINESEGS
    #    3 for TRIANGLES
    #    4 for QUADRILATERALS
    #    4 for TETRAHEDRONS
    #    8 for BRICKS
    #
    #L = 0
    #if zone_type in {3, 4}:
        # FEQUADRILATERAL, FETETRAHEDRON
        #L = 4
    #elif zone_type == 5: # BRICKS
        #L = 8
    #else:
        #raise RuntimeError(zone_type)

    # :--------------------i. For both ordered and fe zones:--------------------
    # FLOAT32
    # Zone marker Value = 299.0
    ndata = 4; data = file_obj.read(ndata); self.n += ndata
    marker, = unpack('<f', data)
    assert marker == 299.0, marker
    log.debug(f'marker = {marker}')

    #  INT32*N
    # Variable data format, N=Total number of vars
    # 1=Float, 2=Double, 3=LongInt,
    # 4=ShortInt, 5=Byte, 6=Bit
    ndata = 4 * nvars; data = file_obj.read(ndata); self.n += ndata
    fmt = '<%ii' % nvars
    data_fmt = unpack(fmt, data)
    #assert data_fmt in {0, 1, 2, 3, 4, 5, 6}
    log.debug(f'data_fmt={data_fmt}')
    #show_ndata(file_obj, 100, types='ifsdq')

    #if version == b'102':
        ## Tecplot 10
        #has_var_sharing = 0
        #connectivity_sharing = -1
    #else:
    #INT32
    # Has passive variables: 0 = no, 1 = yes.
    ndata = 4; data = file_obj.read(ndata); self.n += ndata
    has_passive_vars, = unpack('<i', data)
    assert has_passive_vars in {0, 1}
    #print(f'has_passive_vars={has_passive_vars}')
    if has_passive_vars:
        # INT32*NV
        # Is variable passive: 0 = no, 1 = yes
        # (Omit entirely if 'Has passive variables' is 0).
        ndata = 4 * nvars; data = file_obj.read(ndata); self.n += ndata
        fmt = '<%ii' % nvars
        is_passive = unpack(fmt, data)
        #print(f'is_passive = {is_passive}')

    if version != b'102':  #  suspect...could be something else (like no passive vars) for 10.2
        # trying to parse (0, -1, 0) in POINT_FEBRICK_3D_02, which probably has connectivity_sharing=-1
        # Has variable sharing 0 = no, 1 = yes.
        ndata = 4; data = file_obj.read(ndata); self.n += ndata
        has_var_sharing, = unpack('<i', data)
        assert has_var_sharing in {0, 1}
        #print(f'has_var_sharing={has_var_sharing}')
        if has_var_sharing:
            # p.157 in Tecplot 2015
            ndatai = 4 * nvars; data = file_obj.read(ndatai); self.n += ndatai
            var_sharingi = np.frombuffer(data, dtype='int32')
            assert var_sharingi.min() == -1 and var_sharingi.max() == -1, var_sharingi

            #raise NotImplementedError((has_var_sharing, var_sharingi))

    #INT32
    # Zero based zone number to share connectivity
    # list with (-1 = no sharing). FEPOLYGON and
    # FEPOLYHEDRON zones use this zone number to
    # share face map data.
    ndata = 4; data = file_obj.read(ndata); self.n += ndata
    connectivity_sharing, = unpack('<i', data)
    assert connectivity_sharing in {-1}
    #print(f'connectivity_sharing={connectivity_sharing}')

    fe_poly = is_fe_poly_zone(zone_type)
    if is_ordered_zone(zone_type):
        raise RuntimeError('is_ordered; p.158')

    if is_fe_zone(zone_type):
        if connectivity_sharing != -1:
            raise RuntimeError('is_fe_zone; p.158')

        assert zone_type in {2, 3, 4, 5}, zone_type # fequad
        assert data_packing in {0, 1}, data_packing # block
        assert raw_local == 0, raw_local
        assert n_misc_neighbor_connections > 0, n_misc_neighbor_connections

        if n_misc_neighbor_connections:
            # nodal output
            nnodes = n_misc_neighbor_connections
            #print(f'nelement={nelement} nnodes={nnodes}')
            # NumElements * NumFacesPerElement

        if version != b'102': # definitely right
            # Tecplot 10
            min_max_data, ndata = _read_binary_min_max(file_obj, nvars)
            self.n += ndata

        ## Zone Data. Each variable is in data format as
        ## specified above.
        nvalues = nvars * nnodes
        zone_data0, ndata = _load_binary_data(file_obj, nvalues, data_fmt)
        zone_data = zone_data0.reshape((nvars, nnodes)).T
        self.n += ndata

        assert file_obj.tell() == self.n
        #find_ints(file_obj, self.n)

        if zone_type == 2: # FETRIANGLE
            nints = 3 * nelement
            tris, ndatai = _load_binary_data(file_obj, nints, data_fmt='int32')
            tris = tris.reshape(nelement, 3)
        elif zone_type == 3: # FEQUADRILATERAL
            nints = 4 * nelement
            quads, ndatai = _load_binary_data(file_obj, nints, data_fmt='int32')
            quads = quads.reshape(nelement, 4)
        elif zone_type == 4:  # FETETRAHEDRON
            nints = 4 * nelement
            tets, ndatai = _load_binary_data(file_obj, nints, data_fmt='int32')
            tets = tets.reshape(nelement, 4)
        elif zone_type == 5:  # FEBRICK
            nints = 8 * nelement
            hexas, ndatai = _load_binary_data(file_obj, nints, data_fmt='int32')
            hexas = hexas.reshape(nelement, 8)
        else:
            raise RuntimeError(zone_type)
        self.n += ndatai

    #self.show(10000, types='ifd')
    #self.show(10000, types='if')

    del fe_poly
    return tris, quads, tets, hexas, zone_data
    # :--------------------ii. specific to ordered zones:--------------------
    if is_ordered_zone(zone_type):
        raise RuntimeError('is_ordered; p.158')

    # :--------------------iii. specific to fe zones:--------------------
    if is_fe_zone(zone_type):
        if not is_fe_poly_zone(zone_type):
            N = L*JMax # Note 2
            #raise RuntimeError('is_fe_poly_zone; p.158')

        assert zone_type == 3, zone_type # fequad
        assert data_packing == 0, data_packing # block
        assert raw_local == 0, raw_local
        assert n_misc_neighbor_connections > 0, n_misc_neighbor_connections
        x = 1
    raise RuntimeError('you should have already returned....this should never happen')

def find_ints(file_obj: BinaryIO, n0: int) -> int: # pragma: no cover
    """
    The Tecplot 360 manual is unclear, so to figure out the data
    structure we can just search for known data. It's best to search
    for elements (they're ints).  If I know that there are 100,000
    values prior to the element data, I can figure out that that data
    has to be a lot of results first
    """
    # first_quad is the element we're looking for
    # values are 0 - based, so first_quad -> first_quad2
    first_quad = [38913, 7087, 7094, 7094]
    first_quad2 = [val - 1 for val in first_quad]
    first_val = first_quad2[0]

    #nvalues = 1_000_000
    nvalues = 921_000
    ndata = nvalues * 4; datai = file_obj.read(ndata)
    ints = np.frombuffer(datai, dtype='int32')
    i0s = np.where(ints == first_val)[0]
    for i0 in i0s:
        print('  ', i0, ints[i0:i0+4])
    file_obj.seek(n0)

    #i0 = 920705
    #ndata2 = 3682712
    ndata2 = i0s[0] * 4
    print('i0 =', i0, ndata2)
    return ndata2

def _load_binary_data(file_obj: BinaryIO,
                      nvalues: int,
                      data_fmt: Union[list[int], int, str]) -> tuple[np.ndarray, int]:
    """
    Variable data format
    N=Total number of vars
    1=Float, 2=Double, 3=LongInt,
    4=ShortInt, 5=Byte, 6=Bit

    Parameters
    ----------
    nvalues : int
        number of values to read
    data_fmt : list[int], int, str
        str: int32, int64, float32, float64
        int: 1, 2, 3
        data_fmt: [1, 1, 1, ...], [2, 2, 2, ...], [3, 3, 3, ...]

    TODO: doesn't support mixed data_fmt types
    TODO: there is no official int32 option using the integer format?
    TODO: doesn't support 4, 5, 6 -> (shortInt, Byte, Bit)
    """
    data_fmts = _data_fmt_to_list(data_fmt)
    nvars = len(data_fmts)
    nvalues_per_variable = nvalues // nvars

    ndatas, dtypes = _data_fmt_ndatas_dtype(data_fmts)

    # second check is to force a crash to test
    if max(data_fmts) == min(data_fmts):# and len(data_fmts) == 1:
        data_fmt = data_fmt[0]
        dtype = dtypes[0]
        ndatai = ndatas[0]
        assert nvalues % nvars == 0, (nvalues, nvars, nvalues % nvars)

        ndatas_to_read = (nvalues_per_variable * ndatas).sum()
        ndata = nvalues * ndatai
        assert ndatas_to_read == ndata

        data = file_obj.read(ndata)
        data_array = np.frombuffer(data, dtype=dtype)
    else:
        raise RuntimeError('Mixed types not supported\n'
                           f'data_fmt={data_fmt} -> dtypes={dtypes}')
    return data_array, ndata

def _data_fmt_ndatas_dtype(data_fmts: list[int]) -> tuple[np.ndarray, list[str]]:
    """
    TODO: there is no official int32 option using the integer format?
    TODO: doesn't support 4, 5, 6 -> (shortInt, Byte, Bit)
    """
    ndata_fmts = len(data_fmts)
    ndatas = np.zeros(ndata_fmts, dtype='int32')
    dtypes = []
    for i, data_fmt in enumerate(data_fmts):
        if data_fmt == -1:  # faked
            ndata = 4
            dtype = 'int32'
        elif data_fmt == 1:
            ndata = 4
            dtype = 'float32'
        elif data_fmt == 2:
            dtype = 'float64'
            ndata = 8
        elif data_fmt == 3:
            dtype = 'int64'
            ndata = 8
        else:
            raise NotImplementedError(data_fmt)
        ndatas[i] = ndata
        dtypes.append(dtype)
    return ndatas, dtypes

def _data_fmt_to_list(data_fmt: Union[str, int, list[int]]) -> list[int]:
    """
    puts the data in tecplot form

    TODO: there is no official int32 option using the integer format?
    TODO: doesn't support 4, 5, 6 -> (shortInt, Byte, Bit)
    """
    if isinstance(data_fmt, str):
        if data_fmt == 'int32':
            data_fmt = -1  # faked
        elif data_fmt == 'float32':
            data_fmt = 1
        elif data_fmt == 'float64':
            data_fmt = 2
        elif data_fmt == 'int64':
            data_fmt = 3
        else:
            raise RuntimeError(data_fmt)

    if isinstance(data_fmt, int):
        data_fmts = [data_fmt]
    else:
        data_fmts = data_fmt
    return data_fmts

def _get_dtypes_from_array(data: np.ndarray) -> list[int]:
    """
     INT32*N
    Variable data format, N=Total number of vars
    1=Float, 2=Double, 3=LongInt,
    4=ShortInt, 5=Byte, 6=Bit
    """
    assert len(data.shape) == 2, data.shape
    dtype = data.dtype.name
    nvars = data.shape[1]
    if dtype == 'float32':
        dtypes = [1] * nvars
    elif dtype == 'float64':
        dtypes = [1] * nvars
    else:
        raise RuntimeError(dtype)
    return dtypes

def _read_binary_zone_headers(
        self,
        file_obj: BinaryIO,
        n: int, nvars: int,
        version: bytes) -> [int, list[ZoneTuple]]:
    """TODO: has an artificial cap of 1000 zones"""
    zones = []
    izone = -1
    log = self.log
    while izone < 1000:
        izone += 1
        #self.n = n; self.show_ndata(120, types='ifs')
        data = file_obj.read(4); n += 4
        assert n == file_obj.tell()
        zone_marker = unpack('<f', data)[0]
        if zone_marker != 299.0:
            log.debug('end of zones! %s' % zone_marker)
            #self.n = n; self.show(1000, types='ifs')
            n -= 4
            file_obj.seek(n)
            break
        else:
            log.debug('new zone %s' % zone_marker)
        assert zone_marker == 299.0, zone_marker
        assert n == file_obj.tell()
        #self.n = n; show(file_obj, 100, types='ifsd')
        # Zone name.
        # N = (length of zone name) + 1.
        zone_name, n = _read_string(file_obj, n)
        log.debug(f'zone_name = {zone_name!r}')
        assert n == file_obj.tell()
        #self.n = n; show(file_obj, 80, types='ifsd')

        # ParentZone: Zero-based zone number within this
        # datafile to which this zone is
        # a child.
        #
        # StrandID: -2 = pending strand ID for assignment by Tecplot
        #  -1 = static strand ID
        #  0 <= N < 32700 valid strand ID
        #
        # Solution time.
        #
        # Not used. Set to -1

        # ZoneType
        # 0=ORDERED,
        # 1=FELINESEG,
        # 2=FETRIANGLE,
        # 3=FEQUADRILATERAL,****
        # 4=FETETRAHEDRON,
        # 5=FEBRICK,***
        # 6=FEPOLYGON,
        # 7=FEPOLYHEDRON
        #
        # Data packing.
        # 0 = Block
        # 1 = Point
        #
        # Specify Var Location.
        # 0 = Don’t specify, all data is located at the nodes.
        # 1 = Specify

        if version == b'102':
            # per 2008 data format guide
            # BLOCK
            # ints   = (-1, 5, 0, 0, 0, 125, 64, 0, 0, 0, 0, 1135771648, 1133871104)
            #chrs    = [-1, 5, 0, 0, 0, '}', '@', 0, 0, 0, 0, 1135771648, 1133871104]
            #floats  = (nan, 7.006492321624085e-45, 0.0, 0.0, 0.0, 1.7516230804060213e-43, 8.96831017167883e-44, 0.0, 0.0, 0.0, 0.0, 357.0, 299.0)

            #
            #ints    = (-1, 5, 0, 0, 0, 125, 64, 0, 0, 0, 0, 1135771648, 1133871104, 1, 1, 1, 1, 0, -1, 0)
            #dn = 80; data = file_obj.read(dn); n += dn
            #show_ndata(file_obj, 52, types='ifqd', endian='<')
            #show_data(data[4:], types='qd', endian='<')

            #-1, 5, 0, 0, 0, 125, 64, 0, 0, 0, 0,
            parent_zone = 0
            data = file_obj.read(8); n += 8
            strand_id, zone_type = unpack('<2i', data)
            assert zone_type in {0, 1, 2, 3, 4, 5, 6, 7}, zone_type
            log.debug(f'strand_id={strand_id} zone_type={zone_type}')

            data = file_obj.read(12); n += 12
            a, b, c = unpack('<3i', data)
            assert (a, b, c) == (0, 0, 0), (a, b, c) # POINT_FEBRICK_3D_02.plt
            log.debug('abc = (0, 0, 0)')

            data = file_obj.read(8); n += 8
            nnodes, nelement = unpack('<2i', data)
            log.debug(f'nnodes={nnodes} nelement={nelement}')

            data = file_obj.read(16); n += 16
            a, b, c, d = unpack('<4i', data)
            assert (a, b, c, d) == (0, 0, 0, 0), (a, b, c, d) # POINT_FEBRICK_3D_02.plt
            #show(file_obj, 24, types='if')
            log.debug('abcd = (0, 0, 0, 0)')

            # POINT_FEBRICK_3D_02.plt
            #assert zone_type == 5, zone_type # POINT_FEBRICK_3D_02.plt
            celldim = nnodes
            n_misc_neighbor_connections = nnodes

            specify_var = 0
            data_packing = 0  # BLOCK
            #assert data_packing in {0, 1}, data_packing # ???; POINT_FEBRICK_3D_02.plt

            #file.read(4); n += 4
            #show(file_obj, 92, types='dq')
            #n_misc_neighbor_connections = 0
            raw_local = 0
            solution_time = np.nan
            #aaa
        elif version == b'112':
            data = file_obj.read(32); n += 32
            (parent_zone, strand_id, solution_time,
             unused_a, zone_type, data_packing, specify_var) = unpack('<iid 4i', data)
            print(f'  parent_zone={parent_zone} strand_id={strand_id} solution_time={solution_time:g}')
            print(f'  zone_type={zone_type} data_packing={data_packing} specify_var={specify_var}')

            ordered_zone = is_ordered_zone(zone_type)
            fe_poly = is_fe_poly_zone(zone_type)
            #specify_var = 1
            if specify_var:
                ndatai = 4 * nvars; data = file_obj.read(ndatai); n += ndatai
                specify_varsi = np.frombuffer(data, dtype='int32')

            data = file_obj.read(8); n += 8
            raw_local, n_misc_neighbor_connections = unpack('<2i', data)
            assert raw_local in {0, 1}, raw_local
            assert n_misc_neighbor_connections >= 0, n_misc_neighbor_connections
            print(f'  raw_local={raw_local} n_misc_neighbor_connections={n_misc_neighbor_connections}')

            # Are raw local 1-to-1 face neighbors supplied?
            # (0=FALSE 1=TRUE).
            if is_fe_zone(zone_type):
                show_ndata(file_obj, 100, types='ifs', endian='<')
            if ordered_zone:
                raise RuntimeError('ordered zone; p.153')
            if fe_poly:
                raise RuntimeError(' FEPOLYGON or FEPOLYHEDRON; p.153')

            data = file_obj.read(16); n += 16
            nelement, *celldim = unpack('<4i', data)
            print(f'  nelement={nelement} celldim={celldim}')
            data = file_obj.read(4); n += 4
        else:
            raise NotImplementedError(version)
        assert zone_type in {0, 1, 2, 3, 4, 5, 6, 7}, zone_type
        assert -2 <= strand_id < 32700, strand_id
        assert data_packing in {0, 1}, data_packing
        assert nelement > 0, nelement

        assert zone_type in {2, 3, 4, 5}, zone_type
        #assert data_packing == 0, data_packing
        assert specify_var in {0, 1}, specify_var

        #Nodes=41849, Elements=80778, ZONETYPE=FEQuadrilateral

        # Data packing.
        # 0 = Block
        # 1 = Point
        data_packing_str = 'BLOCK' if data_packing == 0 else 'POINT'

        zone_type_str = zone_type_int_to_str(zone_type)

        zone = ZoneTuple(
            zone_name=zone_name,
            strand_id=strand_id,
            solution_time=solution_time,
            zone_type=zone_type,
            zone_type_str=zone_type_str,
            data_packing=data_packing,
            data_packing_str=data_packing_str,
            raw_local=raw_local,
            n_misc_neighbor_connections=n_misc_neighbor_connections,
            celldim=celldim,
            nelement=nelement)
        zones.append(zone)

        # ----------------------------------------------
        #print()
    return n, zones

def zone_type_int_to_str(zone_type: int) -> str:
    if zone_type == 2:
        zone_type_str = 'FETRIANGLE'
    elif zone_type == 3:
        zone_type_str = 'FEQUADRILATERAL'
    elif zone_type == 4:
        zone_type_str = 'FETETRAHEDRON'
    elif zone_type == 5:
        zone_type_str = 'FEBRICK'
    else:
        raise NotImplementedError(zone_type)
    return zone_type_str

def is_ordered_zone(zone_type: int) -> bool:
    """
    0=ORDERED***
    1=FELINESEG
    2=FETRIANGLE
    3=FEQUADRILATERAL
    4=FETETRAHEDRON
    5=FEBRICK
    6=FEPOLYGON
    7=FEPOLYHEDRON
    """
    ordered_zone = (zone_type == 0)
    return ordered_zone

def is_fe_poly_zone(zone_type: int) -> bool:
    """
    0=ORDERED
    1=FELINESEG
    2=FETRIANGLE
    3=FEQUADRILATERAL
    4=FETETRAHEDRON
    5=FEBRICK
    6=FEPOLYGON***
    7=FEPOLYHEDRON***
    """
    fe_poly = (zone_type in {6, 7})
    return fe_poly

def is_fe_zone(zone_type: int) -> bool:
    """
    0=ORDERED
    1=FELINESEG***
    2=FETRIANGLE***
    3=FEQUADRILATERAL***
    4=FETETRAHEDRON***
    5=FEBRICK***
    6=FEPOLYGON***
    7=FEPOLYHEDRON***
    """
    fe_poly = (zone_type in {1, 2, 3, 4, 5, 6, 7})
    return fe_poly

def _write_string(title: str) -> bytes:
    """we convert to an A to 65 using the ascii system"""
    ints = []

    # verify we're in ascii, so we don't get ord bugs
    title_bytes = title.encode('ascii')

    for char in title:
        ints.append(ord(char))
    ints.append(0)  # end of string flag
    nchars = len(ints)

    data = pack(f'<{nchars}i', *ints)
    #print(f'string_ints = {title} {ints}')
    return data

def _write_binary_header(model: TecplotBinary,
                         teplot_file: BinaryIO,
                         version: str) -> None:
    # write header
    assert version == '102', version
    version_bytes = version.encode('ascii')
    version_marker = b'#!TDV' + version_bytes
    assert version_marker in {b'#!TDV102', b'#!TDV111', b'#!TDV112', }, version_marker

    byte_order = 1 # little endian
    data = pack('<8si', version_marker, byte_order)
    teplot_file.write(data)
    teplot_file.write(_write_string(model.title))

    #header_dict['title'] = title
    #header_dict['variables'] = var_names

    nvars = len(model.variables)
    data = pack(b'<i', nvars)
    teplot_file.write(data)
    assert len(model.variables) > 0
    for var in model.variables:
        teplot_file.write(_write_string(var))

def _read_binary_header(model: TecplotBinary,
                        file_obj: BinaryIO) -> tuple[bytes, CaseInsensitiveDict,
                                                     str, str,
                                                     list[str], list[ZoneTuple]]:
    #version, header_dict, title, file_type, var_names, zones
    """
    HEADER SECTION
    The header section contains: the version number of the file, a title
    of the file, the names of the variables to be plotted, the
    descriptions of all zones to be read in and all text and geometry
    definitions.
    """
    self = model
    log = self.log
    data = file_obj.read(8); self.n += 8

    # i. Magic number, Version number
    version_marker = unpack('<8s', data)[0]
    #b'#!TDV102',  # Tecplot 10; 2003
    #b'#!TDV111',  # Tecplot 360; 2008  <--- not supported
    #b'#!TDV112',  # Tecplot 360; 2015 Fun3d
    #b'#!TDV191',  # Tecplot 360; 2020  <--- not supported
    assert version_marker in {b'#!TDV102', b'#!TDV111', b'#!TDV112', }, version_marker

    version = version_marker[5:]
    assert version in {b'102', b'112'}, version
    log.info(f'version = {version!r}')

    # ii. Integer value of 1.
    data = file_obj.read(4); self.n += 4
    byte_order = unpack('<i', data)[0]
    assert byte_order == 1, byte_order

    if version == b'102':  # Tecplot 10; 2003
        file_type_int = 0
    else:
        # iii. Title and variable names.
        # FileType:
        #   0 = FULL,
        #   1 = GRID,
        #   2 = SOLUTION
        data = file_obj.read(4); self.n += 4
        file_type_int = unpack('<i', data)[0]
        assert file_type_int == 0, file_type_int

    if file_type_int == 0:
        file_type = 'full'
    else:
        raise NotImplementedError(file_type_int)

    title, n = _read_string(file_obj, self.n)
    log.info(f'title = {title!r}')
    self.n = n

    # Number of variables (NumVar) in the datafile.
    data = file_obj.read(4); self.n += 4
    nvars = unpack('<i', data)[0]
    assert 3 <= nvars <= 100, nvars
    print(f'nvars = {nvars}')
    model.log.info(f'nvars = {nvars}')

    # Variable names.
    # N = L[1] + L[2] + .... L[NumVar]
    #where:
    # L[i] = length of the ith variable name + 1
    # (for the terminating 0 value).
    var_names = []
    for unused_ivar in range(nvars):
        var_name, n = _read_string(file_obj, self.n)
        self.n = n
        var_names.append(var_name)
    log.info('var_names = %s' % var_names)
    assert self.n == self.f.tell()

    #------------------------------------------
    # iv. Zones
    #Zone marker. Value = 299.0
    n, zones = _read_binary_zone_headers(
        self, self.f, self.n, nvars, version)
    self.n = n

    data = file_obj.read(4); self.n += 4
    flag = unpack('f', data)[0]

    if flag == 299.0:
        raise RuntimeError('geometry')
        # geometry
        y = 1
        # INT32
        # Position CoordSys
        #   0=Grid, 1=Frame,
        #   2=FrameOffset(not used),
        #   3= OldWindow(not used),
        #   4=Grid3D
        #
        # INT32
        # Scope 0=Global 1=Local
        #
        # INT32
        # DrawOrder 0=After, 1=Before
        #
        # FLOAT64*3
        # (X or Theta),(Y or R),(Z or dummy)
        # i.e. the starting location
        #
        # INT32
        # Zone (0=all)
        #
        # INT32
        # Color
        #
        # INT32
        # FillColor
        #
        # INT32
        # IsFilled (0=no 1=yes)
        #
        # INT32
        # GeomType 0=Line, 1=Rectangle 2=Square,
        # 3=Circle, 4=ellipse
        #
        # INT32
        # LinePattern 0=Solid 1=Dashed 2=DashDot
        # 3=DashDotDot 4=Dotted
        # 5=LongDash
        #
        # FLOAT64
        # Pattern Length
        #
        #
        # FLOAT64
        # Line Thickness
        ndata = (
            4 * 4 + 3 * 8 +
            4 * 5 + 8
        )
        data = file_obj.read(ndata); self.n += ndata
        (position_coord, scope, draw_order, x, y, z, zone,
        color, fillcolor, is_filled, geom_type, linepattern, pattern_length,
        ) = unpack('<3idddi 5id ', data)
        del position_coord, scope, draw_order, x, y, z, zone
        del color, fillcolor, is_filled, geom_type, linepattern, pattern_length
    elif flag == 357.0:
        # EOHMARKER (end of header marker)
        pass
    else:
        self.show_ndata(100, types='if')
        raise RuntimeError(flag)

    header_dict = CaseInsensitiveDict()
    header_dict['title'] = title
    header_dict['variables'] = var_names
    return version, header_dict, title, file_type, var_names, zones


def _get_min_max(results: list[np.ndarray]) -> list[np.ndarray]:
    min_max = []
    for result in results:
        mini = result.min()
        maxi = result.max()
        min_max.append(np.array([mini, maxi]))
    return min_max

def _read_binary_min_max(file_obj: BinaryIO, nvars: int) -> tuple[np.ndarray, int]:
    """
    Compressed list of min/max pairs for each non-shared and
    non-passive variable.  For each non-shared and non-passive
    variable (as specified previously)

    TODO: support shared/non-passive variables
    """
    nvalues_min_max = nvars * 2 # nnodes
    min_max_data, ndata = _load_binary_data(file_obj, nvalues_min_max,
                                            data_fmt='float64')

    min_max_data = min_max_data.reshape((2, nvars)).T
    return min_max_data, ndata

def _read_string(file_obj: BinaryIO, n: int) -> tuple[str, int]:
    """
    The letter 'A' has an ASCII value of 65. The WORD
       written to the data file for the letter 'A' is then
       65. In fortran this could be done by doing the following:
    All character strings are null terminated
    (i.e. terminated by a zero value)
    """
    buffer = 20 * 4
    fmt = Struct('<20i')
    all_ints = []
    while 1:
        data = file_obj.read(buffer)
        ints = fmt.unpack(data)
        #print('ints =', ints)
        try:
            i0 = ints.index(0)
        except ValueError:
            i0 = -1

        if i0 == -1:
            n += buffer
            all_ints.extend(ints)
        else:
            n += 4 * i0 + 4  # we have a +4 to skip the null byte
            all_ints.extend(ints[:i0])  # we leave off the +1 to skip the null byte
            file_obj.seek(n)
            break

    chars = [chr(val) for val in all_ints]
    string = ''.join(chars)
    #print('string =', string)
    return string, n

def zones_to_exclude_to_set(zones_to_exclude: Optional[list[int]]=None) -> set[int]:
    """
    Simplifies zone exclusion

        Parameters
    ----------
    zones_to_exclude : list[int]; default=None -> []
        0-based list of zones to exlcude

    Returns
    -------
    set_zones_to_exclude : set[int]
        0-based list of zones to exlcude

    """
    if zones_to_exclude is None:
        set_zones_to_exclude = set([])
    else:
        set_zones_to_exclude = set(zones_to_exclude)
    return set_zones_to_exclude

def show(file_obj: BinaryIO, n: int, types: str='ifs', endian: str='<'):
    """not safe from positional errors"""
    nold = file_obj.tell()
    nints = n // 4
    data = file_obj.read(4 * nints)
    strings, ints, floats = show_data(data, types=types, endian=endian)
    file_obj.seek(nold)
    return strings, ints, floats

def show_data(data: bytes, types: str='ifs', endian: str='<'):
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
    return _write_data(sys.stdout, data, types=types, endian=endian)

def _write_data(f, data: bytes, types: str='ifs', endian: str='<'):
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
        h/H - signed/unsigned short (16 bit int)
    endian : str; default=None -> auto determined somewhere else in the code
        the big/little endian {>, <}

    """
    n = len(data)
    nshorts = n // 2
    nints = n // 4
    ndoubles = n // 8
    strings = None
    ints = None
    floats = None
    longs = None

    assert endian is not None, endian
    f.write('\nndata = %s:\n' % n)
    for typei in types:
        assert typei in 'sifdq lIL hH', 'type=%r is invalid' % typei

    if 's' in types:
        strings = unpack('%s%is' % (endian, n), data)
        f.write("  strings = %s\n" % str(strings))
    if 'i' in types and nints > 0:
        ints = unpack('%s%ii' % (endian, nints), data)
        chrs = [chr(val) if 32 <= val <= 130 else val
                for val in ints]
        f.write("  ints    = %s\n" % str(ints))
        f.write("  chrs    = %s\n" % str(chrs))
    if 'f' in types and nints > 0:
        floats = unpack('%s%if' % (endian, nints), data)
        f.write("  floats  = %s\n" % str(floats))
    if 'd' in types and ndoubles > 0:
        doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
        f.write("  doubles (float64) = %s\n" % str(doubles))

    if 'l' in types and nints > 0:
        longs = unpack('%s%il' % (endian, nints), data)
        f.write("  long  = %s\n" % str(longs))
    if 'I' in types and nints > 0:
        ints2 = unpack('%s%iI' % (endian, nints), data)
        f.write("  unsigned int = %s\n" % str(ints2))
    if 'L' in types and nints > 0:
        longs2 = unpack('%s%iL' % (endian, nints), data)
        f.write("  unsigned long = %s\n" % str(longs2))
    if 'q' in types and ndoubles > 0:
        longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
        f.write("  long long (int64) = %s\n" % str(longs))

    if 'h' in types and nshorts > 0:
        signed_short = unpack('%s%ih' % (endian, nshorts), data)
        f.write("  signed short (int16) = %s\n" % str(signed_short))
    if 'H' in types and nshorts > 0:
        unsigned_short = unpack('%s%iH' % (endian, nshorts), data)
        f.write("  unsigned short (int16) = %s\n" % str(unsigned_short))
    f.write('\n')
    return strings, ints, floats

def show_ndata(file_obj: BinaryIO, n: int, types: str='ifs', endian: str='<'):
    return _write_ndata(file_obj, sys.stdout, n, types=types, endian=endian)

def _write_ndata(file_obj: BinaryIO, f, n: int, types: str='ifs', endian: str='<'):
    """
    Useful function for seeing what's going on locally when debugging.
    """
    nold = file_obj.tell()
    data = file_obj.read(n)
    file_obj.seek(nold)
    return _write_data(f, data, types=types, endian=endian)

