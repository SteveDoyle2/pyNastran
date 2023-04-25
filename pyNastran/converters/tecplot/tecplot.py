"""
models from:
    http://people.sc.fsu.edu/~jburkardt/data/tec/tec.html
"""
import os
import copy
import itertools
from io import StringIO
from pathlib import PurePath
from typing import TextIO, Union, Optional, Any

import numpy as np
from cpylog import SimpleLogger

from pyNastran.utils import is_binary_file
from pyNastran.converters.tecplot.zone import Zone
from pyNastran.converters.tecplot.tecplot_binary import (
    TecplotBinary, zones_to_exclude_to_set)
from pyNastran.converters.tecplot.write_ascii import (
    write_ascii_header, write_ascii_tecplot_zone)
from pyNastran.converters.tecplot.read_ascii import (
    read_header_lines, header_lines_to_header_dict,
    read_zonetype,
)


def read_tecplot(tecplot_filename: str,
                 use_cols=None, dtype=None,
                 filetype: str='guess',
                 zones_to_exclude: Optional[list[int]] = None,
                 log=None, debug=False):
    """loads a tecplot file

    Parameters
    ----------
    zones_to_exclude : list[int]; default=None -> []
        0-based list of zones to exlcude
    """
    tecplot = Tecplot(log=log, debug=debug)
    if use_cols:
        tecplot.use_cols = use_cols
        tecplot.dtype = dtype
    tecplot.read_tecplot(tecplot_filename,
                         filetype=filetype,
                         zones_to_exclude=zones_to_exclude)
    return tecplot


class Tecplot(TecplotBinary):
    """
    Parses a binary/ASCII Tecplot 360 file.
    Writes a binary/ASCII Tecplot 10 file.

    Supports:
     - title
     - single zone only?
     - binary unstructured:
       - nodal results
       - single element type;
            ZONETYPE = [FETRIANGLE, FEQUADRILATERAL,
                        FETETRAHEDRON, FEBRICK]
       - DATAPACKING = [BLOCK] writing
       - 3D
       - full support for writing
     - ASCII unstructured:
       - nodal results
       - single element type;
           ZONETYPE = [FETRIANGLE, FEQUADRILATERAL,
                       FETETRAHEDRON, FEBRICK]
       - DATAPACKING = [POINT, BLOCK] writing
       - 2D/3D
       - full support for writing
     - ASCII structured:
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

    @property
    def nzones(self):
        """gets the number of zones"""
        return len(self.zones)

    def read_tecplot(self, tecplot_filename: str,
                     filetype: str='guess',
                     zones_to_exclude: Optional[list[int]] = None):
        """
        Reads an ASCII/binary Tecplot file.

        The binary file reader must have ONLY CHEXAs and be Tecplot 360 with
        `rho`, `u`, `v`, `w`, and `p`.
        The ASCII file reader has only been tested with Tecplot 10, but will
        probably work on Tecplot360.  It **should** work with any set of
        variables.
        """
        filetype = filetype.lower()
        assert filetype in ['guess', 'ascii', 'binary'], filetype
        if filetype == 'binary' or (filetype == 'guess' and is_binary_file(tecplot_filename)):
            return self.read_tecplot_binary(tecplot_filename, zones_to_exclude=zones_to_exclude)
        return self.read_tecplot_ascii(tecplot_filename, zones_to_exclude=zones_to_exclude)

    def read_tecplot_ascii(self, tecplot_filename: Union[str, PurePath, StringIO],
                           nnodes=None, nelements=None,
                           zones_to_exclude: Optional[list[int]]=None):
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
        set_zones_to_exclude = zones_to_exclude_to_set(zones_to_exclude)
        self.tecplot_filename = tecplot_filename
        iline = 0
        nnodes = -1
        nelements = -1
        if isinstance(tecplot_filename, StringIO):
            lines = tecplot_filename.readlines()
        else:
            assert os.path.exists(tecplot_filename), tecplot_filename
            with open(tecplot_filename, 'r') as tecplot_file:
                lines = tecplot_file.readlines()

        quads_list: list[np.ndarray] = []
        hexas_list: list[np.ndarray] = []
        tris_list: list[np.ndarray] = []
        tets_list: list[np.ndarray] = []
        zone_data_list: list[np.ndarray] = []

        line = lines[iline].strip()
        iline += 1
        iblock = 0
        izone = 0
        while 1:
            #print('start...')
            iline, title_line, header_lines, line = read_header_lines(
                lines, iline, line, self.log)
            #print('header_lines', header_lines)
            headers_dict = header_lines_to_header_dict(
                title_line, header_lines, self.variables, self.log)
            if headers_dict is None:
                break

            zone = Zone(self.log)
            zone.headers_dict = copy.deepcopy(headers_dict)
            self.variables = headers_dict['VARIABLES']
            if 'NAME' in headers_dict:
                zone.name = headers_dict['NAME']
            #print('self.variables', self.variables)

            #print(headers_dict.keys())
            if 'ZONETYPE' in headers_dict:
                zone_type = zone.zonetype  # febrick
                data_packing = headers_dict['DATAPACKING'].upper() # block
                iline = read_zonetype(
                    self.log,
                    zone, zone_type, lines, iline, iblock, headers_dict, line,
                    nnodes, nelements,
                    zone_data_list, hexas_list, tets_list, quads_list, tris_list,
                    data_packing=data_packing)
                iline -= 1
            elif 'DATAPACKING' in headers_dict:
                data_packing = headers_dict['DATAPACKING'] # FEPoint
                assert isinstance(data_packing, str), headers_dict
                zone_type = data_packing.upper() # FEPoint
                self.log.debug('zone_type = %r' % zone_type[0])
                iline = read_zonetype(
                    self.log,
                    zone, zone_type, lines, iline, iblock, headers_dict, line,
                    nnodes, nelements,
                    zone_data_list, hexas_list, tets_list, quads_list, tris_list,
                    fe=data_packing)
                iline -= 1
            elif (('ZONE' in headers_dict) and
                  (headers_dict['ZONE'] is None) and
                  ('T' in headers_dict)):
                lines2 = itertools.chain((line, ), iter(lines[iline:]))
                A = self._read_table_from_lines(lines2, headers_dict)
                self.A = A
                self.log.debug(f'read_table; A.shape={A.shape}...')
                asdf
                return
            else:
                msg = 'Expected ZONETYPE, F, or "ZONE T"\n'
                msg += 'headers=%s\n' % str(headers_dict)
                msg += 'line = %r' % line.strip()
                raise NotImplementedError(msg)

            #print(zone)
            if izone in set_zones_to_exclude:
                log.warning(f'skipping izone={izone}')
            else:
                self.zones.append(zone)

                #sline = line.split()
                #print('stack...')
                _stack(zone, zone_data_list,
                       quads_list, tris_list,
                       tets_list, hexas_list,
                       self.log)

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
            zone_data_list = []
            izone += 1

    @property
    def nnodes(self) -> int:
        nnodes = sum([zone.nnodes for zone in self.zones])
        return nnodes

    def stack_geometry(self) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                      np.ndarray, np.ndarray, np.ndarray,
                                      list[str]]:
        nnodes = self.nnodes
        nodes = np.zeros((nnodes, 3), dtype='float32')
        inode = 0
        tris = []
        quads = []
        tets = []
        hexas = []
        zone_ids = []
        names = []
        for izone, zone in enumerate(self.zones):
            names.append(zone.name)
            nodes2d = zone.xy
            nodes3d = zone.xyz
            nnodes2d = nodes2d.shape[0]
            nnodes3d = nodes3d.shape[0]
            #assert nnodes2d == 0, zone

            # elements
            quadsi = np.zeros((0, 4), dtype='int32')
            hexasi = np.zeros((0, 8), dtype='int32')
            trisi = np.zeros((0, 3), dtype='int32')
            tetsi = np.zeros((0, 4), dtype='int32')
            if 'I' in zone.headers_dict:
                i = zone.headers_dict['I']
                if 'J' in zone.headers_dict:
                    j = zone.headers_dict['J']
                    if 'K' in zone.headers_dict:
                        k = zone.headers_dict['K']
                        nnodes = i * j * k
                        elements = np.arange(0, nnodes).reshape(k, j, i)
                        #print(elements[0, :, :])
                        #print(elements[1:, :, :])
                        n1 = elements[:-1, :-1, :-1].ravel()
                        n2 = elements[:-1, :-1, 1:].ravel()
                        n3 = elements[:-1, 1:, 1:].ravel()
                        n4 = elements[:-1, 1:, :-1].ravel()

                        n5 = elements[1:, :-1, :-1].ravel()
                        n6 = elements[1:, :-1, 1:].ravel()
                        n7 = elements[1:, 1:, 1:].ravel()
                        n8 = elements[1:, 1:, :-1].ravel()
                        #nhexas = (i - 1) * (j - 1) * (k - 1)
                        hexasi = np.vstack([n1, n2, n3, n4, n5, n6, n7, n8]).T
                    else:
                        nnodes = i * j
                        elements = np.arange(0, nnodes).reshape(j, i)
                        #print('elements:')
                        #print(elements)
                        n1 = elements[:-1, :-1].ravel()
                        n2 = elements[:-1, 1:].ravel()
                        n3 = elements[1:, 1:].ravel()
                        n4 = elements[1:, :-1].ravel()
                        #nquads = (i - 1) * (j - 1)
                        quadsi = np.vstack([n1, n2, n3, n4]).T
                    #trisi = None
                    #tetsi = None
                    #hexasi = None
            else:
                quadsi = zone.quad_elements
                hexasi = zone.hexa_elements
                tetsi = zone.tet_elements
                trisi = zone.tri_elements

            nquadsi = len(quadsi)
            ntrisi = len(trisi)
            ntetsi = len(tetsi)
            nhexasi = len(hexasi)
            nelementsi = nquadsi + ntrisi + ntetsi + nhexasi
            assert nelementsi > 0, str(zone)
            if nquadsi:
                quads.append(inode + quadsi)
            if ntrisi:
                tris.append(inode + trisi)
            if ntetsi:
                tets.append(inode + tetsi)
            if nhexasi:
                hexas.append(inode + hexasi)

            # nodes
            if nnodes2d and nnodes3d:
                raise RuntimeError('2d and 3d nodes is not supported\n'
                                   f'name={zone.name!r} zone.xy={zone.xy} zone.xyz={zone.xyz}')
            elif nnodes2d:
                nodes[inode:inode+nnodes2d, :2] = nodes2d
                inode += nnodes2d
            elif nnodes3d:
                nodes[inode:inode+nnodes3d] = nodes3d
                inode += nnodes3d
            else:  # pragma: no cover
                raise RuntimeError('failed to find 2d/3d nodes')
            #print('inode', inode)
            #print(nodes)
            #print('-------------')
            zone_idsi = np.ones(nelementsi) * izone
            zone_ids.append(zone_idsi)
        #print('stack', len(quads))
        quads = stack(quads)
        tris = stack(tris)
        if len(quads) and not len(tris):
            # convert quads written as tris into tris
            itri = (quads[:, 2] == quads[:, 3])
            tris = quads[itri, :3]
            quads = quads[~itri, :]

        tets = stack(tets)
        hexas = stack(hexas)
        zone_ids = np.hstack(zone_ids).astype('int32')
        return nodes, tris, quads, tets, hexas, zone_ids, names

    def stack_results(self) -> np.ndarray:
        results = []
        for zonei in self.zones:
            results.append(zonei.nodal_results)
        results_array = np.vstack(results)
        return results_array

    def read_table(self, tecplot_file: TextIO,
                   unused_iblock: int,
                   headers_dict: dict[str, Any],
                   line: str) -> tuple[np.ndarray, Any]:
        """
        reads a space-separated tabular data block
        """
        # add on the preceding line to the line "list"
        # that's not a hack at all...
        lines = itertools.chain((line, ), iter(tecplot_file))
        A, blank = self._read_table_from_lines(lines, headers_dict)
        return A, None

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

    def _slice_plane(self, zone: Zone, y: np.ndarray, slice_value: float) -> None:
        """
        - Only works for CHEXA
        - Doesn't remove unused nodes/renumber elements
        """
        slice_value = float(slice_value)
        inodes = np.where(y < slice_value)[0]
        zone._slice_plane_inodes(inodes)

    def write_tecplot(self, tecplot_filename: str,
                      res_types: list[str]=None,
                      adjust_nids: bool=True) -> None:
        """
        Only handles single type writing

        Parameters
        ----------
        tecplot_filename : str
            the path to the output file
        res_types : str; list[str, str, ...]; default=None -> all
            the results that will be written (must be consistent with
            self.variables)
        adjust_nids : bool; default=True
            element_ids are 0-based in binary and must be switched to
            1-based in ASCII

        """
        self.write_tecplot_ascii(tecplot_filename, res_types, adjust_nids)

    def write_tecplot_ascii(self, tecplot_filename: str,
                            res_types: list[str]=None,
                            adjust_nids: bool=True) -> None:
        """writes an ASCII tecplot file"""
        self.log.info('writing tecplot %s' % tecplot_filename)
        msg, ivars = _get_write_header(self.title, self.zones, res_types)
        with open(tecplot_filename, 'w') as tecplot_file:
            tecplot_file.write(msg)
            for zone in self.zones:
                write_ascii_tecplot_zone(
                    tecplot_file, zone, ivars,
                    self.log, adjust_nids)

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

def _get_write_header(title: str,
                      zones: list[Zone],
                      res_types: list[str]) -> tuple[str, np.ndarray]:
    """gets the tecplot header"""
    is_x = True
    is_y = True
    is_z = True
    #is_results = False
    nzones = len(zones)
    assert nzones >= 1, nzones
    for zone in zones:
        variables = zone.variables
        is_x = 'X' in zone.headers_dict['VARIABLES']
        is_y = 'Y' in zone.headers_dict['VARIABLES']
        is_z = 'Z' in zone.headers_dict['VARIABLES']
        res_types_all = []
        if is_x:
            res_types_all.append('X')
        if is_y:
            res_types_all.append('Y')
        if is_z:
            res_types_all.append('Z')

        if res_types is None:
            res_types = zone.variables
        elif isinstance(res_types, str):
            res_types = [res_types]

        for res in res_types:
            if res not in res_types_all:
                res_types_all.append(res)
        break

    assert len(variables), variables
    msg, ivars = write_ascii_header(
        title, is_x, is_y, is_z,
        res_types_all, variables)
    return msg, ivars


def stack(elements: list[np.ndarray]) -> np.ndarray:
    if len(elements) == 0:
        pass
    elif len(elements) == 1:
        elements = elements[0]
        #print(elements)
    else:
        #print('----stack------')
        #for elementsi in elements:
            #print(elementsi)
        #print('----stack------')
        elements = np.vstack(elements)
    return elements

def _stack(zone: Zone,
           zone_data_list: list[np.ndarray],
           quads_list: list[np.ndarray], tris_list: list[np.ndarray],
           tets_list: list[np.ndarray], hexas_list: list[np.ndarray],
           log: SimpleLogger):
    """
    elements are read as a list of lines, so we need to stack them
    and cast them while we're at it.
    """
    log.debug('stacking elements')
    if len(hexas_list):
        zone.hexa_elements = np.vstack(hexas_list)
        zone.headers_dict['zonetype'] = 'FEBRICK'
    if len(tets_list):
        zone.tet_elements = np.vstack(tets_list)
        zone.headers_dict['zonetype'] = 'FETETRAHEDRON'
    if len(quads_list):
        zone.quad_elements = np.vstack(quads_list)
        zone.headers_dict['zonetype'] = 'FEQUADRILATERAL'
    if len(tris_list):
        zone.tri_elements = np.vstack(tris_list)
        zone.headers_dict['zonetype'] = 'FETRIANGLE'
    #if 'zonetype' not in zone.headers_dict:
        #raise RuntimeError('no elements')

    #log.debug('stacking nodes')
    #if len(xyz_list) == 1:
        #xyz = xyz_list[0]
    #else:
        #xyz = np.vstack(xyz_list)

    #self.elements = elements - 1
    #print(self.elements)
    if zone_data_list:
        zone.zone_data = np.vstack(zone_data_list)
        #print('**', len(zone_data_list), 'shape=', zone.zone_data.shape)

    if zone.is_3d:
        zone.log.info('3d zone')
        zone.xyz
        zone.nodal_results
    elif zone.is_2d:
        zone.log.info('2d zone')
        zone.xy
    else:
        zone.nodal_results
    x = 1

def main():  # pragma: no cover
    #plt = Tecplot()
    fnames = os.listdir(r'C:\output\time20000')

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
    fnames = [os.path.join(r'C:\output\time20000', fname)
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
        #except Exception:
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
        plt.write_tecplot(f'processor_{iprocessor:d}.plt')

def main3():  # pragma: no cover
    import sys
    tecplot_filename = sys.argv[1]
    plt = Tecplot()
    plt.read_tecplot(tecplot_filename)
    tecplot_filename_out = tecplot_filename + '.out'
    plt.write_tecplot(tecplot_filename_out, res_types=None, adjust_nids=False)
    x = 1

if __name__ == '__main__':   # pragma: no cover
    main3()
