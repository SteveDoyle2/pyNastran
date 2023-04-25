from __future__ import annotations
import copy
from collections import defaultdict
from typing import TextIO, Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.nptyping_interface import NDArrayN3float, NDArrayN3int, NDArrayN4int

#if TYPE_CHECKING:  # pragma: no cover
from cpylog import SimpleLogger

class CaseInsensitiveDict(dict):
    def __contains__(self, key: str) -> bool:
        val_in = dict.__contains__(self, key.upper())
        return val_in
    def __getitem__(self, key) -> Any:
        val = dict.__getitem__(self, key.upper())
        #log.info("GET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        return val
    def __setitem__(self, key, val) -> None:
        #log.info("SET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        dict.__setitem__(self, key.upper(), val)


class Zone:
    def __init__(self, log: SimpleLogger):
        self.log = log

        self.name = '???'
        self.strand_id = -1
        self.headers_dict = CaseInsensitiveDict()
        self.zone_data = np.zeros((0, 0), dtype='float32')
        self.tet_elements = np.array([], dtype='int32')
        self.hexa_elements = np.array([], dtype='int32')
        self.quad_elements = np.array([], dtype='int32')
        self.tri_elements = np.array([], dtype='int32')

        # TODO: does or does not consider xyz?
        #self.variables: list[str] = []

        # VARLOCATION=([3-7,10]=CELLCENTERED, [11-12]=CELLCENTERED)
        # default=NODAL
        #self.nodal_results = np.array([], dtype='float32')
        self.A = None
        #self.centroidal_results = np.array([], dtype='float32')

    @property
    def nxyz(self) -> int:
        nxyz = 0
        for var in self.headers_dict['variables']:
            if var.lower() in ['x', 'y', 'z']:
                nxyz += 1
        return nxyz

    @classmethod
    def set_zone_from_360(self,
                          log: SimpleLogger,
                          header_dict,
                          variables: list[str],
                          name: str,
                          zone_type: str,
                          data_packing: str,
                          strand_id: int,
                          tris=None,
                          quads=None,
                          tets=None,
                          hexas=None,
                          zone_data=None):
        assert isinstance(log, SimpleLogger), log
        assert isinstance(header_dict, dict), header_dict
        assert isinstance(variables, list), variables
        assert isinstance(name, str), name
        assert isinstance(zone_type, str), zone_type
        assert isinstance(data_packing, str), data_packing
        assert zone_type in {'FETRIANGLE', 'FEQUADRILATERAL', 'FETETRAHEDRON', 'FEBRICK'}, zone_type
        assert data_packing in {'POINT', 'BLOCK'}, data_packing
        for variable in variables:
            assert isinstance(variable, str), variable
        variables = copy.deepcopy(variables)

        zone = Zone(log)
        zone.headers_dict = copy.copy(header_dict)
        zone.name = name
        zone.strand_id = strand_id


        if tris is not None:
            zone.tri_elements = tris
        if quads is not None:
            zone.quad_elements = quads

        if tets is not None:
            zone.tet_elements = tets
        if hexas is not None:
            zone.hexa_elements = hexas

        nxyz = zone.nxyz

        nvars = len(variables)
        nvars_actual = zone_data.shape[1]
        if not nvars == nvars_actual:
            raise RuntimeError(f'variables={variables} nvars={nvars}; '
                               f'nodal_results.shape={nodal_results.shape} nxyz={nxyz}\n'
                               f'nvars={nvars} nvars_actual={nvars_actual}')

        nnodes = zone_data.shape[0]
        if self.variables_exist(variables, ['rho', 'u', 'v', 'w'], ['rhoUVW']):
            irho = variables.index('rho')
            iu = variables.index('u')
            iv = variables.index('v')
            iw = variables.index('w')

            rho = zone_data[:, irho]
            uvw = zone_data[:, [iu, iv, iw]]

            uvw_mag = np.linalg.norm(uvw, axis=1)
            rho_uvw_mag = (rho * uvw_mag).reshape((nnodes, 1))
            rho_uvw = np.abs(rho[:, np.newaxis] * uvw)
            zone_data = np.hstack([zone_data, rho_uvw, rho_uvw_mag])
            variables.extend(['rhoU', 'rhoV', 'rhoW', 'rhoUVW'])

        #assert self.variables_exist(variables, ['mach', 'tt', 'gamma'], ['ttot'])
        if self.variables_exist(variables, ['mach', 'tt', 'gamma'], ['ttot', 'ptot']):
            itt = variables.index('tt')  # translational-rotational temperature
            ip = variables.index('p')  # pressure
            imach = variables.index('mach')
            igamma = variables.index('gamma')

            tt = zone_data[:, itt]
            p = zone_data[:, ip]
            gamma = zone_data[:, igamma]
            mach = zone_data[:, imach]

            # t0 / t = 1 + (gamma - 1) / 2 * mach ^ 2
            gm1 = gamma - 1
            ktemp = 1 + gm1 / 2 * mach ** 2
            ttot = tt * ktemp
            ptot = p * ktemp ** (gamma / gm1) # https://en.wikipedia.org/wiki/Total_pressure

            zone_data = np.hstack([zone_data,
                                   ttot.reshape((nnodes, 1)),
                                   ptot.reshape((nnodes, 1))])
            variables.extend(['ttot', 'ptot'])
        zone.zone_data = zone_data
        header_dict['variables'] = variables

        #zone.headers_dict = copy.copy(header_dict)
        zone.headers_dict['datapacking'] = data_packing
        zone.headers_dict['zonetype'] = zone_type
        zone.variables = variables

        #print(str(zone))
        return zone

    def variables_exist(variables: list[str],
                        variables_to_check: list[str],
                        variables_to_save: list[str]) -> bool:
        all_exist = all(var in variables for var in variables_to_check)
        if all(var in variables for var in variables_to_save):
            all_exist = True
        return all_exist

    @property
    def is_unstructured(self) -> bool:
        """are there unstructured elements"""
        is_unstructured = (
            len(self.tri_elements) or len(self.quad_elements) or
            len(self.tet_elements) or len(self.hexa_elements)
        )
        return is_unstructured

    @property
    def is_structured(self) -> bool:
        """is the model in plot3d style format"""
        return not self.is_unstructured

    @property
    def nnodes(self) -> int:
        """gets the number of nodes in the model"""
        return self.zone_data.shape[0]

    @property
    def nelements(self) -> int:
        """gets the number of elements in the model"""
        return (self.hexa_elements.shape[0] + self.tet_elements.shape[0] +
                self.quad_elements.shape[0] + self.tri_elements.shape[0])

    @property
    def result_names(self) -> list[str]:
        """gets the variables"""
        return self.variables

    @result_names.setter
    def result_names(self, vals: list[str]) -> None:
        """sets the variables"""
        self.variables = vals

    @property
    def ndim(self) -> int:
        is_x = 'X' in self.variables
        is_y = 'Y' in self.variables
        is_z = 'Z' in self.variables
        return is_x + is_y + is_z

    @property
    def title(self) -> str:
        if 'TITLE' not in self.headers_dict:
            self.headers_dict['TITLE'] = 'tecplot geometry and solution file'
        return self.headers_dict['TITLE']
    @title.setter
    def title(self, title: str) -> None:
        self.headers_dict['TITLE'] = title

    @property
    def is_point(self) -> bool:
        if 'F' in self.headers_dict:
            raise RuntimeError(self.headers_dict)
        if 'DATAPACKING' not in self.headers_dict:
            return True
        is_point = (self.headers_dict['DATAPACKING'] == 'POINT')
        is_block = (self.headers_dict['DATAPACKING'] == 'BLOCK')
        assert is_point or is_block
        return is_point

    @property
    def is_block(self) -> bool:
        return not self.is_point

    @property
    def datapacking(self) -> Optional[str]:
        """'POINT', 'BLOCK'"""
        try:
            return self.headers_dict['DATAPACKING']
        except KeyError:
            return None

    @property
    def zone_type_int(self) -> int:
        """
        0=ORDERED
        1=FELINESEG
        2=FETRIANGLE
        3=FEQUADRILATERAL
        4=FETETRAHEDRON
        5=FEBRICK
        6=FEPOLYGON
        7=FEPOLYHEDRON
        """
        zonetype = self.headers_dict['ZONETYPE']
        if zonetype == 'FETRIANGLE':
            zone_type_int = 2
        elif zonetype == 'FEQUADRILATERAL':
            zone_type_int = 3
        elif zonetype == 'FETETRAHEDRON':
            zone_type_int = 4
        elif zonetype == 'FEBRICK':
            zone_type_int = 5
        else:
            raise RuntimeError(zonetype)
        return zone_type_int

    @property
    def zonetype(self) -> Optional[str]:
        """FEBrick,  FETETRAHEDRON,  FETRIANGLE,  FEQUADRILATERAL"""
        #headers_dict['ZONETYPE'].upper() # FEBrick
        #try:
        return self.headers_dict['ZONETYPE'].upper()
        #except KeyError:
            #return None

    def get_xyz(self) -> np.ndarray:
        """turns 2d points into 3d points"""
        nnodes3d = self.xyz.shape[0]
        nnodes2d = self.xy.shape[0]
        if nnodes2d and nnodes3d:
            raise RuntimeError('2d and 3d nodes is not supported')
        elif nnodes2d:
            npoints = self.xy.shape[0]
            xyz = np.zeros((npoints, 3), dtype=self.xy.dtype)
            xyz[:, :2] = self.xy
        elif nnodes3d:
            xyz = self.xyz
        else:  # pragma: no cover
            raise RuntimeError('failed to find 2d/3d nodes')
        return xyz

    def __repr__(self) -> str:
        name = ''
        if hasattr(self, 'name'):
            name = f'  name = {self.name}\n'

        xy_shape = str(self.xy.shape)
        xyz_shape = str(self.xyz.shape)
        is3d = self.is_3d
        if 'I' in self.headers_dict:
            msgi = self.repr_nijk()
        else:
            msgi = (
                f'  datapacking = {self.datapacking!r}\n'
                f'  zonetype = {self.zone_type_int} ({self.zonetype!r})\n'
            )

        title2 = '  T = %r\n' % self.headers_dict['T'] if 'T' in self.headers_dict else ''
        msg = (
            'Zone:\n'
            #f'  filename = {self.tecplot_filename!r}\n'
            #f'  is_mesh = {self.is_mesh}\n'
            + name +
            f'  variables = {self.variables}\n'
            f'  is3d = {is3d}\n'
            f'    xyz.shape = {xyz_shape}\n'
            f'  nodal_results.shape = {str(self.nodal_results.shape)}\n'
            #f'  headers_dict = {self.headers_dict}\n'
            f'  title = {self.title!r}\n{title2}{msgi}'
            #f'  use_cols = {self.use_cols}\n'
            #f'  dtype = {self.dtype}\n'
        )
        return msg

    def repr_nijk(self) -> str:
        ni = self.headers_dict['I']
        if 'J' in self.headers_dict:
            nj = self.headers_dict['J']
            if 'K' in self.headers_dict:
                nk = self.headers_dict['K']
                msgi = f'  nI={ni} nJ={nj} nK={nk}\n'
            else:
                msgi = f'  nI={ni} nJ={nj}\n'
        else:
            assert 'K' not in self.headers_dict, list(self.headers_dict.keys())
            msgi = f'  nI={ni}\n'
        return msgi

    def write_unstructured_zone(self, tecplot_file: TextIO, ivars: list[int], is_points: bool,
                                nnodes: int, nelements: int, zone_type: int, log: SimpleLogger,
                                is_tris: bool, is_quads: bool, is_tets: bool, is_hexas: bool,
                                adjust_nids: bool=True) -> None:
        msg = 'ZONE '
        self.log.debug(f'is_points = {is_points}')
        datapacking = 'POINT' if  is_points else 'BLOCK'
        msg += f' n={nnodes:d}, e={nelements:d}, ZONETYPE={zone_type}, DATAPACKING={datapacking}\n'
        tecplot_file.write(msg)

        self._write_xyz_results(tecplot_file, is_points, ivars)
        self._write_elements(tecplot_file, nnodes,
                             is_tris, is_quads, is_tets, is_hexas,
                             adjust_nids=adjust_nids)

    def determine_element_type(self) -> tuple[bool, bool, bool, Optional[str],
                                              bool, bool, bool, bool]:
        """
        Returns
        -------
        is_structured : bool
           is this a structured grid
        is_unstructured : bool
           is this an unstructured grid
        is_points : bool
           is the zone in POINT/BLOCK format
        zone_type : str / None
           str : FEBrick,  FETETRAHEDRON,  FETRIANGLE,  FEQUADRILATERAL
           None : ???
        is_tris : bool
            are there CTRIA3s
        is_quads : bool
            are there CQUAD4s
        is_tets : bool
            are there CTETRAs
        is_hexas : bool
            are there CHEXAs

        """
        etype_elements = [
            ('HEXA', self.hexa_elements),
            ('TETRA', self.tet_elements),
            ('TRIA3', self.tri_elements),
            ('QUAD4', self.quad_elements),
        ]
        is_point = self.is_point
        is_tets = False
        is_hexas = False
        is_tris = False
        is_quads = False

        is_structured = False
        is_unstructured = False
        zone_type = None
        for etype, elements in etype_elements:
            if not len(elements):
                continue
            if etype == 'HEXA' and len(elements):
                # is_points = False
                is_hexas = True
                #nnodes_per_element = 8
                zone_type = 'FEBrick'
            elif etype == 'TETRA' and len(elements):
                # is_points = False
                is_tets = True
                #nnodes_per_element = 4
                zone_type = 'FETETRAHEDRON'
            elif etype == 'TRIA3' and len(elements):
                # is_points = True
                is_tris = True
                #nnodes_per_element = 3
                zone_type = 'FETRIANGLE'
            elif etype == 'QUAD4' and len(elements):
                # is_points = True
                is_quads = True
                #nnodes_per_element = 4
                zone_type = 'FEQUADRILATERAL'
            else:
                self.log.info('etype=%r' % etype)
                self.log.info(elements)
                continue
            break
        if any([is_tris, is_quads, is_tets, is_hexas]):
            is_unstructured = True
        elif 'I' in self.headers_dict:
            is_structured = True
        else:  # pragma: no cover
            raise RuntimeError('is this structured or unstructured?')
        return is_structured, is_unstructured, is_point, zone_type, is_tris, is_quads, is_tets, is_hexas

    def _write_elements(self, tecplot_file: TextIO, nnodes: int,
                        is_tris: bool, is_quads: bool, is_tets: bool, is_hexas: bool,
                        adjust_nids: bool=True) -> None:
        """Writes the unstructured elements.  Verifies that nodes are sequential."""
        self.log.debug('is_hexas=%s is_tets=%s is_quads=%s is_tris=%s' %
                       (is_hexas, is_tets, is_quads, is_tris))
        if is_hexas:
            efmt = ' %d %d %d %d %d %d %d %d\n'
            elements = self.hexa_elements
        elif is_tets:
            efmt = ' %d %d %d %d\n'
            elements = self.tet_elements
        elif is_quads:
            efmt = ' %d %d %d %d\n'
            elements = self.quad_elements
        elif is_tris:
            efmt = ' %d %d %d\n'
            elements = self.tri_elements
        else:
            raise RuntimeError()

        # we do this before the nid adjustment
        node_min = elements.min()
        node_max = elements.max()
        self.log.debug(f'inode: min={node_min:d} max={node_max:d}')
        assert node_min >= 0, node_min

        if node_max > nnodes:
            msg = 'elements.min()=node_min=%d elements.max()=node_max=%d nnodes=%d' % (
                node_min, node_max, nnodes)
            self.log.error(msg)
            raise RuntimeError(msg)
        # assert elements.min() == 1, elements.min()
        # assert elements.max() == nnodes, elements.max()

        if adjust_nids:
            elements = elements + 1

        for element in elements:
            tecplot_file.write(efmt % tuple(element))

    def _write_xyz_results(self, tecplot_file: TextIO,
                           is_points: bool,
                           ivars: list[int]) -> None:
        """writes XY/XYZs and results in POINT or BLOCK format"""
        # xyz
        xyz = self.xyz
        xy = self.xy

        #print('nnodes', nodes.shape)
        #print('results', self.nodal_results.shape)
        assert self.nnodes > 0, f'nnodes={self.nnodes:d}'
        nresults = len(ivars)
        assert nresults, ivars
        if is_points:
            _write_xyz_results_point(
                tecplot_file, self.zone_data, ivars)
        else:
            _write_xyz_results_block(
                tecplot_file, self.zone_data, ivars)

    def write_structured_zone(self, tecplot_file: TextIO, ivars: list[int],
                              log: SimpleLogger, headers_dict: dict[str, Any],
                              adjust_nids: bool=True) -> None:
        """writes a structured IxJ or IxJxK grid"""
        headers_dict = self.headers_dict
        ni = headers_dict['I']
        if 'J' in headers_dict:
            nj = headers_dict['J']
            if 'K' in headers_dict:
                nk = headers_dict['K']
                msg = f'ZONE I={ni:d}, J={nj:d}, K={nk:d} F=POINT\n'
            else:
                msg = f'ZONE I={ni:d}, J={nj:d}, F=POINT\n'
        else:
            assert 'K' not in headers_dict, list(headers_dict.keys())
            msg = f'ZONE I={ni:d}, F=POINT\n'
        tecplot_file.write(msg)
        is_points = True
        self._write_xyz_results(tecplot_file, is_points, ivars)

    def skin_elements(self) -> tuple[NDArrayN3int, NDArrayN4int]:
        """get the tris/quads from tets/hexas"""
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
            # tris = np.vstack(tris)
            # tris.sort(axis=0)
            # tris = unique_rows(tris)
        # if quads:
            # quads = np.vstack(quads)
            # quads.sort(axis=0)
            # quads = unique_rows(tris)
        return tris, quads

    def get_free_faces(self) -> list[tuple[int, int, int, int]]:
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
                if len(np.unique(face)) >= 3:
                    sort_face = tuple(sorted(face))
                    sort_face_to_element_map[sort_face].append(ie)
                    sort_face_to_face[sort_face] = face

        free_faces = []
        for sort_face, eids in sort_face_to_element_map.items():
            if len(eids) == 1:
                free_faces.append(sort_face_to_face[sort_face])
        self.log.info('finished get_free_faces')
        return free_faces

    def _slice_plane_inodes(self, inodes: list[int]) -> None:
        """TODO: doesn't remove unused nodes/renumber elements"""
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
            boolean_hexa = np.in1d(self.hexa_elements.ravel(), inodes).reshape(nhexas, 8)
            #print(boolean_hexa)
            # assert len(boolean_hexa) == self.hexa_elements.shape[0]
            assert True in boolean_hexa
            irow = np.where(boolean_hexa)[0]
            isave = np.unique(irow)
            nsave = len(isave)
            self.hexa_elements = self.hexa_elements[isave, :]
            #print(self.hexa_elements)
            #self.hexa_elements =

            # vect_lookup(self.hexa_elements) # Reassign the elements you want to change
            self.hexa_elements.reshape(nsave, 8)

        #print(boolean_hexa)
        #for hexa in hexas:
            #if
        #self.hexa_elements

    @property
    def variables(self) -> list[str]:
        return self.headers_dict['variables']

    @variables.setter
    def variables(self, variables: list[str]) -> None:
        self.headers_dict['variables'] = variables

    @property
    def xyz(self):
        if self.is_3d:
            return self.zone_data[:, :3]
        return np.zeros((0, 3), dtype='float32')

    @property
    def xy(self):
        if self.is_2d:
            return self.zone_data[:, :2]
        return np.zeros((0, 2), dtype='float32')

    @property
    def is_3d(self) -> bool:
        variables = self.headers_dict['VARIABLES']
        is_x = 'X' in variables or 'x' in variables
        is_y = 'Y' in variables or 'y' in variables
        is_z = 'Z' in variables or 'z' in variables
        is_3d = is_x and is_y and is_z
        return is_3d

    @property
    def is_2d(self) -> bool:
        variables = self.headers_dict['VARIABLES']
        is_x = 'X' in variables or 'x' in variables
        is_y = 'Y' in variables or 'y' in variables
        is_z = 'Z' in variables or 'z' in variables
        is_2d = is_x and is_y and not is_z
        return is_2d

    @property
    def nodal_results(self) -> np.ndarray:
        if self.is_3d:
            out = self.zone_data[:, 3:]
        elif self.is_2d:
            out = self.zone_data[:, 2:]
        else:
            out = self.zone_data
        return out

def _write_xyz_results_point(tecplot_file: TextIO,
                             zone_data: np.ndarray,
                             ivars: np.ndarray) -> None:
    """writes a POINT formatted result of X,Y,Z,res1,res2 output"""
    assert len(ivars), ivars
    try:
        data = zone_data[:, ivars]
    except IndexError:
        msg = 'Cant access...\n'
        msg += 'ivars = %s\n' % ivars
        #nresults = self.zone_data.shape[1]
        #msg += 'nresults = %s\n' % nresults
        msg += 'results.shape=%s\n' % str(zone_data.shape)
        raise IndexError(msg)
    fmt = ' %15.9E' * data.shape[1]
    # works in numpy 1.15.1
    np.savetxt(tecplot_file, data, fmt=fmt)

def _write_xyz_results_block(tecplot_file: TextIO,
                             zone_results: np.ndarray,
                             ivars: list[int]) -> None:
    """TODO: hasn't been tested for 2d?"""
    #print('nnodes_per_element =', nnodes_per_element)
    # for ivar in range(nnodes_per_element):

    for ivar in ivars:
        #tecplot_file.write('# ivar=%i\n' % ivar)
        vals = zone_results[:, ivar].ravel()
        msg = ''
        for ival, val in enumerate(vals):
            msg += ' %15.9E' % val
            if (ival + 1) % 5 == 0:
                tecplot_file.write(msg)
                msg = '\n'
        tecplot_file.write(msg.rstrip() + '\n')

