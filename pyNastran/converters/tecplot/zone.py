from __future__ import annotations
import copy
from collections import defaultdict
from typing import TextIO, Optional, Any, Union # , TYPE_CHECKING

import numpy as np
from pyNastran.nptyping_interface import NDArrayN3int, NDArrayN4int # NDArrayN3float,

#if TYPE_CHECKING:  # pragma: no cover
from cpylog import SimpleLogger

Face = tuple[int, int, int, int]


#class CaseInsensitiveDict(dict):
    #def __contains__(self, key: str) -> bool:
        #val_in = dict.__contains__(self, key.upper())
        #return val_in
    #def __getitem__(self, key) -> Any:
        #val = dict.__getitem__(self, key.upper())
        ##log.info("GET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        #return val
    #def __setitem__(self, key, val) -> None:
        ##log.info("SET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        #dict.__setitem__(self, key.upper(), val)


class CaseInsensitiveDictAlternates(dict):
    def __init__(self, *args, **kwargs):
        self._key_map: dict[str, str] = {}
        super().__init__(*args, **kwargs)

    def set_key_map(self, key_map: dict[str, str]) -> None:
        _key_map = {}
        for key, value in key_map.items():
            key_upper = key.upper()
            value_upper = value.upper()
            _key_map[key_upper] = value_upper
        self._key_map = _key_map

    def __contains__(self, key: str) -> bool:
        key_upper = key.upper()
        val_in = dict.__contains__(self, key_upper)
        if key:
            return val_in
        return key in self._key_map

    def __getitem__(self, key: str) -> Union[int, float, str]:
        key_upper = key.upper()
        key2 = self._key_map.get(key_upper, key_upper)
        #print(f'key_upper={key_upper} key2={key2}')
        try:
            val = dict.__getitem__(self, key2)
        except KeyError:
            #print(f'key={key!r} key_upper={key_upper!r} key2={key2!r} self._key_map={self._key_map}')
            raise
        #log.info("GET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        return val

    def __setitem__(self, key: str, val: Union[int, float, str, list[str]]) -> None:
        #log.info("SET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        key_upper = key.upper()
        key2 = self._key_map.get(key_upper, key_upper)
        dict.__setitem__(self, key2, val)

class TecplotDict(CaseInsensitiveDictAlternates):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        alternate_map = {
            'data_packing': 'datapacking',
            'ZONE_TYPE': 'ZONETYPE',
        }
        self.set_key_map(alternate_map)

class Zone:
    def __init__(self, log: SimpleLogger):
        self.log = log

        self.name = '???'
        self.strand_id = -1
        self.headers_dict = TecplotDict()
        self.zone_data = np.zeros((0, 0), dtype='float32')
        self.element_data = np.zeros((0, 0), dtype='float32')
        self.tet_elements = np.zeros((0, 4), dtype='int32')
        self.hexa_elements = np.zeros((0, 8), dtype='int32')
        self.quad_elements = np.zeros((0, 4), dtype='int32')
        self.tri_elements = np.zeros((0, 3), dtype='int32')

        # TODO: does or does not consider xyz?
        #self.variables: list[str] = []

        # VARLOCATION=([3-7,10]=CELLCENTERED, [11-12]=CELLCENTERED)
        # default=NODAL
        #self.nodal_results = np.array([], dtype='float32')
        #self.centroidal_results = np.array([], dtype='float32')

    @property
    def nxyz(self) -> int:
        nxyz = 0
        variables: list[str] = self.headers_dict['variables']
        for var in variables:
            if var.lower() in {'x', 'y', 'z'}:
                nxyz += 1
        return nxyz

    @classmethod
    def set_zone_from_360(cls,
                          log: SimpleLogger,
                          header_dict: dict[str, Any],
                          variables: list[str],
                          name: str,
                          zone_type: str,
                          data_packing: str,
                          strand_id: int,
                          tris: Optional[np.ndarray]=None,
                          quads: Optional[np.ndarray]=None,
                          tets: Optional[np.ndarray]=None,
                          hexas: Optional[np.ndarray]=None,
                          zone_data: Optional[np.ndarray]=None,
                          element_data: Optional[np.ndarray]=None):
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
        nvars_node_actual = zone_data.shape[1]
        nvars_element_actual = 0 if element_data is None else element_data.shape[1]
        nvars_actual = nvars_node_actual + nvars_element_actual
        if not nvars == nvars_actual:
            raise RuntimeError(f'variables={variables} nvars={nvars}; '
                               f'nodal_results.shape={zone_data.shape} nxyz={nxyz}\n'
                               f'nvars_node_actual={nvars_node_actual:d} nvars_element_actual={nvars_element_actual:d}\n'
                               f'nvars={nvars} nvars_actual={nvars_actual}')

        #zone_data = _add_rho_uvw(zone_data, variables)
        zone_data = _add_ttot(zone_data, variables)
        zone_data = _add_ptot(zone_data, variables)
        zone_data = _add_rhotot(zone_data, variables)
        zone_data = _add_mixture_R_Cp(zone_data, variables)
        zone_data = _add_mixture_Cv(zone_data, variables)

        zone.zone_data = zone_data
        zone.element_data = element_data
        header_dict['variables'] = variables

        #zone.headers_dict = copy.copy(header_dict)
        zone.headers_dict['datapacking'] = data_packing
        zone.headers_dict['zonetype'] = zone_type
        zone.variables = variables

        #print(str(zone))
        return zone

    def variables_exist(self,
                        variables: list[str],
                        variables_to_check: list[str],
                        variables_to_save: list[str]) -> bool:
        return variables_exist(variables, variables_to_check, variables_to_save)

    @property
    def is_unstructured(self) -> bool:
        """are there unstructured elements"""
        is_unstructured: bool = (
            len(self.tri_elements) > 0 or
            len(self.quad_elements) > 0 or
            len(self.tet_elements) > 0 or
            len(self.hexa_elements) > 0
        )
        return is_unstructured

    @property
    def is_structured(self) -> bool:
        """is the model in plot3d style format"""
        is_structured = not self.is_unstructured
        return is_structured

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
        title = self.headers_dict['TITLE']
        assert isinstance(title, str), title
        return title
    @title.setter
    def title(self, title: str) -> None:
        self.headers_dict['TITLE'] = title

    @property
    def is_point(self) -> bool:
        if 'F' in self.headers_dict:
            raise RuntimeError(self.headers_dict)
        if 'DATAPACKING' not in self.headers_dict:
            return True
        is_point = (self.headers_dict['DATA_PACKING'] == 'POINT')
        is_block = (self.headers_dict['DATA_PACKING'] == 'BLOCK')
        assert is_point or is_block
        return is_point

    @property
    def is_block(self) -> bool:
        return not self.is_point

    @property
    def datapacking(self) -> Optional[str]:
        """'POINT', 'BLOCK'"""
        try:
            data_packing = self.headers_dict['DATAPACKING']
            assert isinstance(data_packing, str), data_packing
            return data_packing
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
        else:  # pragma: no cover
            raise RuntimeError(zonetype)
        return zone_type_int

    def quads_to_tris(self) -> None:
        zonetype = self.headers_dict['ZONETYPE']
        if zonetype == 'FEQUADRILATERAL':
            n3 = self.quad_elements[:, 2]
            n4 = self.quad_elements[:, 3]
            if np.all(n3 == n4):
                self.tri_elements = self.quad_elements[:, :3]
                assert self.tri_elements.shape[1] == 3, self.tri_elements.shape
                self.quad_elements = np.zeros((0, 4), self.quad_elements.dtype)
                self.headers_dict['ZONETYPE'] = 'FETRIANGLE'

        #if zonetype == 'FETRIANGLE',  'FETETRAHEDRON', 'FEBRICK'


    def split_elements(self, ntri_nodes: int=1) -> None:
        r"""
        Splits elements and linearly interpolates the data.

        Supports:
         - 2d/3d
         - unused nodes

        Falls apart if:
         - ???

            n=0; N=2         n=1;N=3;Ne=4  n=2;N=4;Ne=9=1+3+5
              *                *3              *        1
             / \              / \             / \     1
            /   \            / C \           *---*      2
           /     \    ->   6*-----*5        / \ / \   3
          /       \        / \ D / \       *---*---*    3
         /         \      / A \ / B \     / \ / \ / \ 5
        *-----------*   1*-----*-----*2  *---*---*---*  4
                               4

        n  = [0, 1, 2,  3,  4,  5,  6,  7, 8,    9]  # number of intermediate nodes
        N  = [2, 3, 4,  5,  6,  7,  8,  9, 10,  11]  # total number of nodes
        Ne = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]  # number of elements
        Nn = [1, 3, 6, 10, 15, 21, 28, 36, 45,  55]  # number of nodes (ideally)

        N = n + 2
        Ne = (n+1)^2
        n = Nn / (n + 1) * 2 - 2
        Nn = (n + 1)(n + 2) // 2

        Refinements
        -----------
        N =  [2, 3, 5,   9,   33]
        n =  [0, 1, 3,   7,   31]
        Ne = [1, 4, 16, 64, 1024]


                        n=2                           n=5; N=7; Ne=36
                          *                             *
                         / \                           / \           1
                        /   \                         *---*
                       /     \                       / \ / \         3
                      *-------*                     *---*---*
                     / \     / \                   / \ / \ / \       5
                    /   \   /   \                 *---*---*---*
                   /     \ /     \               / \ / \ / \ / \     7
                  *-------*-------*             *---*---*---*---*
                 / \     / \     / \           / \ / \ / \ / \ / \   9
                /   \   /   \   /   \         *---*---*---*---*---*
               /     \ /     \ /     \       / \ / \ / \ / \ / \ / \ 11
              *-------*-------*-------*     *---*---*---*---*---*---*

        Using the n=2 (or more) method would create fewer duplicate nodes
        N=1 is something like 4*n^2
        nnodes0=4 nelements0=1; dt=0.00
        nnodes1=7 nelements1=4; dt=0.00
        nnodes2=19 nelements2=16; dt=0.00
        nnodes3=67 nelements3=64; dt=0.00
        nnodes4=259 nelements4=256; dt=0.00
        nnodes5=1027 nelements5=1024; dt=0.00
        nnodes6=4099 nelements6=4096; dt=0.00
        nnodes7=16387 nelements7=16384; dt=0.00
        nnodes8=65539 nelements8=65536; dt=0.02
        nnodes9=262147 nelements9=262144; dt=0.05
        nnodes10=1048579 nelements10=1048576; dt=0.12
        nnodes11=4194307 nelements11=4194304; dt=0.57
        nnodes12=16777219 nelements12=16777216; dt=2.16
        """
        zonetype = self.headers_dict['ZONETYPE']
        if zonetype == 'FETRIANGLE':
            assert ntri_nodes > 0, ntri_nodes
            #ntri_nodes0 = ntri_nodes  # for debugging

            # check for no unused nodes
            #nnodes0 = self.zone_data.shape[0]
            #unodes = np.unique(self.tri_elements.ravel())
            #assert nnodes0 == len(unodes), 'run nodal_equivalence on the zone first'

            import time
            t0 = time.time()
            # not super fast, but it's a lot more general
            for irefine in range(ntri_nodes):
                t1 = time.time()
                nnodes0 = self.zone_data.shape[0]
                nelements0 = self.tri_elements.shape[0]

                assert self.tri_elements.min() == 0, self.tri_elements.min()
                n1 = self.tri_elements[:, 0]
                n2 = self.tri_elements[:, 1]
                n3 = self.tri_elements[:, 2]

                # we're lying a bit with the xyz1 name and
                # working with all the result data
                xyz1 = self.zone_data[n1, :]
                xyz2 = self.zone_data[n2, :]
                xyz3 = self.zone_data[n3, :]
                #xyz = self.zone_data[:, :3]

                # create nodes 456
                #
                # to simplify creating matched nodes, we'll
                # create all the 4s, then 5s, and then 6s
                #
                # we're going to create a xyz4 that matches xyz5 or xyz6
                # for the neighboring element to make stacking easier
                # we waste some memory, but it's fine
                assert self.tri_elements.min() == 0
                n4 = np.arange(nnodes0, nnodes0 + nelements0)
                n5 = n4 + nelements0
                n6 = n5 + nelements0
                xyz4 = (xyz1 + xyz2) / 2
                xyz5 = (xyz2 + xyz3) / 2
                xyz6 = (xyz3 + xyz1) / 2
                self.zone_data = np.vstack([self.zone_data, xyz4, xyz5, xyz6])

                # follow the pattern up above
                tria = np.column_stack([n1, n4, n6])
                trib = np.column_stack([n4, n2, n5])
                tric = np.column_stack([n6, n5, n3])
                trid = np.column_stack([n4, n5, n6])
                self.tri_elements = np.vstack([tria, trib, tric, trid])
                dt = time.time() - t1
                self.log.debug(f'nnodes{irefine}={nnodes0} nelements{irefine}={nelements0}; dt={dt:.2f}')

            dt = time.time() - t0
            nnodes_final = self.zone_data.shape[0]
            nelements_final = self.tri_elements.shape[0]
            self.log.debug(f'nnodes_final={nnodes_final} nelements_final={nelements_final} dt={dt:.0f}')
            self.log.debug('split triangles')
            # TODO: this would be the place to do a nodal equivalence

        elif zonetype in {'FEQUADRILATERAL', 'FETETRAHEDRON', 'FEBRICK'}:
            # not supported
            self.log.warning(f'splitting {zonetype} is not supported')
            pass
        else:  # pragma: no cover
            raise RuntimeError(zonetype)

    def demote_elements(self):
        """
        Tecplot is wasting memory, which is no big deal,
        but mistyping elements makes some things more difficult
         - accurate centroids
         - splitting elements (and not creating mixed quad/tri meshes)
        """
        zonetype = self.headers_dict['ZONETYPE']
        if zonetype in {'FETRIANGLE', 'FETETRAHEDRON'}:
            # no demotion possible
            pass
        elif zonetype == 'FEQUADRILATERAL':
            quads = self.quad_elements
            #n1 = quads[:, 0]
            #n2 = quads[:, 0]
            n3 = quads[:, 0]
            n4 = quads[:, 0]
            itri = (n3 == n4)
            demote_to_tris = itri.all()
            if demote_to_tris:
                zonetype = self.headers_dict['ZONETYPE'] ='FETRIANGLE'
                self.tri_elements = self.quad_elements[:, :3]
                assert self.tri_elements.shape[1] == 3, self.tri_elements.shape
                self.quad_elements = np.zeros((0, 4), dtype='int32')
        elif zonetype == 'FEBRICK':
            self.log.warning('demoting FEBRICK is not supported')
        else:  # pragma: no cover
            raise RuntimeError(zonetype)

    @property
    def zonetype(self) -> str:
        """FEBrick,  FETETRAHEDRON,  FETRIANGLE,  FEQUADRILATERAL"""
        #headers_dict['ZONETYPE'].upper() # FEBrick
        #try:
        zonetype = self.headers_dict['ZONETYPE']
        assert isinstance(zonetype, str), zonetype
        return zonetype.upper()
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

        unused_xy_shape = str(self.xy.shape)
        xyz_shape = str(self.xyz.shape)
        is3d = self.is_3d
        if 'I' in self.headers_dict:
            msgi = repr_nijk(self.headers_dict)
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

    def write_unstructured_zone(self, tecplot_file: TextIO, ivars: list[int], is_points: bool,
                                nnodes: int, nelements: int, zone_type: int, log: SimpleLogger,
                                is_tris: bool, is_quads: bool, is_tets: bool, is_hexas: bool,
                                adjust_nids: bool=True) -> None:
        msg = 'ZONE '
        self.log.debug(f'is_points = {is_points}')
        datapacking = 'POINT' if  is_points else 'BLOCK'
        nvar_node = self.zone_data.shape[1]
        nvar_element = 0 if self.element_data is None else self.element_data.shape[1]

        # varlocation=([1,2]=nodal,[3]=cellcentered)
        nodaL_var_list = list(range(1, nvar_node+1))
        element_var_list = list(range(nvar_node+1, nvar_node+nvar_element+1))
        if nvar_node and nvar_element:
            var_location = f'VARLOCATION=({nodaL_var_list}=NODAL,{element_var_list}=CELLCENTERED)'
        elif nvar_node:
            var_location = ''
        elif nvar_element:
            var_location = f'VARLOCATION=({element_var_list}=CELLCENTERED)'
        else:  # pramga: no cover
            raise RuntimeError((nvar_node, nvar_element))

        msg += (
                f' T=\"{self.title}\", n={nnodes:d}, e={nelements:d}, '
            f'ZONETYPE={zone_type}, DATAPACKING={datapacking}{var_location}\n'
        )
        tecplot_file.write(msg)

        ivars_array = np.asarray(ivars, dtype='int32')
        self._write_xyz_results(tecplot_file, is_points, ivars_array)
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
                           ivars: np.ndarray) -> None:
        """writes XY/XYZs and results in POINT or BLOCK format"""
        assert self.nnodes > 0, f'nnodes={self.nnodes:d}'
        nresults = len(ivars)
        assert nresults, ivars
        if is_points:
            _write_xyz_results_point(
                tecplot_file, self.zone_data, ivars)
        else:
            _write_xyz_results_block(
                tecplot_file, self.zone_data, self.element_data, ivars)

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
                msg = f'ZONE T=\"{self.title}\", I={ni:d}, J={nj:d}, K={nk:d} F=POINT\n'
            else:
                msg = f'ZONE T=\"{self.title}\", I={ni:d}, J={nj:d}, F=POINT\n'
        else:
            assert 'K' not in headers_dict, list(headers_dict.keys())
            msg = f'ZONE T=\"{self.title}\", I={ni:d}, F=POINT\n'
        tecplot_file.write(msg)
        is_points = True

        ivars_array = np.asarray(ivars, dtype='int32')
        self._write_xyz_results(tecplot_file, is_points, ivars_array)

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
            for face in (btm, top, left, right, front, back):
                if len(np.unique(face)) >= 3:
                    sort_face = tuple(sorted(face))
                    sort_face_to_element_map[sort_face].append(ie)
                    sort_face_to_face[sort_face] = face

        free_faces: list[Face] = []
        for sort_face, eids in sort_face_to_element_map.items():
            if len(eids) == 1:
                facei: Face = sort_face_to_face[sort_face]
                free_faces.append(facei)
        self.log.info('finished get_free_faces')
        return free_faces

    def _slice_plane_inodes(self, inodes: np.ndarray | list[int]) -> None:
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
            boolean_hexa = np.isin(self.hexa_elements.ravel(), inodes).reshape(nhexas, 8)
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
        variables = self.headers_dict['variables']
        assert isinstance(variables, list), variables
        return variables

    @variables.setter
    def variables(self, variables: list[str]) -> None:
        self.headers_dict['variables'] = variables

    @property
    def xyz(self) -> np.ndarray:
        if self.is_3d:
            return self.zone_data[:, :3]
        return np.zeros((0, 3), dtype='float32')

    @property
    def xy(self) -> np.ndarray:
        if self.is_2d:
            return self.zone_data[:, :2]
        return np.zeros((0, 2), dtype='float32')

    @property
    def is_3d(self) -> bool:
        variables: list[str] = self.variables
        is_x = 'X' in variables or 'x' in variables
        is_y = 'Y' in variables or 'y' in variables
        is_z = 'Z' in variables or 'z' in variables
        is_3d = is_x and is_y and is_z
        return is_3d

    @property
    def is_2d(self) -> bool:
        variables: list[str] = self.variables
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


def repr_nijk(headers_dict: TecplotDict) -> str:
    ni = headers_dict['I']
    if 'J' in headers_dict:
        nj = headers_dict['J']
        if 'K' in headers_dict:
            nk = headers_dict['K']
            msgi = f'  nI={ni} nJ={nj} nK={nk}\n'
        else:
            msgi = f'  nI={ni} nJ={nj}\n'
    else:
        assert 'K' not in headers_dict, list(headers_dict.keys())
        msgi = f'  nI={ni}\n'
    return msgi

def variables_exist(variables: list[str],
                    variables_to_check: list[str],
                    variables_to_save: list[str]) -> bool:
    """
    Checks to see if the required variable exist.  If so, and the
    new variable doesn't exist, then make it.
    """
    all_exist = all(var in variables for var in variables_to_check)
    if all(var in variables for var in variables_to_save):
        all_exist = False
    return all_exist

def _add_ptot(zone_data: np.ndarray,
              variables: list[str]) -> np.ndarray:
    """create total pressure"""
    if not variables_exist(variables, ['mach', 'tt', 'gamma'], ['ptot']):
        return zone_data

    ip = variables.index('p')  # pressure
    imach = variables.index('mach')
    igamma = variables.index('gamma')

    p = zone_data[:, ip]
    gamma = zone_data[:, igamma]
    mach = zone_data[:, imach]
    nnodes = len(mach)

    # https://en.wikipedia.org/wiki/Total_pressure
    gm1 = gamma - 1
    ktemp = 1 + gm1 / 2 * mach ** 2
    exp = gamma / gm1
    ptot = (p * ktemp ** exp).reshape((nnodes, 1))

    zone_data = np.hstack([zone_data, ptot])
    variables.append('ptot')
    return zone_data


def _add_ttot(zone_data: np.ndarray,
              variables: list[str]) -> np.ndarray:
    """create total temperature"""
    if not variables_exist(variables, ['mach', 'tt', 'gamma'], ['ttot']):
        return zone_data

    itt = variables.index('tt')  # translational-rotational temperature
    imach = variables.index('mach')
    igamma = variables.index('gamma')

    tt = zone_data[:, itt]
    gamma = zone_data[:, igamma]
    mach = zone_data[:, imach]
    nnodes = len(mach)

    # t0 / t = 1 + (gamma - 1) / 2 * mach ^ 2
    gm1 = gamma - 1
    ktemp = 1 + gm1 / 2 * mach ** 2
    ttot = (tt * ktemp).reshape((nnodes, 1))

    zone_data = np.hstack([zone_data, ttot])
    variables.append('ttot')
    return zone_data

def _add_rhotot(zone_data: np.ndarray,
                variables: list[str]) -> np.ndarray:
    """
    Create total density
    https://www.hkdivedi.com/2019/01/stagnation-properties-in-fluid.html#:~:text=Stagnation%20density-,The%20point%2C%20at%20which%20resultant%20velocity%20of%20the%20fluid%20becomes,be%20called%20as%20stagnation%20density.
    """
    if not variables_exist(variables, ['mach', 'rho', 'gamma'], ['rhotot']):
        return zone_data
    irho = variables.index('rho')
    imach = variables.index('mach')
    igamma = variables.index('gamma')

    rho = zone_data[:, irho]
    gamma = zone_data[:, igamma]
    mach = zone_data[:, imach]
    nnodes = len(mach)

    # t0 / t = 1 + (gamma - 1) / 2 * mach ^ 2
    gm1 = gamma - 1
    ktemp = 1 + gm1 / 2 * mach ** 2
    exp = 1 / gm1
    ttot = (rho * ktemp ** exp).reshape((nnodes, 1))

    zone_data = np.hstack([zone_data, ttot])
    variables.append('rhotot')
    return zone_data

def _add_rho_uvw(zone_data: np.ndarray,
                 variables: list[str]) -> np.ndarray:
    if not variables_exist(variables, ['rho', 'u', 'v', 'w'], ['rhoUVW']):
        return zone_data
    irho = variables.index('rho')
    iu = variables.index('u')
    iv = variables.index('v')
    iw = variables.index('w')

    rho = zone_data[:, irho]
    uvw = zone_data[:, [iu, iv, iw]]
    nnodes = len(rho)

    uvw_mag = np.linalg.norm(uvw, axis=1)
    rho_uvw_mag = (rho * uvw_mag).reshape((nnodes, 1))
    rho_uvw = np.abs(rho[:, np.newaxis] * uvw)
    zone_data = np.hstack([zone_data, rho_uvw, rho_uvw_mag])
    variables.extend(['rhoU', 'rhoV', 'rhoW', 'rhoUVW'])
    return zone_data

def _add_mixture_R_Cp(zone_data: np.ndarray,
                      variables: list[str]) -> np.ndarray:
    """https://www.tau.ac.il/~tsirel/dump/Static/knowino.org/wiki/Specific_heat_ratio.html."""
    if not variables_exist(variables, ['mixture_mol_weight', 'gamma'], ['mix_Cp']):
        return zone_data

    imolecular_weight = variables.index('mixture_mol_weight')
    igamma = variables.index('gamma')

    molecular_weight = zone_data[:, imolecular_weight] # g/mol
    Runv = 8314 # J/(kg*K)
    R: np.ndarray = Runv / molecular_weight  # J/(kg*K)

    gamma = zone_data[:, igamma]
    nnodes = len(gamma)

    # Cp = gamma*R/(gamma-1)
    # Cv = R / (gamma-1)
    #
    Cp = gamma * R / (gamma-1)
    zone_data = np.hstack([
        zone_data,
        R.reshape((nnodes, 1)),
        Cp.reshape((nnodes, 1)),
    ])
    variables.extend(['R', 'mixture_Cp'])
    return zone_data

def _add_mixture_Cv(zone_data: np.ndarray,
                    variables: list[str]) -> np.ndarray:
    """https://www.tau.ac.il/~tsirel/dump/Static/knowino.org/wiki/Specific_heat_ratio.html."""
    if not variables_exist(variables, ['R', 'gamma'], ['mix_Cp']):
        return zone_data

    iR = variables.index('R')  # gas constant
    igamma = variables.index('gamma')

    R = zone_data[:, iR]
    gamma = zone_data[:, igamma]
    nnodes = len(gamma)

    # Cp = gamma*R/(gamma-1)
    # Cv = R / (gamma-1)
    #
    Cv = R / (gamma-1)
    zone_data = np.hstack([zone_data, Cv.reshape((nnodes, 1))])
    variables.append('mixture_Cv')
    return zone_data

def _write_xyz_results_point(tecplot_file: TextIO,
                             zone_data: np.ndarray,
                             ivars: np.ndarray | list[int]) -> None:
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
                             element_results: np.ndarray,
                             ivars: np.ndarray) -> None:
    """TODO: hasn't been tested for 2d?"""
    #print('nnodes_per_element =', nnodes_per_element)
    # for ivar in range(nnodes_per_element):

    nvar_node = zone_results.shape[1]
    nvar_element = element_results.shape[1] if element_results is not None else 0
    for ivar in ivars:
        #tecplot_file.write('# ivar=%i\n' % ivar)
        if ivar < nvar_node:
            vals = zone_results[:, ivar].ravel()
        else:
            vals = element_results[:, ivar-nvar_node].ravel()
        msg = ''
        for ival, val in enumerate(vals):
            msg += ' %15.9E' % val
            if (ival + 1) % 5 == 0:
                tecplot_file.write(msg)
                msg = '\n'
        tecplot_file.write(msg.rstrip() + '\n')
