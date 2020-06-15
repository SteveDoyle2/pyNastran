from __future__ import annotations
from collections import defaultdict
from typing import Tuple, List, Dict, TextIO, Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.nptyping import NDArrayN3float, NDArrayN3int, NDArrayN4int

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

class CaseInsensitiveDict(dict):
    def __getitem__(self, key):
        val = dict.__getitem__(self, key.upper())
        #log.info("GET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        return val
    def __setitem__(self, key, val):
        #log.info("SET %s['%s'] = %s" % str(dict.get(self, 'name_label')), str(key), str(val)))
        dict.__setitem__(self, key.upper(), val)


class Zone:
    def __init__(self, log: SimpleLogger):
        self.log = log

        self.headers_dict = CaseInsensitiveDict()
        self.xy = np.array([], dtype='float32')
        self.xyz = np.array([], dtype='float32')
        self.tet_elements = np.array([], dtype='int32')
        self.hexa_elements = np.array([], dtype='int32')
        self.quad_elements = np.array([], dtype='int32')
        self.tri_elements = np.array([], dtype='int32')
        self.variables = []

        # VARLOCATION=([3-7,10]=CELLCENTERED, [11-12]=CELLCENTERED)
        # default=NODAL
        self.nodal_results = np.array([], dtype='float32')
        self.A = None
        #self.centroidal_results = np.array([], dtype='float32')

    @property
    def nnodes(self) -> int:
        """gets the number of nodes in the model"""
        return self.xyz.shape[0] + self.xy.shape[0]

    @property
    def nelements(self) -> int:
        """gets the number of elements in the model"""
        return (self.hexa_elements.shape[0] + self.tet_elements.shape[0] +
                self.quad_elements.shape[0] + self.tri_elements.shape[0])

    @property
    def result_names(self) -> List[str]:
        """gets the variables"""
        return self.variables

    @result_names.setter
    def result_names(self, vals: List[str]) -> None:
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
    def datapacking(self) -> Optional[str]:
        """'POINT', 'BLOCK'"""
        try:
            return self.headers_dict['DATAPACKING']
        except KeyError:
            return None
    @property
    def zonetype(self) -> Optional[str]:
        """FEBrick,  FETETRAHEDRON,  FETRIANGLE,  FEQUADRILATERAL"""
        try:
            return self.headers_dict['ZONETYPE']
        except KeyError:
            return None

    @property
    def results(self) -> np.ndarray:
        self.log.warning('depecrecated tecplot.results for tecplot.nodal_results')
        return self.nodal_results
    @results.setter
    def results(self, nodal_results: np.ndarray) -> None:
        self.log.warning('depecrecated tecplot.results for tecplot.nodal_results')
        self.nodal_results = nodal_results

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
        xy_shape = str(self.xy.shape)
        xyz_shape = str(self.xyz.shape)
        a_shape = str(self.A.shape) if self.A is not None else None
        is3d = is_3d(self.headers_dict)
        if 'I' in self.headers_dict:
            msgi = self.repr_nijk()
        else:
            msgi = (
                '  datapacking = %r\n'
                '  zonetype = %r\n' % (self.datapacking, self.zonetype)
            )


        title2 = '  T = %r\n' % self.headers_dict['T'] if 'T' in self.headers_dict else ''
        msg = (
            'Zone:'
            #'  filename = %r\n'
            #'  is_mesh = %s\n'
            '  variables = %s\n'
            '  is3d = %s\n'
            '    xy.shape = %s\n'
            '    xyz.shape = %s\n'
            '  A.shape = %s\n'
            '  nodal_results.shape = %s\n'
            #'  headers_dict = %s\n'
            '  title = %r\n%s%s'
            #'  use_cols = %s\n'
            #'  dtype = %s\n'
            % (
                #self.tecplot_filename, # self.is_mesh,
                self.variables,
                is3d, xy_shape, xyz_shape,
                a_shape, str(self.nodal_results.shape),
                #str(self.headers_dict),
                self.title, title2, msgi,
                #self.use_cols,
                #self.dtype,
            )
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

    def write_unstructured_zone(self, tecplot_file: TextIO, ivars: List[int], is_points: bool,
                                nnodes: int, nelements: int, zone_type: int, log: SimpleLogger,
                                is_tris: bool, is_quads: bool, is_tets: bool, is_hexas: bool,
                                adjust_nids: bool=True) -> None:
        msg = 'ZONE '
        self.log.info('is_points = %s' % is_points)
        datapacking = 'POINT' if  is_points else 'BLOCK'
        msg += f' n={nnodes:d}, e={nelements:d}, ZONETYPE={zone_type}, DATAPACKING={datapacking}\n'
        tecplot_file.write(msg)

        self._write_xyz_results(tecplot_file, is_points, ivars)
        self._write_elements(tecplot_file, nnodes,
                             is_tris, is_quads, is_tets, is_hexas,
                             adjust_nids=adjust_nids)

    def determine_element_type(self) -> Tuple[bool, bool, bool, Optional[str],
                                              bool, bool, bool, bool]:
        """
        Returns
        -------
        is_structured : bool
           is this a structured grid
        is_unstructured : bool
           is this an unstructured grid
        is_points : bool
           ???
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

        is_structured = False
        is_unstructured = False
        zone_type = None
        for etype, elements in etype_elements:
            if not len(elements):
                continue
            if etype == 'CHEXA' and len(elements):
                # is_points = False
                is_hexas = True
                #nnodes_per_element = 8
                zone_type = 'FEBrick'
            elif etype == 'CTETRA' and len(elements):
                # is_points = False
                is_tets = True
                #nnodes_per_element = 4
                zone_type = 'FETETRAHEDRON'
            elif etype == 'CTRIA3' and len(elements):
                # is_points = True
                is_tris = True
                #nnodes_per_element = 3
                zone_type = 'FETRIANGLE'
            elif etype == 'CQUAD4' and len(elements):
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
        return is_structured, is_unstructured, is_points, zone_type, is_tris, is_quads, is_tets, is_hexas

    def _write_elements(self, tecplot_file: TextIO, nnodes: int,
                        is_tris: bool, is_quads: bool, is_tets: bool, is_hexas: bool,
                        adjust_nids: bool=True) -> None:
        """Writes the unstructured elements.  Verifies that nodes are sequential."""
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
            self.log.error(msg)
            raise RuntimeError(msg)
        # assert elements.min() == 1, elements.min()
        # assert elements.max() == nnodes, elements.max()

        if adjust_nids:
            elements += 1

        for element in elements:
            tecplot_file.write(efmt % tuple(element))

    def _write_xyz_results(self, tecplot_file: TextIO, is_points: bool,
                           ivars: List[int]) -> None:
        """writes XY/XYZs and results in POINT or BLOCK format"""
        # xyz
        xyz = self.xyz
        xy = self.xy

        nnodes3d = self.xyz.shape[0]
        nnodes2d = self.xy.shape[0]
        if nnodes2d and nnodes3d:
            raise RuntimeError('2d and 3d nodes is not supported')
        elif nnodes2d:
            nodes = xy
            word = 'xy'
            ndim = 2
        elif nnodes3d:
            nodes = xyz
            word = 'xyz'
            ndim = 3
        else:  # pragma: no cover
            raise RuntimeError('failed to find 2d/3d nodes')

        #print('nnodes', nodes.shape)
        #print('results', self.nodal_results.shape)
        assert self.nnodes > 0, 'nnodes=%s' % self.nnodes
        nresults = len(ivars)
        if is_points:
            _write_xyz_results_point(tecplot_file, nodes, self.nodal_results, nresults, ivars,
                                     ndim=ndim, word=word)

        else:
            _write_xyz_results_block(tecplot_file, nodes, self.nodal_results, nresults, ivars,
                                     ndim=ndim, word=word)

    def write_structured_zone(self, tecplot_file: TextIO, ivars: List[int],
                              log: SimpleLogger, headers_dict: Dict[str, Any],
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

    def skin_elements(self) -> Tuple[NDArrayN3int, NDArrayN4int]:
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

    def get_free_faces(self) -> List[Tuple[int, int, int, int]]:
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

    def _slice_plane_inodes(self, inodes: List[int]) -> None:
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



def is_3d(headers_dict: Dict[str, Any]) -> bool:
    #print(headers_dict)
    variables = headers_dict['VARIABLES']
    is_3d = 'Z' in variables # or 'z' in variables
    return is_3d

def _write_xyz_results_point(tecplot_file, nodes, nodal_results, nresults, ivars, ndim=3, word='xyz'):
    """writes a POINT formatted result of X,Y,Z,res1,res2 output"""
    if nresults:
        try:
            res = nodal_results[:, ivars]
        except IndexError:
            msg = 'Cant access...\n'
            msg += 'ivars = %s\n' % ivars
            msg += 'nresults = %s\n' % nresults
            msg += 'results.shape=%s\n' % str(nodal_results.shape)
            raise IndexError(msg)

        try:
            data = np.hstack([nodes, res])
        except ValueError:
            msg = 'Cant hstack...\n'
            msg += '%s.shape=%s\n' % (word, str(nodes.shape))
            msg += 'results.shape=%s\n' % str(nodal_results.shape)
            raise ValueError(msg)
        fmt = ' %15.9E' * (ndim + nresults)
    else:
        data = nodes
        fmt = ' %15.9E %15.9E %15.9E'

    # works in numpy 1.15.1
    np.savetxt(tecplot_file, data, fmt=fmt)

def _write_xyz_results_block(tecplot_file: TextIO, nodes: NDArrayN3float,
                             nodal_results, nresults, ivars: List[int],
                             ndim: int=3, word: str='xyz') -> None:
    """TODO: hasn't been tested for 2d?"""
    #nvalues_per_line = 5
    for ivar in range(ndim):
        #tecplot_file.write('# ivar=%i\n' % ivar)
        vals = nodes[:, ivar].ravel()
        msg = ''
        for ival, val in enumerate(vals):
            msg += ' %15.9E' % val
            if (ival + 1) % 3 == 0:
                tecplot_file.write(msg)
                msg = '\n'
        tecplot_file.write(msg.rstrip() + '\n')

    if nresults:
        #print('nnodes_per_element =', nnodes_per_element)
        # for ivar in range(nnodes_per_element):
        for ivar in ivars:
            #tecplot_file.write('# ivar=%i\n' % ivar)
            vals = nodal_results[:, ivar].ravel()
            msg = ''
            for ival, val in enumerate(vals):
                msg += ' %15.9E' % val
                if (ival + 1) % 5 == 0:
                    tecplot_file.write(msg)
                    msg = '\n'
            tecplot_file.write(msg.rstrip() + '\n')
