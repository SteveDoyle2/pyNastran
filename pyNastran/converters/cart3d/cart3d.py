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
from collections import defaultdict
from typing import Union

import numpy as np
from cpylog import SimpleLogger

from pyNastran.utils import is_binary_file
from pyNastran.converters.cart3d.cart3d_reader_writer import Cart3dReaderWriter, _write_cart3d_ascii, _write_cart3d_binary

class Cart3D(Cart3dReaderWriter):
    """Cart3d interface class"""
    model_type = 'cart3d'
    is_structured = False
    is_outward_normals = True

    def __init__(self, log=None, debug=False):
        Cart3dReaderWriter.__init__(self, log=log, debug=debug)

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

        #ax, ax0 = self._get_ax(axis, self.log)
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
        ax, ax0 = _get_ax(axis, self.log)

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

    def remove_elements(self, ielements_to_remove=None, ielements_to_keep=None) -> np.ndarray:
        if ielements_to_remove is not None and ielements_to_keep is not None:
            raise RuntimeError('Either ielements_to_remove or ielements_to_keep must not be None')
        remove = False
        keep = False
        if ielements_to_remove is not None and len(ielements_to_remove):
            remove = True
        if ielements_to_keep is not None and len(ielements_to_keep):
            keep = True

        if not(remove or keep):
            raise RuntimeError('Either ielements_to_remove or ielements_to_keep must be None')

        if remove:
            iall_elements = np.arange(self.nelements)
            ielements_to_keep = np.setdiff1d(iall_elements, ielements_to_remove)
        assert len(ielements_to_keep) > 0, ielements_to_keep
        self.elements = self.elements[ielements_to_keep, :]
        self.regions = self.regions[ielements_to_keep]
        return self.elements

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
            self._read_cart3d_binary(infilename, self._endian)
        else:
            self._read_cart3d_ascii(infilename, self._encoding, result_names)

        self.log.debug(f'npoints={self.npoints:d} nelements={self.nelements:d} nresults={self.nresults:d}')
        assert self.npoints > 0, f'npoints={self.npoints:d}'
        assert self.nelements > 0, f'nelements={self.nelements:d}'

    def write_cart3d(self, outfilename, is_binary=False, float_fmt='%6.7f'):
        """
        writes a cart3d file

        """
        assert len(self.points) > 0, 'len(self.points)=%s' % len(self.points)

        self.log.info(f'---writing cart3d...{outfilename!r}---')
        min_e = self.elements.min()
        assert min_e >= 0, 'min(elements)=%s' % min_e

        if is_binary:
            _write_cart3d_binary(outfilename, self.points, self.elements, self.regions,
                                 self.loads, self._endian)
        else:
            _write_cart3d_ascii(outfilename, self.points, self.elements, self.regions,
                                self.loads, float_fmt)
        return
        #file_fmt = 'wb' if is_binary else 'w'
        #with open(outfilename, file_fmt) as outfile:
            #int_fmt = self._write_header(outfile, self.points, self.elements, is_loads, is_binary)
            #self._write_points(outfile, self.points, is_binary, float_fmt)
            #self._write_elements(outfile, self.elements, is_binary, int_fmt)
            #self._write_regions(outfile, self.regions, is_binary)

            #if is_loads:
                #assert is_binary is False, 'is_binary=%r is not supported for loads' % is_binary
                #self._write_loads(outfile, self.loads, is_binary, float_fmt)

    def get_area(self):
        """
        Gets the element area

        Returns
        -------
        area : (n,) ndarray
            the element areas
        """

        ne = self.elements.shape[0]
        n = _get_area_vector(self.nodes, self.elements)
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
        n = _get_area_vector(self.nodes, self.elements)
        ni = np.linalg.norm(n, axis=1)
        assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)

        assert ni.min() > 0.0, ni[np.where(ni <= 0.0)[0]]
        n /= ni[:, None]  # normal vector
        return n

    def get_area_centroid_normals(self):
        """
        Gets the area, centroid, and centroidal normals

        Returns
        -------
        area : (n,) ndarray
            the element areas
        centroid : (n, 3) ndarray
            centroid of the element
        cnormals : (n, 3) ndarray
            normalized centroidal normal vectors

        """
        area, centroid, normals = get_area_centroid_normals(self.nodes, self.elements)
        return area, centroid, normals

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


def read_cart3d(cart3d_filename, log=None, debug=False, result_names=None) -> Cart3D:
    """loads a Cart3D file"""
    model = Cart3D(log=log, debug=debug)
    model.read_cart3d(cart3d_filename, result_names)
    return model


def comp2tri(in_filenames, out_filename,
             is_binary=False, float_fmt='%6.7f',
             log=None, debug=False) -> Cart3D:
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

def get_area_centroid_normals(nodes, elements):
    """
    Gets the area, centroid, and centroidal normals

    Returns
    -------
    area : (n,) ndarray
        the element areas
    centroid : (n, 3) ndarray
        centroid of the element
    cnormals : (n, 3) ndarray
        normalized centroidal normal vectors

    """
    ne = elements.shape[0]
    assert ne > 0, elements.shape
    normal = _get_area_vector(nodes, elements)
    ni = np.linalg.norm(normal, axis=1)

    area = 0.5 * ni
    assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)

    assert ni.min() > 0.0, ni[np.where(ni <= 0.0)[0]]
    normal /= ni[:, None]  # normal vector

    xyz1 = nodes[elements[:, 0], :]
    xyz2 = nodes[elements[:, 1], :]
    xyz3 = nodes[elements[:, 2], :]
    centroid = (xyz1 + xyz2 + xyz3) / 3.
    return area, centroid, normal

def _get_area_vector(nodes: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """
    Gets the area vector (unnormalized normal vector)
    Returns
    -------
    normals : (n, 3) ndarray
        unnormalized centroidal normal vectors

    """
    p1 = nodes[elements[:, 0], :]
    p2 = nodes[elements[:, 1], :]
    p3 = nodes[elements[:, 2], :]

    ne = elements.shape[0]
    avec = p2 - p1
    bvec = p3 - p1
    n = np.cross(avec, bvec)
    assert len(n) == ne, 'len(n)=%s ne=%s' % (len(n), ne)

    return n

def _get_ax(axis: Union[str, int], log: SimpleLogger) -> tuple[int, int]:
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
    log.debug("axis=%r ax=%s" % (axis, ax))

    # shift ax to the actual column index in the nodes array
    ax0 = ax if ax in [0, 1, 2] else ax - 3
    return ax, ax0

