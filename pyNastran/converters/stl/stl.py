#pylint:  disable=C0111
from __future__ import print_function
import copy
from struct import unpack, Struct, pack

from six import iteritems
from six.moves import zip, range

import numpy as np
import scipy


from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils import is_binary_file
from pyNastran.utils.log import get_logger

def read_stl(stl_filename, log=None, debug=False):
    """

    Reads an STL file

    Parameters
    ----------
    stl_filename : str
        the filename to read
    """
    stl = STL(log=log, debug=debug)
    stl.read_stl(stl_filename)
    return stl


class STL(object):
    #model_type = 'stl'
    #isStructured = False
    #isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')

        self.nodes = None
        self.elements = None
        self.header = ''
        self.infilename = None

    def write_stl(self, stl_out_filename, is_binary=False, float_fmt='%6.12f',
                  normalize_normal_vectors=False, stop_on_failure=True):
        """
        Writes an STL file

        Parameters
        ----------
        stl_out_filename : str
            the filename to write
        is_binary : bool; default=False
            should a binary file be written
        float_fmt : str; default='%6.12f'
            the format to use if an ASCII file is used
        normalize_normal_vectors : bool; default=False
            should the vectors be normalized
        """
        self.log.info("---writing STL file...%r---" % stl_out_filename)
        assert len(self.nodes) > 0
        solid_name = 'dummy_name'
        if is_binary:
            self.write_binary_stl(stl_out_filename,
                                  normalize_normal_vectors=normalize_normal_vectors,
                                  stop_on_failure=stop_on_failure)
        else:
            self.write_stl_ascii(stl_out_filename, solid_name, float_fmt=float_fmt,
                                 normalize_normal_vectors=normalize_normal_vectors,
                                 stop_on_failure=stop_on_failure)

    def read_stl(self, stl_filename):
        """
        Reads an STL file

        Parameters
        ----------
        stl_filename : str
            the filename to read
        """
        self.infilename = stl_filename
        self.log.info("---starting reading STL file...%r---" % self.infilename)

        if is_binary_file(stl_filename):
            self.read_binary_stl(stl_filename)
        else:
            self.read_ascii_stl(stl_filename)

        #self.log.info("nodes=%s  nelements=%s" % (self.nodes, self.nelements))
        self.log.info("---finished reading STL file...%r---" % self.infilename)
        #assert self.nodes > 0, 'nodes=%s' % self.nodes
        #assert self.nelements > 0, 'nelements=%s' % self.nelements


    def write_binary_stl(self, stl_filename, normalize_normal_vectors=False,
                         stop_on_failure=True):
        """
        Write an STL binary file

        Parameters
        ----------
        stl_filename : str
            the filename to read
        normalize_normal_vectors : bool; default=False
            should the vectors be normalized
        """
        with open(stl_filename, 'wb') as infile:
            if hasattr(self, 'header'):
                self.header.ljust(80, '\0')
                header = '%-80s' % self.header[:80]
            else:
                header = '%-80s' % stl_filename[:80]
            infile.write(pack(b'80s', header.encode('ascii')))
            #avector = [0., 0., 0.]
            #bvector = [0., 0., 0.]
            #cvector = [0., 0., 0.]
            nelements = self.elements.shape[0]
            infile.write(pack('i', nelements))
            elements = self.elements

            p1 = self.nodes[elements[:, 0], :]
            p2 = self.nodes[elements[:, 1], :]
            p3 = self.nodes[elements[:, 2], :]
            avector = p2 - p1
            bvector = p3 - p1
            n = np.cross(avector, bvector)
            del avector, bvector
            if normalize_normal_vectors:
                n /= np.linalg.norm(n, axis=1)[:, np.newaxis]

            s = Struct('12fH')
            for eid, element in enumerate(elements):
                data = s.pack(n[eid, 0], n[eid, 1], n[eid, 2],
                              p1[eid, 0], p1[eid, 1], p1[eid, 2],
                              p2[eid, 0], p2[eid, 1], p2[eid, 2],
                              p3[eid, 0], p3[eid, 1], p3[eid, 2], 0)
                infile.write(data)

    def read_binary_stl(self, stl_filename):
        """
        Read an STL binary file

        Parameters
        ----------
        stl_filename : str
            the filename to read
        """
        with open(stl_filename, 'rb') as infile:
            data = infile.read()

        self.header = data[:80]
        nelements, = unpack('i', data[80:84])
        j = 84

        inode = 0
        nodes_dict = {}
        assert nelements > 0, 'nelements=%s' % nelements
        elements = np.zeros((nelements, 3), 'int32')

        s = Struct('12fH')
        for ielement in range(nelements):
            (nx, ny, nz, ax, ay, az, bx, by, bz,
             cx, cy, cz, i) = s.unpack(data[j:j+50])

            t1 = (ax, ay, az)
            t2 = (bx, by, bz)
            t3 = (cx, cy, cz)
            if t1 in nodes_dict:
                i1 = nodes_dict[t1]
            else:
                i1 = inode
                nodes_dict[t1] = inode
                inode += 1

            if t2 in nodes_dict:
                i2 = nodes_dict[t2]
            else:
                i2 = inode
                nodes_dict[t2] = inode
                inode += 1

            if t3 in nodes_dict:
                i3 = nodes_dict[t3]
            else:
                i3 = inode
                nodes_dict[t3] = inode
                inode += 1
            elements[ielement] = [i1, i2, i3]
            j += 50
        assert inode > 0, inode
        nnodes = inode + 1 # accounting for indexing
        self.elements = elements

        nodes = np.zeros((nnodes, 3), 'float64')
        for node, inode in iteritems(nodes_dict):
            nodes[inode] = node
        self.nodes = nodes


    def _get_normals_data(self, elements, nodes=None):
        """
        This is intended as a submethod to help handle the problem of bad normals
        """
        if nodes is None:
            nodes = self.nodes
        self.log.debug("get_normals...elements.shape %s" % str(elements.shape))
        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        v12 = p2 - p1
        v13 = p3 - p1
        v123 = np.cross(v12, v13)
        normals_norm = np.linalg.norm(v123, axis=1)
        inan = np.where(normals_norm == 0)[0]
        return v123, normals_norm, inan

    def remove_elements_with_bad_normals(self, elements, nodes=None):
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=nodes)
        if len(inan):
            inotnan = np.where(normals_norm != 0)[0]
            self.elements = elements[inotnan, :]
            normals_norm = normals_norm[inotnan]
            v123 = v123[inotnan]
            self.log.info('removing %i elements with coincident nodes' % len(inan))

        normals = v123
        normals[:, 0] /= normals_norm
        normals[:, 1] /= normals_norm
        normals[:, 2] /= normals_norm
        return normals

    def get_area(self, elements, stop_on_failure=True):
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=self.nodes)

        if stop_on_failure:
            msg = 'Failed Elements: %s\n' % inan
            if len(inan) > 0:
                for inani in inan:
                    msg += '  eid=%s nodes=%s\n' % (inani, elements[inani, :])
                    for ni in elements[inani]:
                        msg += '    nid=%s node=%s\n' % (ni, self.nodes[ni, :])
                raise RuntimeError(msg)
        return 0.5 * normals_norm

    def get_normals(self, elements, nodes=None, stop_on_failure=True):
        """
        Parameters
        ----------
        elements : (n, 3) ndarray ints
            the elements to get the normals for
        nodes : (n, ) ndarray; default=None -> all
            a subset of the nodes
        stop_on_failure : bool (default=True)
            True:  crash if there are coincident points
            False: delete elements
        """
        if nodes is None:
            nodes = self.nodes
        v123, normals_norm, inan = self._get_normals_data(elements, nodes=nodes)

        if stop_on_failure:
            msg = 'Failed Elements: %s\n' % inan
            if len(inan) > 0:

                for ifail, inani in enumerate(inan):
                    msg += '  eid=%s nodes=%s\n' % (inani, elements[inani, :])
                    for ni in elements[inani]:
                        msg += '    nid=%s node=%s\n' % (ni, nodes[ni, :])
                    if ifail > 10:
                        break
                msg += 'Failed Elements: %s; n=%s\n' % (inan, len(inan))
                raise RuntimeError(msg)
            # we need to divide our (n,3) array in 3 steps
            normals = v123 # / normals_norm
            normals[:, 0] /= normals_norm
            normals[:, 1] /= normals_norm
            normals[:, 2] /= normals_norm
        else:
            inotnan = np.where(normals_norm != 0.)[0]
            inan = np.where(normals_norm == 0.)[0]
            normals_norm[inan] = np.array([1., 0., 0.])
            #elements = elements[inotnan, :]
            #normals_norm = normals_norm[inotnan]
            #v123 = v123[inotnan]

            # we need to divide our (n,3) array in 3 steps
            if 0:
                normals = v123 # / normals_norm
                normals[:, 0] /= normals_norm
                normals[:, 1] /= normals_norm
                normals[:, 2] /= normals_norm
                normals[inan, :] = -1.01
            else:
                normals = v123 # / normals_norm
                normals[inotnan, 0] /= normals_norm[inotnan]
                normals[inotnan, 1] /= normals_norm[inotnan]
                normals[inotnan, 2] /= normals_norm[inotnan]

        return normals


    def flip_normals(self, i=None):
        """
        Flips the normals of the specified elements.

        Parameters
        ----------
        i : (n, ) ndarray ints; default=None -> all
            the indicies to flip
        """
        self.log.info("---flip_normals---")
        if i is None:
            elements = self.elements
        else:
            elements = self.elements[i, :]

        n0, n1, n2 = elements[:, 0], elements[:, 1], elements[:, 2]
        elements2 = elements.copy()
        elements2[:, 0] = n0
        elements2[:, 1] = n2
        elements2[:, 2] = n1
        if i is None:
            self.elements = elements2
        else:
            self.elements[i, :] = elements2 #[i, :]

    def get_normals_at_nodes(self, elements, normals=None, nid_to_eid=None):
        """
        Calculates the normal vector of the nodes based on the average
        element normal.

        Parameters
        ----------
        elements : ????
            The elements...should be removed
        normals : (n, 3) ndarray floats
            The elemental normals
        nid_to_eid : Dict[int] = [int, int, ... ]
            key = node_id
            value = list of element_ids

        Returns
        -------
        normals_at_nodes : (nnodes, 3) ndarray ints
            the normals
        """
        if normals is None:
            nodes = self.nodes
            normals = self.get_normals(elements, nodes=self.nodes)

        if nid_to_eid is None:
            from collections import defaultdict
            nid_to_eid = defaultdict(list)
            eid = 0
            for (n1, n2, n3) in elements:
                nid_to_eid[n1].append(eid)
                nid_to_eid[n2].append(eid)
                nid_to_eid[n3].append(eid)
                eid += 1
            del eid, n1, n2, n3

        normals_at_nodes = np.zeros(nodes.shape, 'float64')
        eid = 0
        for nid, elementsi in iteritems(nid_to_eid):
            pe = normals[elementsi]
            m = pe.mean(axis=0)
            normals_at_nodes[nid] = m / np.linalg.norm(m)
            eid += 1
        return normals_at_nodes

    def equivalence_nodes(self, tol=1e-5):
        nnodes = self.nodes.shape[0]

        # build the kdtree
        kdt = scipy.spatial.KDTree(self.nodes)

        # find the node ids of interest
        nids_new = np.unique(self.elements.ravel())
        nids_new.sort()

        # check the closest 10 nodes for equality
        deq, ieq = kdt.query(self.nodes[nids_new, :], k=10, distance_upper_bound=tol)

        # get the ids of the duplicate nodes
        slots = np.where(ieq[:, 1:] < nnodes)
        replacer = np.unique(ieq[slots])

        # update the duplcated node id with it's partner
        # we'll pick the minimum ID
        for r in replacer:
            ip = np.where(ieq[r, :] < nnodes)[0]
            possible = ieq[r, ip]

            # node 11 can become node 10, but node 10 cannot become node 11
            ip2 = np.where(r > possible)[0]

            if len(ip2):
                # replace the node ids
                possible2 = possible[ip2]
                r_new_nid = possible2.min()
                ireplace = np.where(self.elements == r)
                self.elements[ireplace] = r_new_nid

    def project_boundary_layer(self, nodes, elements, volume_bdfname):
        """
        Create a boundary layer mesh.

        Parameters
        ----------
        nodes : (n, 3) ndarray floats
            The nodes on the surface.
        elements : (n, 3) ndarray ints
            The elements on the surface.
        volume_bdfname : str
            The CPENTA bdf file to write.

        Returns
        -------
        nodes2 : (n, 3) ndarray floats
            The boundary layer nodes
        elements2 : (n, 6) ndarray ints
            The boundary layer elements
        """
        self.log.info("project_mesh...")

        normals_at_nodes = self.get_normals_at_nodes(elements, normals=None, nid_to_eid=None)

        #print "normals_at_nodes[4]", normals_at_nodes[4]
        #----------- make boundary layer---------------
        # deltaN = a^N * delta
        delta = 0.1
        b = 1.0
        a = 1.1
        N = 13
        r = np.array(range(N))
        #r = 1100

        delta_ns = b * a**r * delta
        self.log.info('delta_ns = %s' % delta_ns)
        if not isinstance(delta_ns, np.ndarray):
            delta_ns = np.array([delta_ns])
        N = len(delta_ns)

        nid = 0
        nnodes = nodes.shape[0]
        nodes2 = np.zeros((nnodes * (N+1), 3), 'float64')
        nodes2[:nnodes, :] = nodes

        nelements = elements.shape[0]
        elements2 = np.zeros((nelements * (N+1), 3), 'int32')
        elements2[:nelements, :] = elements

        ni = 0
        cid = None
        nid = 1
        eid2 = 1
        pid = 100
        mid = 100
        with open(volume_bdfname, 'wb') as bdf:
            bdf.write('CEND\nBEGIN BULK\n')
            bdf.write('$NODES in Layer=0\n')
            for (x, y, z) in nodes:
                card = ['GRID', nid, cid, x, y, z]
                bdf.write(print_card_8(card))
                nid += 1

            for deltaN in delta_ns:
                outer_points = nodes + normals_at_nodes * deltaN
                nodes2[ni*nnodes : (ni+1)*nnodes, :] = outer_points

                nnbase = ni * nnodes
                nnshift = (ni+1) * nnodes

                nebase = (ni) * nelements
                neshift = (ni + 1) * nelements
                elements2[neshift : neshift + nelements, :] = elements + nnodes * (ni + 1)

                self.log.info('nodes = %s' % str(nodes))
                self.log.info('deltaN = %s' % str(deltaN))
                self.log.info('normals_at_nodes = %s' % str(normals_at_nodes))
                self.log.info('outer_points = %s' % str(outer_points))
                bdf.write('$NODES in Layer=%i\n' % (ni + 1))
                for x, y, z in outer_points:
                    card = ['GRID', nid, cid, x, y, z]
                    bdf.write(print_card_8(card))
                    nid += 1

                bdf.write('$SOLID ELEMENTS in Layer=%i\n' % (ni + 1))
                for eid in range(nelements):
                    (n1, n2, n3) = elements2[nebase  + eid] + 1
                    (n4, n5, n6) = elements2[neshift + eid] + 1
                    card = ['CPENTA', eid2, pid, n1, n2, n3, n4, n5, n6]
                    bdf.write(print_card_8(card))
                    eid2 += 1

                card = ['PSOLID', pid, mid]
                bdf.write(print_card_8(card))

                E = 1e7
                G = None
                nu = 0.3
                card = ['MAT1', mid, E, G, nu]
                bdf.write(print_card_8(card))

                pid += 1
                mid += 1
                ni += 1

        #print(elements2)
        #for node in nodes:
            #normal = normals_elements[nid]
            #nid += 1
        #print(deltaN)

        #----------- make far field---------------
        nodes3 = nodes2[nnbase:, :]
        nodes3 = nodes2[nnbase:, :]

        elements3 = elements2[nebase:, :]
        elements3 = elements2[nebase:, :]
        self.log.debug("done projecting...")
        return nodes2, elements2


    def write_stl_ascii(self, out_filename, solid_name, float_fmt='%.6f',
                        normalize_normal_vectors=False, stop_on_failure=True):
        """
        solid solid_name
         facet normal -6.601157e-001 6.730213e-001 3.336009e-001
          outer loop
            vertex 8.232952e-002 2.722531e-001 1.190414e+001
            vertex 8.279775e-002 2.717848e-001 1.190598e+001
            vertex 8.557653e-002 2.745033e-001 1.190598e+001
          endloop
         endfacet
        end solid
        """
        nodes = self.nodes
        elements = self.elements
        self.log.info("---write_stl_ascii...%r---" % out_filename)
        msg = ''
        node_format = ' facet normal %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        vertex_format = '     vertex %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        msg += 'solid %s\n' % solid_name
        with open(out_filename, 'w') as out:
            out.write(msg)

            nelements = elements.shape[0]
            normals = self.get_normals(elements, stop_on_failure=stop_on_failure)
            for element, normal in zip(elements, normals):
                try:
                    p1, p2, p3 = nodes[element]
                except IndexError:
                    print(element)
                    raise

                #msg  += '  facet normal -6.601157e-001 6.730213e-001 3.336009e-001\n'
                msg = node_format % tuple(normal)
                msg += '   outer loop\n'
                msg += vertex_format % tuple(p1)
                msg += vertex_format % tuple(p2)
                msg += vertex_format % tuple(p3)
                msg += '   endloop\n'
                msg += ' endfacet\n'
                #print(msg)
                out.write(msg)
            msg = 'endsolid\n'
            out.write(msg)


    def read_ascii_stl(self, stl_filename):
        infile = open(stl_filename, 'r')
        line = infile.readline()
        inode = 0
        ielement = 0
        nodes_dict = {}
        elements = []
        while line:
            if 'solid' in line[:6].lower():
                line = infile.readline()  # facet
                while 'facet' in line.strip()[:5].lower():
                    #facet normal -6.665299e-001 6.795624e-001 3.064844e-001
                    #   outer loop
                    #      vertex 8.142845e-002 2.731541e-001 1.190024e+001
                    #      vertex 8.186898e-002 2.727136e-001 1.190215e+001
                    #      vertex 8.467505e-002 2.754588e-001 1.190215e+001
                    #   endloop
                    #endfacet

                    infile.readline() # outer loop
                    L1 = infile.readline().strip()
                    L2 = infile.readline().strip()
                    L3 = infile.readline().strip()

                    v1 = L1.split()[1:]
                    v2 = L2.split()[1:]
                    v3 = L3.split()[1:]
                    infile.readline() # endloop
                    infile.readline() # endfacet
                    t1 = tuple(v1)
                    t2 = tuple(v2)
                    t3 = tuple(v3)

                    assert len(v1) == 3, '%r' % L1
                    assert len(v2) == 3, '%r' % L2
                    assert len(v3) == 3, '%r' % L3

                    if t1 in nodes_dict:
                        i1 = nodes_dict[t1]
                    else:
                        i1 = inode
                        nodes_dict[t1] = inode
                        inode += 1

                    if t2 in nodes_dict:
                        i2 = nodes_dict[t2]
                    else:
                        i2 = inode
                        nodes_dict[t2] = inode
                        inode += 1

                    if t3 in nodes_dict:
                        i3 = nodes_dict[t3]
                    else:
                        i3 = inode
                        nodes_dict[t3] = inode
                        inode += 1
                    element = [i1, i2, i3]
                    elements.append(element)
                    ielement += 1
                    line = infile.readline()  # facet
                #print "end of solid..."
            elif 'endsolid' in line.lower():
                line = infile.readline()
            elif line.strip() == '':
                line = infile.readline()
            else:
                self.log.error(line)
                #line = f.readline()
                raise NotImplementedError('multiple solids are not supported; line=%r' % line)
                #break
        infile.close()

        assert inode > 0, inode
        nnodes = inode + 1 # accounting for indexing
        self.elements = np.array(elements, 'int32')
        nodes = np.zeros((nnodes, 3), 'float64')

        for node, inode in iteritems(nodes_dict):
            nodes[inode] = node
        self.nodes = nodes

    def scale_nodes(self, xscale, yscale, zscale):
        self.nodes[:, 0] *= xscale
        self.nodes[:, 1] *= yscale
        self.nodes[:, 2] *= zscale

    def shift_nodes(self, xshift, yshift, zshift):
        self.nodes[:, 0] += xshift
        self.nodes[:, 1] += yshift
        self.nodes[:, 2] += zshift

    def flip_axes(self, axes, scale):
        if axes == 'xy':
            x = copy.deepcopy(self.nodes[:, 0])
            y = copy.deepcopy(self.nodes[:, 1])
            self.nodes[:, 0] = y * scale
            self.nodes[:, 1] = x * scale
        elif axes == 'yz':
            y = copy.deepcopy(self.nodes[:, 1])
            z = copy.deepcopy(self.nodes[:, 2])
            self.nodes[:, 1] = z * scale
            self.nodes[:, 2] = y * scale
        elif axes == 'xz':
            x = copy.deepcopy(self.nodes[:, 0])
            z = copy.deepcopy(self.nodes[:, 2])
            self.nodes[:, 0] = z * scale
            self.nodes[:, 2] = x * scale

    def create_mirror_model(self, xyz, tol):
        """
        Creates a mirror model.

        Parameters
        ----------
        xyz : str {x, y, z}
            the direction of symmetry
        tol: float
            the tolerance for symmetry plane nodes

        .. note:: All elements on the symmetry plane will be removed
        """
        assert xyz in ['x', 'y', 'z'], 'xyz=%r' % xyz
        assert tol >= 0.0, 'tol=%s' % tol

        nnodes = self.nodes.shape[0]
        if xyz == 'x':
            xyzi = 0
        elif xyz == 'y':
            xyzi = 1
        elif xyz == 'z':
            xyzi = 2
        else:
            raise RuntimeError(xyz)

        # the nodes on the symmetry plane
        i = np.where(self.nodes[:, xyzi] < tol)[0]

        # smash the symmetry nodes to 0.0
        self.nodes[i, xyzi] = 0.
        nodes_sym = copy.deepcopy(self.nodes)
        nodes_sym[:, xyzi] *= -1.

        # we're lazy and duplicating all the nodes
        # but will only write out a subset of them
        nodes = np.vstack([self.nodes, nodes_sym])

        # create the symmetrical elements
        elements2 = []
        elements3 = []
        for element in self.elements:
            epoints = nodes[element, xyzi][0]
            je = np.where(epoints <= tol)[0]
            if len(je) < 3:  # not a symmetry element, so we save it
                elements2.append(element)

                # duplicate the node if it's not on the symmetry plane
                element3 = [elementi if elementi in i else (elementi + nnodes)
                            for elementi in element]

                # the normal is now backwards, so we flip it
                element3.reverse()
                elements3.append(element3)

        self.nodes = nodes
        self.elements = np.array(elements2 + elements3, dtype='int32')


def _rotate_model(stl):  # pragma: no cover
    nodes = stl.nodes
    elements = stl.elements
    if 0:
        # rotate the model
        x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
        #i = where(y > 0.0)[0]
        R = x**2 + y**2
        theta = np.arctan2(y, x)
        iRz = np.where(R == 0)[0]
        theta[iRz] = 0.0

        min_theta = min(theta)
        dtheta = max(theta) - np.pi / 4
        theta2 = theta + min_theta

        x2 = R * np.cos(theta2)
        y2 = R * np.sin(theta2)
        #print("x.shape", x.shape, y2.shape)
        nodes_rotated = np.transpose(np.vstack([x2, y2, z]))
        #print("nodes.shape", nodes_rotated.shape)
        #print(nodes_rotated)

    if 0:
        # project the volume
        (nodes2, elements2) = stl.project_mesh(nodes_rotated, elements)

    # write the model
    stl_geom_out = 'rotated.stl'
    stl.write_stl_ascii(stl_geom_out, 'sphere')


if __name__ == '__main__':  # pragma: no cover
    from pyNastran.converters.stl.stl_reshape import main
    main()

