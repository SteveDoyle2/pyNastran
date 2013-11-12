#pylint:  disable=C0103,C0111
import os
import sys
#from struct import pack
from itertools import izip

from numpy import array, zeros, ndarray, cross
from numpy.linalg import norm

#from struct import unpack
#from pyNastran.bdf.fieldWriter import print_card
from pyNastran.utils import is_binary
from pyNastran.utils.log import get_logger


class STLReader(object):
    #modelType = 'cart3d'
    #isStructured = False
    #isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')

        #self.cartType = None  # grid, results
        #self.nPoints = None
        #self.nElements = None
        #self.infile = None
        #self.infilename = None
        #self.log = get_logger(log, 'debug' if debug else 'info')

    def write_stl(self, outfilename, points, elements, regions, loads=None, is_binary=False, float_fmt='%6.7f'):
        assert len(points) > 0

        if loads is None or loads == {}:
            loads = {}
            is_loads = False
            print "no loads"
        else:
            is_loads = True

        self.log.info("---writing STL file...%r---" % outfilename)
        f = open(outfilename, 'wb')

        int_fmt = self.write_header(f, points, elements, is_loads, is_binary)
        #print "int_fmt =", int_fmt
        self.write_points(f, points, is_binary, float_fmt)
        self.write_elements(f, elements, is_binary, int_fmt)
        self.write_regions(f, regions, is_binary)

        if is_loads:
            assert is_binary is False, 'is_binary=%r is not supported for loads' % is_binary
            self.write_loads(f, loads, is_binary, float_fmt)
        f.close()

    def read_stl(self, stl_filename):
        self.infilename = stl_filename
        self.log.info("---starting reading STL file...%r---" % self.infilename)

        if is_binary(stl_filename):
            print "***is_binary"
            aaa
            self.infile = open(stl_filename, 'rb')
            (self.nPoints, self.nElements) = self.read_header_binary()
            points = self.read_points_binary(self.nPoints)
            elements = self.read_elements_binary(self.nElements)
            regions = self.read_regions_binary(self.nElements)
            loads = {}
        else:
            self.infile = open(stl_filename, 'r')
            nodes, elements = self.read_ascii()

        self.infile.close()
        #self.log.info("nPoints=%s  nElements=%s" % (self.nPoints, self.nElements))
        self.log.info("---finished reading STL file...%r---" % self.infilename)
        #assert self.nPoints > 0, 'nPoints=%s' % self.nPoints
        #assert self.nElements > 0, 'nElements=%s' % self.nElements
        return nodes, elements

    def get_normals(self, nodes, elements):
        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)
        normals = v123 / norm(v123)
        return normals

    def project_mesh(self, nodes, elements):
        print "project_mesh..."

        from collections import defaultdict
        nid_to_eid = defaultdict(list)

        eid = 0
        for (n1, n2, n3) in elements:
            #n1, n2, n3 = element
            #print n1, n2, n3
            nid_to_eid[n1].append(eid)
            nid_to_eid[n2].append(eid)
            nid_to_eid[n3].append(eid)
            eid += 1
        del eid, n1, n2, n3
        print "nid_to_eid[0] =", nid_to_eid[0]

        normals = self.get_normals(nodes, elements)
        normals_at_nodes = zeros(nodes.shape, 'float64')
        print "normals_elements =", normals_at_nodes.shape
        eid = 0
        for nid, elementsi in nid_to_eid.iteritems():
            pe = normals[elementsi]
            m = pe.mean(axis=0)
            normals_at_nodes[nid] = m/norm(m)
            eid += 1
            #print "elementsi", elementsi
            #print "pe =", pe
            #print "m =", m

        print "normals_at_nodes[4]", normals_at_nodes[4]
        #----------- make boundary layer---------------
        # deltaN = a^N * delta
        delta = 0.01
        N = 5
        b = 1.0
        a = 1.01
        #r = array(range(N))
        r = 1100

        deltaNs = b * a**r * delta
        if not isinstance(deltaNs, ndarray):
            deltaNs = array([deltaNs])
        N = len(deltaNs)

        nid = 0
        print "deltaNs =", deltaNs
        nnodes, three = nodes.shape
        nodes2 = zeros((nnodes * (N+1), 3), 'float64')
        nodes2[:nnodes, :] = nodes

        nelements, three = elements.shape
        elements2 = zeros((nelements * (N+1), 3), 'int32')
        elements2[:nelements, :] = elements
        print "elements2.shape =", elements2.shape

        ni = 1
        print "nelements =", nelements
        for deltaN in deltaNs:
            print "  deltaNi =", deltaN
            outer_points = nodes + normals_at_nodes * deltaN
            print "  outer_points.shape", outer_points.shape
            nodes2[   ni*nnodes   :(ni+1)*nnodes, :] = outer_points
            print "  %s:%s" % (ni*nelements, (ni+1)*nelements)
            elements2[ni*nelements:(ni+1)*nelements, :] = elements + nelements
            ni += 1
        print elements2
        #for node in nodes:
            #normal = normals_elements[nid]
            #nid += 1
        print deltaN

        #----------- make far field---------------
        print "done projecting..."


    def _write_stl_ascii(self, out_filename, solid_name, nodes, elements):
        """
        facet normal -6.601157e-001 6.730213e-001 3.336009e-001
           outer loop
              vertex 8.232952e-002 2.722531e-001 1.190414e+001
              vertex 8.279775e-002 2.717848e-001 1.190598e+001
              vertex 8.557653e-002 2.745033e-001 1.190598e+001
           endloop
        endfacet
        """
        msg = ''
        form = '%.6e'
        nFormat = '  facet normal %s %s %s\n' % (form, form, form)
        vFormat = '    vertex %s %s %s\n' % (form, form, form)
        msg += 'solid %s\n' % solid_name
        out = open(out_filename, 'wb')
        out.write(msg)
        for element in elements:
            #msg = ''
            p1, p2, p3 = nodes[element]
            v12 = p2 - p1
            v13 = p3 - p1
            v123 = cross(v12, v13)
            normal = v123 / norm(v123)

            #msg  += '  facet normal -6.601157e-001 6.730213e-001 3.336009e-001\n'
            msg = nFormat % tuple(normal)
            #print "p1 =", p1
            #print "p2 =", p2
            #print "p3 =", p3
            msg += '    outer loop\n'
            msg += vFormat % tuple(p1)
            msg += vFormat % tuple(p2)
            msg += vFormat % tuple(p3)
            msg += '    endloop\n'
            msg += '  endfacet\n'
            #print msg
            out.write(msg)
        msg = 'endsolid'
        out.write(msg)


    def read_ascii(self):
        f = self.infile
        line = f.readline()
        inode = 0
        ielement = 0
        nodes_dict = {}
        elements = []
        while line:
            if 'solid ' in line[:6]:
                line = f.readline()  # facet
                #print "line = %r" % line
                while 'facet' in line.strip()[:5]:
                    #facet normal -6.665299e-001 6.795624e-001 3.064844e-001
                    #   outer loop
                    #      vertex 8.142845e-002 2.731541e-001 1.190024e+001
                    #      vertex 8.186898e-002 2.727136e-001 1.190215e+001
                    #      vertex 8.467505e-002 2.754588e-001 1.190215e+001
                    #   endloop
                    #endfacet

                    #f.readline() # facet
                    f.readline() # outer loop
                    L1 = f.readline().strip()
                    L2 = f.readline().strip()
                    L3 = f.readline().strip()

                    v1 = L1.split()[1:]
                    v2 = L2.split()[1:]
                    v3 = L3.split()[1:]
                    f.readline() # endloop
                    f.readline() # endfacet
                    #v1 = array(v1)
                    #v2 = array(v2)
                    #v3 = array(v3)
                    t1 = tuple(v1)
                    t2 = tuple(v2)
                    t3 = tuple(v3)
                    #print "%s\n%s\n%s" % (t1, t2, t3)
                    assert len(v1) == 3, '%r' % L1
                    assert len(v2) == 3, '%r' % L2
                    assert len(v3) == 3, '%r' % L3
                    #print "--------------"

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
                    #print "element =", inode, ielement, element
                    elements.append(element)
                    ielement += 1
                    line = f.readline()  # facet
                #print "ielement =", ielement
                nnodes = inode + 1 # accounting for indexing
                #nelements = ielement
                elements = array(elements, 'int32')
                nodes = zeros((nnodes, 3), 'float64')
                #print "end of solid..."
                #print "*line = %r" % line
            else:
                break



        print "end of while looop..."
        #print "nnodes = ", nnodes
        #print "nodes =", nodes.shape
        #print "elements =", elements.shape
        for node, inode in nodes_dict.iteritems():
            nodes[inode] = node
        #print "***"
        return nodes, elements


if __name__ == '__main__':
    import time
    t0 = time.time()

    # binary
    #cart3dGeom = os.path.join('flat.tri')  # half
    #full_model = os.path.join('flat_full.tri')

    # ascii
    stl_geom = 'spw_half.STL'

    log = None
    debug = False
    model = STLReader(log, debug)  # ascii/binary

    (nodes, elements) = model.read_stl(stl_geom)
    model.project_mesh(nodes, elements)

    model._write_stl_ascii('spw_half.out.STL', 'spw_half', nodes, elements)
