#pylint:  disable=C0103,C0111
import os
import sys
#from struct import pack
from itertools import izip

from docopt import docopt

from numpy import array, zeros, ndarray, cross, where
from numpy.linalg import norm

#from struct import unpack
import pyNastran
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.utils import is_binary
from pyNastran.utils.log import get_logger


class STLReader(object):
    #modelType = 'stl'
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

    def write_stl(self, out_filename, is_binary=False, float_fmt='%6.7f'):
        self.log.info("---writing STL file...%r---" % out_filename)
        assert len(self.nodes) > 0
        solid_name = 'dummy_name'
        if is_binary:
            raise NotImplementedError('binary writing is not supported')
        else:
            self._write_stl_ascii(out_filename, solid_name, float_fmt=float_fmt)

    def read_stl(self, stl_filename):
        self.infilename = stl_filename
        self.log.info("---starting reading STL file...%r---" % self.infilename)

        if is_binary(stl_filename):
            self.log.debug("***is_binary")
            aaa
            self.infile = open(stl_filename, 'rb')
            (self.nPoints, self.nElements) = self.read_header_binary()
            points = self.read_points_binary(self.nPoints)
            elements = self.read_elements_binary(self.nElements)
            regions = self.read_regions_binary(self.nElements)
            loads = {}
        else:
            self.infile = open(stl_filename, 'r')
            self.read_ascii()

        self.infile.close()
        #self.log.info("nPoints=%s  nElements=%s" % (self.nPoints, self.nElements))
        self.log.info("---finished reading STL file...%r---" % self.infilename)
        #assert self.nPoints > 0, 'nPoints=%s' % self.nPoints
        #assert self.nElements > 0, 'nElements=%s' % self.nElements

    def get_normals(self, elements, nodes=None):
        if nodes is None:
            nodes = self.nodes
        self.log.debug("get_normals...elements.shape %s" % str(elements.shape))
        p1 = nodes[elements[:, 0]]
        p2 = nodes[elements[:, 1]]
        p3 = nodes[elements[:, 2]]
        v12 = p2 - p1
        v13 = p3 - p1
        v123 = cross(v12, v13)

        # verifies that y normals don't point down (into the model)
        #if min(v123[1]) < 0.0:
            #iy_zero = where(v123[1] < 0)[0]
            #raise RuntimeError(iy_zero+1)
        #print "v123", v123
        n = norm(v123, axis=1)
        inan = where(n==0)[0]
        if len(inan) > 0:
            raise RuntimeError(inan)
        #from numpy import divide

        # we need to divide our (n,3) array in 3 steps
        normals = v123 # /n
        normals[:, 0] /= n
        normals[:, 1] /= n
        normals[:, 2] /= n
        #divide(v123, n)
        return normals

    def flip_normals(self):
        self.log.info("---flip_normals...%r---" % self.infilename)
        elements = self.elements
        n0, n1, n2 = elements[:, 0], elements[:, 1], elements[:, 2]
        elements2 = elements.copy()
        elements2[:, 0] = n0
        elements2[:, 1] = n2
        elements2[:, 2] = n1
        self.elements = elements

    def get_normals_at_nodes(self, elements, normals=None, nid_to_eid=None):
        if normals is None:
            nodes = self.nodes
            normals = self.get_normals(elements, nodes=self.nodes)

        if nid_to_eid is None:
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

        normals_at_nodes = zeros(nodes.shape, 'float64')
        eid = 0
        for nid, elementsi in nid_to_eid.iteritems():
            pe = normals[elementsi]
            m = pe.mean(axis=0)
            normals_at_nodes[nid] = m/norm(m)
            eid += 1
            #print "elementsi", elementsi
            #print "pe =", pe
            #print "m =", m
        return normals_at_nodes

    def project_boundary_layer(self, nodes, elements, volume_bdfname):
        """
        create a boundary layer mesh
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
        r = array(range(N))
        #r = 1100

        deltaNs = b * a**r * delta
        print 'deltaNs =', deltaNs
        if not isinstance(deltaNs, ndarray):
            deltaNs = array([deltaNs])
        N = len(deltaNs)
        #print "N = ", N

        nid = 0
        #print "deltaNs =", deltaNs
        nnodes, three = nodes.shape
        nodes2 = zeros((nnodes * (N+1), 3), 'float64')
        nodes2[:nnodes, :] = nodes

        nelements, three = elements.shape
        elements2 = zeros((nelements * (N+1), 3), 'int32')
        elements2[:nelements, :] = elements
        #print "nodes.shape =", nodes.shape
        #print "nodes2.shape =", nodes2.shape
        #print "elements2.shape =", elements2.shape

        ni = 0
        #print "nelements =", nelements
        cid = None
        nid = 1
        eid2 = 1
        pid = 100
        mid = 100
        bdf = open(volume_bdfname, 'wb')
        bdf.write('CEND\nBEGIN BULK\n')
        bdf.write('$NODES in Layer=0\n')
        for (x, y, z) in nodes:
            card = ['GRID', nid, cid, x, y, z]
            bdf.write(print_card(card))
            nid += 1

        for deltaN in deltaNs:
            outer_points = nodes + normals_at_nodes * deltaN
            nodes2[   ni*nnodes   :(ni+1)*nnodes, :] = outer_points

            nnbase = ni * nnodes
            nnshift = (ni+1) * nnodes

            nebase  = (ni) * nelements
            neshift = (ni + 1) * nelements
            elements2[neshift : neshift + nelements, :] = elements + nnodes * (ni + 1)

            print('nodes = %s' % str(nodes))
            print('deltaN = %s' % str(deltaN))
            print('normals_at_nodes = %s' % str(normals_at_nodes))
            print('outer_points = %s' % str(outer_points))
            bdf.write('$NODES in Layer=%i\n' % (ni + 1))
            for (x, y, z) in outer_points:
                card = ['GRID', nid, cid, x, y, z]
                bdf.write(print_card(card))
                nid += 1

            bdf.flush()
            bdf.write('$SOLID ELEMENTS in Layer=%i\n' % (ni + 1))
            for eid in xrange(nelements):
                (n1, n2, n3) = elements2[nebase  + eid] + 1
                (n4, n5, n6) = elements2[neshift + eid] + 1
                card = ['CPENTA', eid2, pid, n1, n2, n3, n4, n5, n6]
                bdf.write(print_card(card))
                eid2 += 1
                bdf.flush()

            card = ['PSOLID', pid, mid]
            bdf.write(print_card(card))

            E = 1e7
            G = None
            nu = 0.3
            card = ['MAT1', mid, E, G, nu]
            bdf.write(print_card(card))

            pid += 1
            mid += 1
            ni += 1
        bdf.close()
        #print elements2
        #for node in nodes:
            #normal = normals_elements[nid]
            #nid += 1
        #print deltaN

        #----------- make far field---------------
        nodes3 = nodes2[nnbase : , :]
        nodes3 = nodes2[nnbase : , :]

        elements3 = elements2[nebase : , :]
        elements3 = elements2[nebase : , :]
        self.log.debug("done projecting...")
        return nodes2, elements2


    def _write_stl_ascii(self, out_filename, solid_name, float_fmt='%.6f'):
        """
        facet normal -6.601157e-001 6.730213e-001 3.336009e-001
           outer loop
              vertex 8.232952e-002 2.722531e-001 1.190414e+001
              vertex 8.279775e-002 2.717848e-001 1.190598e+001
              vertex 8.557653e-002 2.745033e-001 1.190598e+001
           endloop
        endfacet
        """
        nodes = self.nodes
        elements = self.elements
        self.log.info("---_write_stl_ascii...%r---" % out_filename)
        msg = ''
        nFormat = '  facet normal %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        vFormat = '    vertex %s %s %s\n' % (float_fmt, float_fmt, float_fmt)
        msg += 'solid %s\n' % solid_name
        out = open(out_filename, 'wb')
        out.write(msg)

        nelements, three = elements.shape
        #print "write nelements=%s" % (nelements)
        for element in elements:
            #msg = ''
            #p1 = nodes[element[0]]
            #p2 = nodes[element[1]]
            #p3 = nodes[element[2]]
            try:
                p1, p2, p3 = nodes[element]
            except IndexError:
                print element
                raise
            #print p1, p2, p3
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
            print "**line = %r" % line
            if 'solid' in line[:6]:
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
                #nelements = ielement
                #print "end of solid..."
                #print "*line = %r" % line
            else:
                break

        assert inode > 0
        nnodes = inode + 1 # accounting for indexing
        self.elements = array(elements, 'int32')
        nodes = zeros((nnodes, 3), 'float64')

        #print "end of while loop..."
        #print "nnodes = ", nnodes
        #print "nodes =", nodes.shape
        #print "elements =", elements.shape
        for node, inode in nodes_dict.iteritems():
            nodes[inode] = node
        self.nodes = nodes
        #print "***"


def run_arg_parse():
    msg  = 'This program flips the normal of an STL model.\n'
    msg += 'Usage:\n'
    msg += '  stl_reader INPUT [-o OUTPUT]\n'
    msg += '             [-n] [-q]\n'
    msg += '  stl_reader -h | --help\n'
    msg += '  stl_reader -v | --version\n'
    msg += "  INPUT      path to input file\n"
    msg += "\n"
    msg += "Options:\n"
    msg += "  -h, --help                  show this help message and exit\n"
    #msg += "  -f FORMAT, --format FORMAT  format type (panair, cart3d,\n"
    #msg += "                                           nastran, lawgs, stl)\n"
    msg += "  -o OUTPUT, --output OUTPUT  path to output file\n"
    msg += "  -n, --normal                flip the element normals\n"
    #msg += "  -r XYZ, --rotation XYZ      [x, y, z, -x, -y, -z] default is ???\n"

    msg += "  -q, --quiet                 prints debug messages (default=True)\n"
    msg += "  -v, --version               show program's version number and exit\n"

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    #print data

    #format  = data['--format']
    input = data['INPUT']
    print "input =", input
    output = data['--output']
    if output is None:
        input_base, ext = os.path.splitext(input)
        output = input_base + '_out.stl'

    reverse_normals = data['--normal']
    assert reverse_normals in [True, False]
    quiet = data['--quiet']
    return (input, output, reverse_normals, quiet)

def main():
    (stl_geom_in, stl_geom_out, reverse_normals, quiet) = run_arg_parse()
    import time
    t0 = time.time()
    # binary
    #stl_geom = None

    # ascii
    #stl_geom = 'spw_half.STL'

    log = None
    debug = False
    model = STLReader(log, debug)  # ascii/binary

    (nodes, elements) = model.read_stl(stl_geom_in)
    if reverse_normals:
        elements = model.flip_normals(elements)


    if 0:
        # rotate the model
        from numpy import where, arctan2, cos, sin, hstack, vstack, concatenate, transpose, pi
        x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
        #i = where(y > 0.0)[0]
        R = x**2 + y**2
        theta = arctan2(y, x)
        iRz = where(R == 0)[0]
        theta[iRz] = 0.0

        min_theta = min(theta)
        dtheta = max(theta) - pi/4
        theta2 = theta + min_theta

        x2 = R * cos(theta2)
        y2 = R * sin(theta2)
        #print "x.shape", x.shape, y2.shape
        nodes_rotated = transpose(vstack([x2, y2, z]))
        #print "nodes.shape", nodes_rotated.shape
        #print nodes_rotated

    if 0:
        # project the volume
        (nodes2, elements2) = model.project_mesh(nodes_rotated, elements)

    # write the model
    model._write_stl_ascii(stl_geom_out, 'sphere', nodes, elements)

if __name__ == '__main__':
    main()