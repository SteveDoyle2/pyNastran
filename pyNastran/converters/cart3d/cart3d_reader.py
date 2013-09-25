#pylint:  disable=C0103,C0111
import os
import sys
from struct import pack
from math import ceil#, sqrt
from itertools import izip

from numpy import array, zeros, where, savetxt, sqrt, abs, amax, amin
from numpy import arange, searchsorted, vstack, unique, hstack, ravel

from struct import unpack
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.op2.fortranFile import FortranFile
from pyNastran.utils import is_binary
from pyNastran.utils.log import get_logger


def convert_to_float(svalues):
    values = []
    for value in svalues:
        values.append(float(value))
    return values


class Cart3DReader(object):
    def __init__(self, log=None, debug=False):
        pass

    def read_results(self, i, infile, result_names=None):
        raise NotImplementedError()

    def get_y0_nodes(self, nodes, ax):
        raise RuntimeError('removed...')
        assert ax is 1, 'ax=%s' % ax
        yZeroNodes = {}
        outerNodes = {}
        isInnerNode = {}

        nnodes, three = nodes.shape
        for inode in xrange(nnodes):
            node = nodes[inode, :]
            if node[1] == 0:
                yZeroNodes[inode+1] = node
                isInnerNode[inode+1] = True
            else:
                outerNodes[inode+1] = [node[0], -node[1], node[2]]
                isInnerNode[inode+1] = False
        return (yZeroNodes, outerNodes, isInnerNode)

    def make_mirror_model(self, nodes, elements, regions, loads, axis='y', tol=0.000001):
        """
        Makes a full cart3d model about the given axis.
        
        :param nodes:    the nodes;     (nnodes,    3) ndarray
        :param elements: the elmements; (nelements, 3) ndarray
        :param regions:  the regions;   (nelements)    ndarray
        :param loads:    not supported; dictionary of (nnodes) ndarray
        :param axis:     a string of "x", "y", "z", "-x", "-y", "-z"
        :param tol:      the tolerance for the centerline points (default=0.000001)
        """
        self.log.info('---starting make_mirror_model---')
        assert tol >= 0, 'tol=%r' % tol #  prevents hacks to the axis

        nnodes, three = nodes.shape
        assert nnodes > 0

        nelements, three = elements.shape
        assert nelements > 0

        ax = self._get_ax(axis)
        if ax in [0, 1, 2]:  # positive x, y, z values; mirror to -side
            iy0 = where(nodes[:, ax] > tol)[0]
            ax2 = ax
        elif ax in [3, 4, 5]:  # negative x, y, z values; mirror to +side
            iy0 = where(nodes[:, ax-3] < -tol)[0]
            ax2 = ax - 3 # we create ax2 in order to generalize the data updating
        else:
            raise NotImplementedError(axis)

        # the nodes to be duplicated are the nodes that aren't below the tolerance
        nodes_upper = nodes[iy0]
        nodes_upper[:, ax2] *= -1.0  # flip the nodes about the axis

        nodes2 = vstack([nodes, nodes_upper])
        nnodes2, three = nodes2.shape
        assert nnodes2 > nnodes

        nnodes_upper, three = nodes_upper.shape
        elements_upper = elements.copy()
        nelements, three = elements.shape
        
        # remap the mirrored nodes with the new node ids
        for eid in xrange(nelements):
            element = elements_upper[eid, :]
            for i, eidi in enumerate(element):
                if eidi in iy0:
                    elements_upper[eid][i] = nnodes_upper + eidi
                
                # we need to reverse the element in order to get
                # the proper normal vector
                elements_upper[eid] = elements_upper[eid, ::-1]
            
        elements2 = vstack([elements, elements_upper])
        nelements2, three = elements2.shape
        assert nelements2 > nelements
        
        nregions = len(unique(regions))
        regions_upper = regions.copy() + nregions
        regions2 = hstack([regions, regions_upper])

        loads2 = {}
        for key, data in loads.iteritems():
        
            # flip the sign on the flipping axis terms
            if ((key in ['U', 'rhoU'] and ax2==0) or
                (key in ['V', 'rhoV'] and ax2==1) or
                (key in ['W', 'rhoW'] and ax2==2)):
                data_upper = -data[iy0]
            else:
                data_upper = data[iy0]
            loads2[key] = hstack([data, data_upper])

        self.log.info('---finished make_mirror_model---')
        return (nodes2, elements2, regions2, loads2)

    def _get_ax(self, axis):
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
        else:
            raise NotImplementedError('axis=%r' % axis)
        self.log.info("axis=%r ax=%s" % (axis, ax))
        return ax

    def make_half_model(self, nodes, elements, regions, loads=None, axis='y', remap_nodes=True):
        """
        Makes a half model from a full model
        
        ...note:: Cp is really loads['Cp'] and was meant for loads analysis only
        """
        if loads is None:
            loads = {}

        nnodes, three = nodes.shape
        assert nnodes > 0, 'nnodes=%s'  % nnodes

        nelements, three = elements.shape
        assert nelements > 0, 'nelements=%s'  % nelements

        inodes_remove = set([])
        self.log.info('---starting make_half_model---')
        ax = self._get_ax(axis)

        if ax in [0, 1, 2]:
            inodes_save = where(nodes[:, ax] >= 0.0)[0]
        elif ax in [3, 4, 5]:
            inodes_save = where(nodes[:, ax-3] <= 0.0)[0]
        else:
            raise NotImplementedError('axis=%r ax=%s' % (axis, ax))
        inodes_save.sort()

        inodes_map = arange(len(inodes_save))
        assert 0 < len(inodes_save) < nnodes, 'len(inodes_save)=%s nnodes=%s'  % (len(inodes_save), nnodes)

        nodes2 = nodes[inodes_save, :]
        nnodes2, three = nodes2.shape
        assert 0 < nnodes2 < nnodes, 'nnodes=%s nnodes2=%s'  % (nnodes, nnodes2)
        #print 'nnodes=%s nnodes2=%s'  % (nnodes, nnodes2)
        #print 'inodes_save=%s'  % (inodes_save), len(inodes_save)
        
        inodes_save += 1  # +1 is so we don't have to shift inode
        # ..todo:: still need to handle element's node id renumbering
        ielements_save = set([])
        for ielement in xrange(nelements):
            save_element = True
            element = elements[ielement, :]
            for inode in element:
                # something like this should be faster since we know inodes_save is sorted
                #if inode != inodes_save[searchsorted(inodes_save, inode)]:
                    #print "bad...", ielement, element, inode
                    #aaa
                    #save_element = False
                    #break
                if inode not in inodes_save:
                    #print "bad...", ielement, element, inode
                    save_element = False
                    break

            if save_element:
                #print "saving %s" % ielement
                ielements_save.add(ielement)
            #if ielement == 73560:
                #print element                
                #aaa

        ielements_save = list(ielements_save)
        ielements_save.sort()

        elements2 = elements[ielements_save]
        regions2 = regions[ielements_save]

        # renumbers mesh
        nelements2, three = elements2.shape
        assert 0 < nelements2 < nelements, 'nelements=%s nelements2=%s'  % (nelements, nelements2)
        
        remap_nodes = False
        if amax(elements2) > len(inodes_save):
            aa
            # build a dictionary of old node ids to new node ids
            nodes_map = {}
            for i in xrange(1, len(inodes_save) + 1):
                nid = inodes_save[i - 1]
                nodes_map[nid] = i

            # update the node ids
            for ielement in xrange(nelements2):
                element = elements2[ielement, :]
                elements[ielement, :] = [nodes_map[nid] for nid in element]

        loads2 = {} # 'Cp', 'Mach', 'U', etc.
        for key, load in loads.iteritems():
            loads2[key] = load[inodes_save]

        self.log.info('---finished make_half_model---')
        return (nodes2, elements2, regions2, loads2)

    def renumber_mesh(self, nodes, elements, regions, loads):
        raise NotImplementedError('disabled...')
        iNodeCounter = 1
        iElementCounter = 1

        NodeOldToNew = {}
        #ElementOldToNew = {}
        nodes2 = {}
        elements2 = {}
        regions2 = {}
        
        nnodes, three = nodes.shape
        for inode in xrange(nnodes):
            node = nodes[inode]
            NodeOldToNew[inode] = iNodeCounter
            nodes2[iNodeCounter] = node
            iNodeCounter += 1

        for (ielement, element) in sorted(elements.iteritems()):
            element2 = []
            for nID in element:
                nID2 = NodeOldToNew[nID]
                element2.append(nID2)
            elements2[iElementCounter] = element2
            regions2[iElementCounter] = regions[ielement]
            iElementCounter += 1
        return (nodes, elements, regions, loads)

    def write_cart3d(self, outfilename, points, elements, regions, loads=None, is_binary=False):
        assert len(points) > 0

        if loads is None or loads == {}:
            loads = {}
            is_loads = False
            print "no loads"
        else:
            is_loads = True

        self.log.info("---writing cart3d file...%r---" % outfilename)
        f = open(outfilename, 'wb')

        float_fmt = '%6.6f'
        int_fmt = self.write_header(f, points, elements, is_loads, is_binary)
        #print "int_fmt =", int_fmt
        self.write_points(f, points, is_binary, float_fmt)
        self.write_elements(f, elements, is_binary, int_fmt)
        self.write_regions(f, regions, is_binary)
        
        if is_loads:
            assert is_binary is False, 'is_binary=%r is not supported for loads' % is_binary
            self.write_loads(f, loads, is_binary, float_fmt)
        f.close()

    def write_header(self, f, points, elements, nloads, is_binary=False):
        npoints, three = points.shape
        nelements, three = elements.shape
        
        if is_binary:
            if nloads == 6:
                msg = pack('>iiii', 3*4, npoints, nelements, 6)
            else:
                msg = pack('>iii', 2*4, npoints, nelements)
            
            int_fmt = None
        else:
            if nloads == 6:
                msg = "%s %s 6\n" % (npoints, nelements)
            else:
                msg = "%s %s\n" % (npoints, nelements)

            # take the max value, string it, and length it
            # so 123,456 is length 6
            int_fmt = '%%%si' % len(str(nelements))
        f.write(msg)
        return int_fmt

    def write_points(self, f, points, is_binary, float_fmt='%6.6f'):
        if is_binary:
            four = pack('>i', 4)
            #f.write(four)

            npoints, three = points.shape
            #points2 = ravel(points)
            #print points2
            floats = pack('>' + 'f'*npoints*3, *ravel(points) )

            f.write(floats)
            f.write(four)
        else:
            savetxt(f, points, float_fmt)

    def write_elements(self, f, elements, is_binary, int_fmt='%6i'):
        if is_binary:
            four = pack('>i', 4)
            f.write(four)
            nelements, three = elements.shape
            ints = pack('>' + 'i'*nelements*3, *ravel(elements) )

            f.write(ints)
            f.write(four)
        else:
            savetxt(f, elements, int_fmt)

    def write_regions(self, f, regions, is_binary):
        if is_binary:
            four = pack('>i', 4)
            f.write(four)

            nregions = len(regions)
            ints = pack('>' + 'i'*nregions, *regions)
            f.write(ints)

            f.write(four)
        else:
            savetxt(f, regions, '%i')

    def write_loads(self, f, loads, is_binary, float_fmt='%6.6f'):
        if is_binary:
            raise NotImplementedError('is_binary=%s' % is_binary)
        else:
            Cp = loads['Cp']
            rho = loads['rho']
            rhoU = loads['rhoU']
            rhoV = loads['rhoV']
            rhoW = loads['rhoV']
            E = loads['E']

            nrows = len(Cp)
            fmt = '%s\n%s%s%s%s%s' % (float_fmt, float_fmt, float_fmt, float_fmt, float_fmt, float_fmt)
            savetxt(f, hstack([Cp, rho, rhoU, rhoV, rhowW, E]), fmt)


class Cart3DAsciiReader(Cart3DReader):
    modelType = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        Cart3DReader.__init__(self, log=log, debug=debug)
        self.isHalfModel = True
        self.cartType = None  # grid, result
        self.nPoints = None
        self.nElements = None
        self.infile = None
        self.infilename = None
        self.readHalf = False
        self.log = get_logger(log, 'debug' if debug else 'info')

    def read_cart3d(self, infilename):
        """extracts the points, elements, and Cp"""
        self.infilename = infilename
        self.log.info("---starting reading cart3d file...%r---" % self.infilename)
        self.infile = open(self.infilename, 'r')
        self.read_header()
        points = self.read_points()
        elements = self.read_elements(bypass=False)
        regions = self.read_regions(bypass=False)
        loads = self.read_results(0, self.infile)

        self.log.info("nPoints=%s  nElements=%s" % (self.nPoints, self.nElements))
        self.log.info("---finished reading cart3d file...%r---" % self.infilename)
        assert self.nPoints > 0, 'nPoints=%s' % self.nPoints
        assert self.nElements > 0, 'nElements=%s' % self.nElements
        return (points, elements, regions, loads)

    def read_header(self):
        line = self.infile.readline()
        sline = line.strip().split()
        if len(sline) == 2:
            self.nPoints, self.nElements = int(sline[0]), int(sline[1])
            self.cartType = 'grid'
        elif len(sline) == 3:
            self.nPoints = int(sline[0])
            self.nElements = int(sline[1])
            self.nResults = int(sline[2])
            self.cartType = 'result'
        else:
            raise ValueError('invalid result type')

        if self.readHalf:
            self.nElementsRead = self.nElements // 2
            self.nElementsSkip = self.nElements // 2
        else:
            self.nElementsRead = self.nElements
            self.nElementsSkip = 0

    def read_points(self):
        """
        A point is defined by x,y,z and the ID is the location in points.
        """
        p = 0
        data = []
        assert self.nPoints > 0, 'nPoints=%s' % self.nPoints
        points = zeros((self.nPoints, 3), 'float64')
        while p < self.nPoints:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                x = data.pop(0)
                y = data.pop(0)
                z = data.pop(0)
                points[p] = [x, y, z]
                p += 1

        maxX = self.get_max(points, 0)
        maxY = self.get_max(points, 1)
        maxZ = self.get_max(points, 2)

        minX = self.get_min(points, 0)
        minY = self.get_min(points, 1)
        minZ = self.get_min(points, 2)

        self.log.info("X  max=%g min=%g" % (maxX, minX))
        self.log.info("Y  max=%g min=%g" % (maxY, minY))
        self.log.info("Z  max=%g min=%g" % (maxZ, minZ))
        return points

    def get_min(self, points, i):
        return amin(points[:, i])

    def get_max(self, points, i):
        return amax(points[:, i])

    def read_elements(self, bypass=False):
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.
        """
        assert bypass == False
        assert self.nElementsRead > 0, 'nPoints=%s' % self.nPoints
        elements = zeros((self.nElementsRead, 3), 'int64')

        self.log.info("nElementsRead=%s nElementsSkip=%s" % (self.nElementsRead, self.nElementsSkip))

        e = 0
        data = []
        while e < self.nElementsRead:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                n1 = int(data.pop(0))
                n2 = int(data.pop(0))
                n3 = int(data.pop(0))
                elements[e] = [n1, n2, n3]
                e += 1

        e = 0
        while e < self.nElementsSkip:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                data.pop()  # popping from the end is faster
                data.pop()
                data.pop()
                e += 1
        return elements

    def read_regions(self, bypass=True):
        regions = zeros(self.nElementsRead, 'int32')
        if bypass:
            for i in xrange(self.nElements):
                self.infile.readline()
        else:
            r = 0
            data = []
            while r < self.nElementsRead:
                data += self.infile.readline().strip().split()
                while len(data) > 0:
                    regions[r] = int(data.pop(0))
                    r += 1

            r = 0
            while r < self.nElementsSkip:
                data += self.infile.readline().strip().split()
                while len(data) > 0:
                    data.pop()
                    r += 1
        return regions

    def read_results(self, i, infile, result_names=None):
        """
        Reads the Cp results.
        Cp = (p - 1/gamma) / (0.5*M_inf*M_inf)
        (pV)^2 = (pu)^2+(pv)^2+(pw)^2
        M^2 = (pV)^2/rho^2

        # ???
        p = (gamma-1)*(e- (rhoU**2+rhoV**2+rhoW**2)/(2.*rho))

        # ???
        rho,rhoU,rhoV,rhoW,rhoE
        
        :param result_names: the results to read; default=None -> All
            result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                            'Mach', 'U', 'V', 'W', 'E']
        """
        loads = {}
        if self.cartType == 'grid':
            return loads

        if result_names is None:
            result_names = ['Cp', 'rho', 'rhoU', 'rhoV', 'rhoW', 'rhoE',
                            'Mach', 'U', 'V', 'W', 'E']
        self.log.info('---starting read_results---')

        results = zeros((self.nPoints, 6), 'float64')

        nResultLines = int(ceil(self.nResults / 5.)) - 1
        for ipoint in xrange(self.nPoints):
            #print "pointNum = ", pointNum
            # rho rhoU,rhoV,rhoW,pressure/rhoE/E
            sline = infile.readline().strip().split()
            i += 1
            for n in xrange(nResultLines):
                sline += infile.readline().strip().split()  # Cp
                i += 1
                #gamma = 1.4
                #else:
                #    p=0.
            sline = self._get_list(sline)
            
            # Cp
            # rho       rhoU      rhoV      rhoW      E
            # 0.416594 
            # 1.095611  0.435676  0.003920  0.011579  0.856058  
            results[ipoint, :] = sline

            #p=0
            #cp = sline[0]
            #rho = float(sline[1])
            #if(rho > abs(0.000001)):
                #rhoU = float(sline[2])
                #rhoV = float(sline[3])
                #rhoW = float(sline[4])
                #rhoE = float(sline[5])
                #mach2 = (rhoU) ** 2 + (rhoV) ** 2 + (rhoW) ** 2 / rho ** 2
                #mach = sqrt(mach2)
                #if mach > 10:
                    #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (pointNum, cp, mach, rho, rhoU, rhoV, rhoW))
            #print "pt=%s i=%s Cp=%s p=%s" %(pointNum,i,sline[0],p)
        del sline
        
        #print "shpae =", results.shape
        Cp = results[:, 0]
        rho = results[:, 1]
        rhoU = results[:, 2]
        rhoV = results[:, 3]
        rhoW = results[:, 4]
        E = results[:, 5]

        ibad = where(rho <= 0.000001)[0]
        if len(ibad) > 0:
        
            if 'Mach' in result_names:
                Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2)# / rho
                Mach[ibad] = 0.0
            if 'U' in result_names:
                U = rhoU / rho
                U[ibad] = 0.0
            if 'U' in result_names:
                V = rhoV / rho
                V[ibad] = 0.0
            if 'W' in result_names:
                W = rhoW / rho
                W[ibad] = 0.0
            #if 'rhoE' in result_names:
                #rhoE = rhoE / rho
                #e[ibad] = 0.0
            
            is_bad = True
            n = 0
            #for i in ibad:
                #print("nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (i, Cp[i], Mach[i], rho[i],
                #                                                                  rhoU[i], rhoV[i], rhoW[i]))
                #Mach[i] = 0.0
                #n += 1
                #if n > 10:
                #    break
        else:
            is_bad = False
        
        loc = locals()
        if 'Cp' in result_names:
            loads['Cp'] = Cp
        if 'rhoU' in result_names:
            loads['rhoU'] = rhoU
        if 'rhoV' in result_names:
            loads['rhoV'] = rhoV
        if 'rhoW' in result_names:
            loads['rhoW'] = rhoW
        #if 'rhoE' in result_names:
            #loads['rhoE'] = rhoE

        if 'rho' in result_names:
            loads['rho'] = rho

        if 'Mach' in result_names:
            if not is_bad:
                #Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2) / rho
                Mach = sqrt(rhoU**2 + rhoV**2 + rhoW**2)
                
            loads['Mach'] = Mach
        if 'U' in result_names:
            if not is_bad:
                U = rhoU / rho
            loads['U'] = U
        if 'V' in result_names:
            if not is_bad:
                V = rhoV / rho
            loads['V'] = V
        if 'W' in result_names:
            if not is_bad:
                W = rhoW / rho
            loads['W'] = W
        if 'E' in result_names:
            #if not is_bad:
                #E = rhoE / rho
            loads['E'] = E
        
        #i = where(Mach == max(Mach))[0][0]
        #self.log.info("i=%s Cp=%s rho=%s rhoU=%s rhoV=%s rhoW=%s Mach=%s" % (i, Cp[i], rho[i], rhoU[i],
        #              rhoV[i], rhoW[i], Mach[i]))
        self.log.info('---finished read_results---')
        return loads

    def _get_list(self, sline):
        """Takes a list of strings and converts them to floats."""
        try:
            sline2 = convert_to_float(sline)
        except ValueError:
            print("sline = %s" % sline)
            raise SyntaxError('cannot parse %s' % sline)
        return sline2

    def export_to_nastran(self, fname, points, elements, regions):
        f = open(fname, 'wb')
        nnodes, three = points.shape
        for nid in xrange(nnodes):
            (x, y, z) = nodes[nid, :]
            f.write(print_card(['GRID', nid + 1, '', x, y, z]))

        e = 1e7
        nu = 0.3
        g = ''
        thickness = 0.25
        setRegions = list(set(regions))
        for pidMid in setRegions:
            f.write(print_card(['MAT1', pidMid, e, g, nu]))
            f.write(print_card(['PSHELL', pidMid, pidMid, thickness]))

        for eid, (nodes, region) in enumerate(1, izip(elements, regions)):
            (n1, n2, n3) = nodes
            f.write(print_card(['CTRIA3', eid, region, n1, n2, n3]))
        f.close()

    def get_segments(self, nodes, elements):
        segments = {}  # key=eid,
        lengths = {}
        for eid, e in elements.iteritems():
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))
            print(eid, a, b, c)

            #a.sort()
            #b.sort()
            #c.sort()
            if a in segments:
                segments[a].append(eid)
            else:
                segments[a] = [eid]
                #print(a)
                #print nodes[1]
                lengths[a] = nodes[a[1]] - nodes[a[0]]

            if b in segments:
                segments[b].append(eid)
            else:
                segments[b] = [eid]
                lengths[b] = nodes[b[1]] - nodes[b[0]]

            if c in segments:
                segments[c].append(eid)
            else:
                segments[c] = [eid]
                lengths[c] = nodes[c[1]] - nodes[c[0]]

        return segments, lengths

    def make_quads(self, nodes, elements):
        segments, lengths = self.get_segments(nodes, elements)
        for eid, e in elements.iteritems():
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))
            #a.sort()
            #b.sort()
            #c.sort()
            print(eid, e)
            print(segments[a])
            print(lengths[a])
            print(len(segments[a]))

            eidA = self.get_segment(a, eid, segments)
            eidB = self.get_segment(b, eid, segments)
            eidC = self.get_segment(c, eid, segments)
            print("eidA=%s eidB=%s eidC=%s" % (eidA, eidB, eidC))
            if eidA:
                i = 0
                e2 = elements[eidA]
                self.check_quad(nodes, eid, eidA, e, e2, a, b, c, i)
                del segments[a]
            if eidB:
                i = 1
                e2 = elements[eidB]
                self.check_quad(nodes, eid, eidB, e, e2, a, b, c, i)
                del segments[b]
            if eidC:
                i = 2
                e2 = elements[eidC]
                self.check_quad(nodes, eid, eidC, e, e2, a, b, c, i)
                del segments[c]

            print("------")
            #break
        #for segment in segments:
        asdf

    def check_quad(self, nodes, eid, eidA, e, e2, a, b, c, i):
        r"""
        ::
          A----B
          | \ e|
          |e2 \|
          C----D

        two tests
           1.  folding angle A-B x A-C
           2a. abs(A-C) - abs(B-D)  = 0  (abs to prevent 2L)
           2b. abs(A-B) - abs(C-D)  = 0
        """

        iplus1 = i + 1
        iplus2 = i + 2

        if iplus1 > 2:
            iplus1 -= 3
        if iplus2 > 2:
            iplus2 -= 3
        print(i, iplus1)
        print(iplus1, iplus2)
        print(iplus2, i)
        AD = nodes[e[i]] - nodes[e[iplus1]]
        AB = nodes[e[iplus1]] - nodes[e[iplus2]]
        BD = nodes[e[iplus2]] - nodes[e[i]]

        print(AD)
        print(AB)
        print(BD)
        print(e2)
        j = e2.index(e[i])

        jplus1 = j + 1
        jplus2 = j + 2
        if jplus1 > 2:
            jplus1 -= 3
        if jplus2 > 2:
            jplus2 -= 3

        print("DA = ", e[j], e[jplus1])
        DA = nodes[e[j]] - nodes[e[jplus1]]
        print(DA)

        asdf

    def get_segment(self, a, eid, segments):
        if a in segments:
            aElems = segments[a]
            print(aElems)
            i = aElems.index(eid)
            #print i
            aElems.pop(i)
            #print aElems
            eidA = aElems[0]
            #eidA = elements[a]
            print("eidA = ", eidA)
            return eidA
        return None


class Cart3DBinaryReader(FortranFile, Cart3DAsciiReader):
    modelType = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        Cart3DReader.__init__(self)
        FortranFile.__init__(self)

        self.isHalfModel = True
        self.cartType = None  # grid, result
        self.nPoints = None
        self.nElements = None
        self.infile = None
        self.infilename = None
        self.readHalf = False
        self.isNodeArray = True

        self.make_op2_debug = False
        self.n = 0
        self.log = get_logger(log, 'debug' if debug else 'info')

    def read_cart3d(self, infileName):
        self.infilename = infileName
        self.op2 = open(infileName, 'rb')
        (self.nPoints, self.nElements) = self.read_header()
        points = self.read_points(self.nPoints)
        elements = self.read_elements(self.nElements)
        regions = self.read_regions(self.nElements)

        loads = {}
        return (points, elements, regions, loads)

    def read_header(self):
        data = self.op2.read(4)
        size, = unpack('>i', data)
        print "size = ",size

        data = self.op2.read(size)
        so4 = size // 4  # size over 4
        if so4 == 3:
            (nPoints, nElements, nResults) = unpack('>iii', data)
            self.log.info("nPoints=%s nElements=%s nResults=%s" % (nPoints, nElements, nResults))
        elif so4 == 2:
            (nPoints, nElements) = unpack('>ii', data)
            self.log.info("nPoints=%s nElements=%s" % (nPoints, nElements))
        else:
            raise RuntimeError('in the wrong spot...endian...')
        self.op2.read(8)  # end of first block, start of second block

        return (nPoints, nElements)

    def read_points(self, npoints):
        #print "starting Nodes"
        #isBuffered = True
        size = npoints * 12  # 12=3*4 all the points
        Format = '>' + 'f' * 3000  # 3000 floats; 1000 points

        #nodes = {}
        np = 0
        points = zeros((npoints, 3), 'float64')
        while size > 12000:  # 12k = 4 bytes/float*3 floats/point*1000 points
            #print "size = ",size
            n = np + 999
            #print "nStart=%s np=%s" % (n ,np)
            data = self.op2.read(4 * 3000)
            nodeXYZs = list(unpack(Format, data))
            while nodeXYZs:
                z = nodeXYZs.pop()
                y = nodeXYZs.pop()
                x = nodeXYZs.pop()
                print x, y, z
                #assert n not in nodes, 'nid=%s in nodes' % n
                points[n, :] = [x, y, z]
                n -= 1
                np += 1
            size -= 4 * 3000

        assert size >= 0
        #print "size = ",size
        #while size>0: # 4k is 1000 points
        n = npoints-1

        if size > 0:
            #leftover = size-(nPoints-np)*12
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover, size//4)
            Format = '>' + 'f' * (size // 4)

            #print "len(data) = ",len(data)
            try:
                nodeXYZs = list(unpack(Format, data))
            except:
                print "len(data)=%s" % len(data)
                raise

            while nodeXYZs:
                z = nodeXYZs.pop()
                y = nodeXYZs.pop()
                x = nodeXYZs.pop()
                #assert n not in nodes, 'nid=%s in nodes' % n
                points[n, :] = [x, y, z]
                n -= 1
                np += 1
            size = 0

        #for p,point in sorted(nodes.iteritems()):
            #print("%s %s %s" % tuple(point))

        #if isBuffered:
            #pass
        #else:
            #raise RuntimeError('unBuffered')

        #for nid in xrange(nPoints):
            #assert nid in points, 'nid=%s not in points' % nid
        self.op2.read(8)  # end of second block, start of third block
        print points
        return points

    def read_elements(self, nElements):
        self.nElementsRead = nElements
        self.nElementsSkip = 0
        #print "starting Elements"
        #isBuffered = True
        size = nElements * 12  # 12=3*4 all the elements
        Format = '>' + 'i' * 3000

        elements = zeros((self.nElements, 3), 'int32')

        ne = 0
        while size > 12000:  # 4k is 1000 elements
            #print "size = ",size
            e = ne + 1000
            data = self.op2.read(4 * 3000)
            nodes = list(unpack(Format, data))
            while nodes:
                n3 = nodes.pop()
                n2 = nodes.pop()
                n1 = nodes.pop()
                element = [n1, n2, n3]
                elements[e, :] = element
                e -= 1
                ne += 1
            size -= 4 * 3000

        assert size >= 0
        #print "size = ",size
        #while size>0: # 4k is 1000 elements
        if size > 0:
            #leftover = size-(nElements-ne)*12
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover,size//4)
            Format = '>' + 'i' * (size // 4)

            #print "len(data) = ",len(data)
            try:
                nodes = list(unpack(Format, data))
            except:
                print "len(data)=%s" % len(data)
                raise

            e = nElements-1
            while nodes:
                n3 = nodes.pop()
                n2 = nodes.pop()
                n1 = nodes.pop()
                print elements
                elements[e, :] = [n1, n2, n3]
                e -= 1
                ne += 1
            size = 0

        #for p,point in sorted(nodes.iteritems()):
            #print "%s %s %s" %(tuple(point))

        #if isBuffered:
            #pass
        #else:
            #raise RuntimeError('unBuffered')

        self.op2.read(8)  # end of third (element) block, start of regions (fourth) block
        #print "finished Elements"
        return elements

    def read_regions(self, nElements):
        #print "starting Regions"
        #isBuffered = True
        size = nElements * 4  # 12=3*4 all the elements
        Format = '>' + 'i' * 3000

        regions = {}

        nr = 0
        while size > 12000:  # 12k is 3000 elements
            #print "size = ",size
            data = self.op2.read(4 * 3000)
            regionData = list(unpack(Format, data))

            r = nr + 3000
            #print "len(regionData) = ",len(regionData)
            while regionData:
                regions[r] = regionData.pop()
                r -= 1
                nr += 1
            size -= 4 * 3000
            #print "size = ",size

        assert size >= 0
        #print "size = ",size

        if size > 0:
            #leftover = size-(nElements-nr)*4
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover,size//4)
            Format = '>' + 'i' * (size // 4)

            #print "len(data) = ",len(data)

            try:
                regionData = list(unpack(Format, data))
            except:
                print "len(data)=%s" % len(data)
                raise
            r = nElements
            while regionData:
                regions[r] = regionData.pop()
                r -= 1
                nr += 1
            size = 0

        self.op2.read(4)  # end of regions (fourth) block
        #print "finished Regions"
        return regions

    def read_results(self, i, infile, result_names=None):
        pass


def generic_cart3d_reader(infileName, log=None, debug=False):
    #print("infileName = ", infileName)
    if is_binary(infileName):
        obj = Cart3DBinaryReader(log, debug)
    else:
        obj = Cart3DAsciiReader(log, debug)
    return obj


if __name__ == '__main__':
    import time
    t0 = time.time()
    
    # binary
    #cart3dGeom = os.path.join('flat.tri')  # half
    #full_model = os.path.join('flat_full.tri')
    
    cart3dGeom = 'Cart3d_bwb.i.tri'
    rewrite    = 'Cart3d_bwb_rewrite.tri'
    half_model = 'Cart3d_bwb_half.tri'
    full_model = 'Cart3d_bwb_full.tri'
    
    cart = generic_cart3d_reader(cart3dGeom)  # ascii/binary

    (points, elements, regions, loads) = cart.read_cart3d(cart3dGeom)
    cart.write_cart3d(rewrite, points, elements, regions)

    (points, elements, regions, loads) = cart.make_half_model(points, elements, regions, loads)
    cart.write_cart3d(half_model, points, elements, regions)

    (points2, elements2, regions2, loads2) = cart.make_mirror_model(points, elements, regions, loads)
    cart.write_cart3d(full_model, points2, elements2, regions2)
    
    print "dt = ", time.time() - t0
    sys.exit('made half model')

    # ascii
    cart3dGeom = os.path.join('models', 'threePlugs.a.tri')
    outfilename = os.path.join('models', 'threePlugs2.a.tri')
    cart2 = generic_cart3d_reader(cart3dGeom)
    (points, elements, regions, loads) = cart2.read_cart3d(cart3dGeom)
    cart2.write_cart3d(outfilename, points, elements, regions)
    #print points


if 0:
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath = os.path.join(basepath, 'outputsFinal')
    #bdfModel   = os.path.join(configpath,'aeroModel.bdf')
    #assert os.path.exists(bdfModel),'%r doesnt exist' % bdfModel
    os.chdir(workpath)
    print("basepath", basepath)

    cart3dGeom = os.path.join(configpath, 'Cart3d_bwb2.i.tri')
    cart3dGeom2 = os.path.join(workpath, 'Cart3d_half.i.tri')
    cart3dGeom3 = os.path.join(workpath, 'Cart3d_full.i.tri')
    #cart3dGeom4 = os.path.join(workpath,'Cart3d_full2.i.tri')

    cart = Cart3DAsciiReader()
    (points, elements, regions, loads) = cart.read_cart3d(cart3dGeom)
    (points, elements, regions, loads) = cart.make_half_model(points, elements, regions, loads, axis='y', remap_nodes=True)
    (points, elements, regions, loads) = cart.make_mirror_model(points, elements, regions, loads, axis='y')
    cart.write_cart3d(cart3dGeom2, points, elements, regions)

    #cart = Cart3DAsciiReader(cart3dGeom)  # bJet.a.tri
    #(cartPoints,elements,regions,loads) = cart.read()
    #points2 = {}
    #for (iPoint,point) in sorted(points.iteritems()):
        #(x,y,z) = point
        #print "p=%s x=%s y=%s z=%s  z2=%s" %(iPoint,x,y,z,z+x/10.)
        #points2[iPoint] = [x,y,z+x/10.]
    #(points, elements, regions, loads) = cart.make_mirror_model(points2, elements, regions, loads)

    cart2 = Cart3DAsciiReader(cart3dGeom2)
    (points, elements, regions, loads) = cart2.read_cart3d()

    #makeFullModel
    (points, elements, regions, loads) = cart2.make_mirror_model(points, elements, regions, loads)  # half model values going in
    cart2.write_cart3d(cart3dGeom3, points, elements, regions)

    #cart3 = Cart3DAsciiReader(cart3dGeom2)
    #(points, elements, regions, loads) = cart3.read_cart3d()
    #(points, elements, regions, loads) = cart3.make_mirror_model(points, elements, regions, loads)


    #print "loads = ",list_print(loads), len(loads)
    #cartOutfile = os.path.join(workpath, 'bJet.a.tri_test')
    #cart.writeInfile(cartOutfile, cartPoints, elements, regions)
