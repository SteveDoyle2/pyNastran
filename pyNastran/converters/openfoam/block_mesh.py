"""
defines:
 - model = read_block_mesh(block_mesh_dict_filename, log=None, debug=False)
"""
from __future__ import print_function
import os
from copy import deepcopy
from codecs import open
from collections import defaultdict, OrderedDict
from itertools import count
from six import iteritems
from six.moves import range, zip

import numpy as np
from numpy import array, cross, unique, allclose, zeros, arange, ravel, ones, argsort
from numpy.linalg import norm  # type: ignore

from pyNastran.converters.openfoam.openfoam_parser import FoamFile, convert_to_dict #, write_dict
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.log import get_logger2
from pyNastran.utils import print_bad_path

def read_block_mesh(block_mesh_dict_filename, log=None, debug=False):
    """functional interface to BlockMesh"""
    model = BlockMesh(log=log, debug=debug)
    model.read_openfoam(block_mesh_dict_filename)
    return model

def area_centroid(n1, n2, n3, n4):
    """calculates the area and centroid of a quad"""
    centroid = (n1 + n2 + n3 + n4) / 4.
    n = cross(n3 - n1, n4 - n2)
    area = 0.5 * norm(n)
    return area, centroid

def volume8(n1, n2, n3, n4, n5, n6, n7, n8):
    """calculates the volume of a hex"""
    (A1, c1) = area_centroid(n1, n2, n3, n4)
    (A2, c2) = area_centroid(n5, n6, n7, n8)
    V = (A1 + A2) / 2. * norm(c1 - c2)
    return V

def get_not_indexes(a, indices):
    """only works for 1D"""
    ia = np.indices(a.shape)
    not_indices = np.setxor1d(ia, indices)
    return not_indices


class FaceFile(object):
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_face_file(self, face_filename, ifaces_to_read=None):
        #p = FoamFile(face_filename)
        #lines = p.read_foam_file()
        #print('converting')
        with open(face_filename, 'r') as face_file:
            self._read_face_file(face_file, ifaces_to_read=ifaces_to_read)

    def _read_face_file(self, face_file, ifaces_to_read=None):
        i = 0
        nfaces = 0
        while nfaces == 0:
            line = face_file.readline()
            i += 1
            try:
                nfaces = int(line)
            except ValueError:
                pass
        line = face_file.readline()
        i += 1

        #print('nfaces = %s' % nfaces)
        #print('lineA = %r' % line)

        #print('building faces')
        assert nfaces > 0, nfaces

        if ifaces_to_read is None:
            faces = ones((nfaces, 4), dtype='int32') * -1

            for j in range(nfaces):
                # 3(a b c) to [a, b, c]
                # 4(a b c d) to [a, b, c, d]
                face_line = face_file.readline()
                i += 1
                sline = face_line[1:].strip('( )\n\r').split()
                try:
                    if len(sline) == 3:
                        faces[j, :3] = sline[:3]
                    elif len(sline) == 4:
                        faces[j, :] = sline[:4]  # TODO: quads only
                    else:
                        msg = 'The sline is the wrong length (3/4 required)\n'
                        msg += 'sline = %s\n' % str(sline)
                        msg += 'line = %r\n' % face_line.rstrip()
                        raise RuntimeError(msg)
                except:
                    print('face = %r' % face_line)
                    print(sline, i, len(sline))
                    raise
        else:
            # we need to sort our faces as we read the faces sequentially
            # however the output needs to be in unsorted order so we save the sorting
            # order and back apply it at the end
            isort = argsort(ifaces_to_read)
            ifaces_to_read_sorted = ifaces_to_read[isort]

            nfaces2 = len(ifaces_to_read)  # must be sorted
            assert len(isort) == nfaces2
            faces = ones((nfaces2, 4), dtype='int32') * -1
            ni = 0
            for j in range(nfaces):
                face_line = face_file.readline()
                i += 1
                try:
                    ni_face = ifaces_to_read_sorted[ni]
                except:
                    raise
                #if j % 5000 == 0 or j == ni_face:
                    #print('j=%s ni=%s' % (j, ni))
                if j == ni_face:
                    # 3(a b c) to [a, b, c]
                    # 4(a b c d) to [a, b, c, d]
                    sline = face_line[1:].strip('( )\n\r').split()
                    try:
                        if len(sline) == 3:
                            faces[ni, :3] = sline[:3]
                        elif len(sline) == 4:
                            faces[ni, :] = sline[:4]
                        else:
                            msg = 'The sline is the wrong length (3/4 required)\n'
                            msg += 'sline = %s\n' % str(sline)
                            msg += 'line = %r\n' % face_line
                            raise RuntimeError(msg)
                    except:
                        print('face = %r' % face_line.rstrip())
                        print(sline, i, j, ni)
                        raise
                    ni += 1
            #assert faces[:, 0].min() > 0, 'where-1 = %s' % where(faces[:, 0] == -1)[0]
            self.log.info(faces)
            faces[isort, :] = deepcopy(faces[:, :])
        self.log.info('faces.shape = %s' % str(faces.shape))
        return faces


class PointFile(object):
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_point_file(self, point_filename, ipoints_to_read=None):
        #p = FoamFile(face_filename)
        #lines = p.read_foam_file()
        #self.log.info('converting')
        i = 0
        npoints = 0
        with open(point_filename, 'r') as points_file:
            while npoints == 0:
                line = points_file.readline()
                i += 1
                try:
                    npoints = int(line)
                except ValueError:
                    pass
            line = points_file.readline()
            i += 1

            self.log.info('npoints = %s' % npoints)
            #print('lineA = %r' % line)

            self.log.info('building points')
            assert npoints > 0, npoints
            if ipoints_to_read is not None:
                ipoints_to_read.sort()
                npoints2 = len(ipoints_to_read)
                points = zeros((npoints2, 3), dtype='float32')

                self.log.info('npoints2 = %s' % npoints2)
                ni = 0
                for j in range(npoints):
                    try:
                        ni_point = ipoints_to_read[ni]
                    except:
                        print('j=%s ni=%s' % (j, ni))
                        raise
                    if j == ni_point:
                        point = points_file.readline()
                        i += 1
                        sline = point.strip('( )\n\r').split()
                        try:
                            points[ni, :] = sline
                        except:
                            print('point = %r' % point)
                            print(sline, i, j, ni)
                            raise
                        ni += 1
            else:
                points = zeros((npoints, 3), dtype='float32')
                for j in range(npoints):
                    point = points_file.readline()
                    i += 1
                    sline = point.strip('( )\n\r').split()
                    try:
                        points[j, :] = sline
                    except:
                        print('point = %r' % point)
                        print(sline, i)
                        raise
        self.log.info('points.shape = %s' % str(points.shape))
        return points


class BoundaryFile(object):
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_boundary_file(self, boundary_filename):
        #p = FoamFile(face_filename)
        #lines = p.read_foam_file()
        #self.log.info('converting')
        i = 0
        nboundaries = 0

        with open(boundary_filename, 'r') as boundary_file:
            while nboundaries == 0:
                line = boundary_file.readline()
                i += 1
                try:
                    nboundaries = int(line)
                except ValueError:
                    pass
            line = boundary_file.readline()
            i += 1

            self.log.info('nboundaries = %s' % nboundaries)
            #self.log.info('lineA = %r' % line)

            self.log.info('building boundaries')
            boundaries = OrderedDict()
            boundaries = self._read_boundaries(boundary_file, i, nboundaries, boundaries)
        return boundaries

    def _read_boundaries(self, boundary_file, i, nboundaries, boundaries, basename=''):
        assert nboundaries > 0, nboundaries
        #read_next = 0
        for unused_j in range(nboundaries):
            # 3(a b c) to [a, b, c]
            # 4(a b c d) to [a, b, c, d]
            nameline = boundary_file.readline()
            i += 1
            name = nameline.strip()
            self.log.info('name = %r' % (basename + name))

            unused_openline = boundary_file.readline()
            i += 1
            typeline = boundary_file.readline()
            i += 1
            word, boundary_type = typeline.strip('\n\r\t ;').split()

            nfacesline = boundary_file.readline()
            i += 1
            sline = nfacesline.strip('\n\r\t ;').split()
            word = sline[0]
            if word == 'inGroups' and len(sline) == 1:
                #nboundary_line = boundary_file.readline(); i+=1
                #nboundaries2 = int(nboundary_line)
                #openline = boundary_file.readline(); i+=1
                for unused_ii in range(7):
                    groupline = boundary_file.readline()
                    i += 1
                    #self.log.info(ii, groupline)
                sline = groupline.strip('\n\r\t ;').split()

                word, nfaces = sline
                nfaces = int(nfaces)
                self.log.info('nfaces = %r' % nfaces)

                startfacesline = boundary_file.readline()
                i += 1
                self.log.info('startfacesline = %r' % startfacesline)
                word, startfaces = startfacesline.strip('\n\r\t ;').split()
                startfaces = int(startfaces)
                unused_closeline = boundary_file.readline()
                i += 1

                boundary_name = basename + name
                if boundary_name in boundaries:
                    msg = ('boundary_name=%r is already defined...'
                           'boundaries must have unique names' % boundary_name)
                    raise KeyError(msg)
                boundaries[boundary_name] = [boundary_type, nfaces, startfaces]

                #self._read_boundaries(boundary_file, i, nboundaries2, boundaries, basename + name)
            else:
                if word == 'inGroups' and len(sline) == 2:
                    unused_ingroups = nfaces
                    nfacesline = boundary_file.readline()
                    i += 1
                    word, nfaces = nfacesline.strip('\n\r\t ;').split()
                else:
                    word, nfaces = sline

                nfaces = int(nfaces)
                self.log.info('nfaces = %r' % (nfaces))

                startfacesline = boundary_file.readline()
                i += 1
                self.log.info('startfacesline = %r' % startfacesline)
                word, startfaces = startfacesline.strip('\n\r\t ;').split()
                startfaces = int(startfaces)

                unused_closeline = boundary_file.readline()
                i += 1

                boundary_name = basename + name
                if boundary_name in boundaries:
                    raise KeyError('boundary_name=%r is already defined...'
                                   'boundaries must have unique names' % boundary_name)
                boundaries[boundary_name] = [boundary_type, nfaces, startfaces]

        for name, boundary in iteritems(boundaries):
            self.log.info('name=%s boundary=%s' % (name, boundary))
        return boundaries


class Boundary(object):
    """defines Boundary"""
    def __init__(self, log=None, debug=False):
        """creates a Boundary object"""
        self.debug = False
        #log = None
        self.log = get_logger2(log, debug=debug)

    def read_openfoam(self, point_filename, face_filename, boundary_filename):
        """reads a Boundary file"""
        assert os.path.exists(face_filename), print_bad_path(face_filename)
        assert os.path.exists(point_filename), print_bad_path(point_filename)
        assert os.path.exists(boundary_filename), print_bad_path(boundary_filename)

        #self.log.info('face_filename = %r' % face_filename)
        #self.log.info('point_filename = %r' % point_filename)
        #self.log.info('boundary_filename = %r' % boundary_filename)

        assert 'faces' in face_filename, face_filename
        assert 'points' in point_filename, point_filename
        assert 'boundary' in boundary_filename, boundary_filename

        #print('starting Boundary')
        point_file = PointFile(log=self.log, debug=self.debug)
        #from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
        #self.log.info(dir(f))

        face_file = FaceFile(log=self.log, debug=self.debug)

        boundary_file = BoundaryFile(log=self.log, debug=False)
        boundaries = boundary_file.read_boundary_file(boundary_filename)

        #if 0:
            #foam = FoamFile(boundary_filename, log=p.log)

            #print('getting lines')
            #blines = foam.read_foam_file()
            #print('converting')
            #bd = convert_to_dict(foam, blines, debug=True)
            #del blines


        self.log.info('getting npoints')
        #pself.log.info(write_dict(d))

        #-------------------------------------------
        # count number of faces by looking at the boundary info
        # so we can allocate faces2
        nfaces2 = 0
        ifaces_to_read = []
        #f_boundary_faces = open('boundary_faces.py', 'wb')
        for name, boundary in iteritems(boundaries):
            # type            patch;  # 0
            # nFaces          nFaces; # 1
            # startFace       777700; # 2
            self.log.info('boundary[%s] = %s' % (name, boundary))
            nfacesi = boundary[1]
            startface = int(boundary[2])
            nfaces2 += nfacesi
            new_faces = list(arange(nfacesi, dtype='int32') + startface)
            #f_boundary_faces.write('boundary_faces[%s, %s] = %s\n' % (
                #name, len(new_faces), new_faces))
            ifaces_to_read += new_faces

        self.log.info('nfaces2 = %s' % nfaces2)
        ifaces_to_read = ravel(ifaces_to_read)
        if len(ifaces_to_read) != nfaces2:
            raise RuntimeError('len(ifaces_to_read)=%s nfaces2=%s' % (
                ifaces_to_read.shape, nfaces2))
        self.log.info(ifaces_to_read)

        faces = face_file.read_face_file(face_filename, ifaces_to_read=ifaces_to_read)
        #faces = f.read_face_file(face_filename, ifaces_to_read=None)
        del ifaces_to_read

        if 0:  # pragma: no cover
            # doesn't work for some reason...
            # we want to only plot a subset of faces to reduce the data set
            # that works, but we also need to decrease the number of nodes
            # (they take wayyy too long)

            # so we take our faces, get the unique nodes
            # sort them so they're consistent with the order in the file
            # using the same block of code that works in the face reader,
            #but it still fails for some reason...

            # after this step, we renumber the faces with the adjusted node ids
            ipoints_to_read = unique(faces.ravel())
            self.log.info('nnodes = %s' % len(ipoints_to_read))
            ipoints_to_read.sort()
            self.log.info('ipoints_to_read = %s' % ipoints_to_read)
        else:
            ipoints_to_read = None
        nodes = point_file.read_point_file(point_filename, ipoints_to_read=ipoints_to_read)

        if ipoints_to_read is not None:
            nid_to_ipoint = {}
            for inid, nid in enumerate(ipoints_to_read):
                nid_to_ipoint[nid] = inid

            self.log.info('%s %s' % (faces, faces.max()))
            for iface, unused_face in enumerate(faces):
                #print('face      = %s' % face)
                faces[iface, 0] = nid_to_ipoint[faces[iface, 0]]
                faces[iface, 1] = nid_to_ipoint[faces[iface, 1]]
                faces[iface, 2] = nid_to_ipoint[faces[iface, 2]]
                #print('faces[%i] = %s' % (i, faces[i, :]))
            self.log.info('%s %s' % (faces, faces.max()))
            self.log.info('done...')
            del ipoints_to_read
            del nid_to_ipoint
        #-------------------------------------------
        # keep only the required faces
        iface = 0
        #faces2 = zeros((nfaces2, 4), dtype='int32')
        names = zeros(nfaces2, dtype='int32')
        iname = 1
        snames = [None] * (len(boundaries) + 1)
        self.log.info('')
        for name, boundary in iteritems(boundaries):
            self.log.info('iname=%s name=%s boundary=%s' % (iname, name, boundary))
            # type            patch;
            # nFaces          nFaces;
            # startFace       777700;
            try:
                unused_type = boundary[0]
                nfacesi = int(boundary[1])
                startface = int(boundary[2])
            except:
                print(boundary.keys())
                raise
            #faces2[iface:iface+nfacesi] = faces[startface:startface + nfacesi]
            names[iface:iface+nfacesi] = iname
            snames[iname] = name
            iface += nfacesi
            iname += 1
        #del faces
        quads = faces


        #if 0:
            #f_boundary_faces.write('\n\n---Faces----\n')
            #for iface, face in enumerate(faces):
                #pid = names[iface]
                #name = snames[pid]
                #f_boundary_faces.write('%i (%i %i %i %i) pid=%s name=%s\n' % (
                    #iface, face[0], face[1], face[2], face[3], pid, name))

            #f_boundary_faces.write('\n\n---First Faces----\n')
            #pid_save = set([])
            #for iface, face in enumerate(faces):
                #pid = names[iface]
                #if pid not in pid_save:
                    #name = snames[pid]
                    #f_boundary_faces.write('%i (%i %i %i %i) pid=%s name=%s\n' % (
                        #iface, face[0], face[1], face[2], face[3], pid, name))
                    #pid_save.add(pid)

        # only save the unique nodes
        # ...
        #unodes = unique(quads.ravel())
        #unodes.sort()
        #nodes = nodes[unodes, :]

        # renumber the nodes on the faces
        # ...
        self.log.debug('names=%s; max=%s min=%s' % (names, names.max(), names.min()))

        print('done with Boundary')
        #self.nodes = nodes
        return nodes, quads, names


class BlockMesh(object):
    """defines BlockMesh"""
    def __init__(self, log=None, debug=False):
        """creates BlockMesh"""
        debug = False
        #log = None
        self.log = get_logger2(log, debug=debug)

        # arrays
        self.nodes = None
        self.quads = None
        self.hexas = None
        self.npoints = None
        self.grading = None

        self.iname_to_quads = None
        self.inames = None
        self.bcs = None
        self.iname_to_name = None
        self.iname_to_type = None

    def make_hex_bar(self, unused_bias, ncells):
        """
        bias = Llast/Lfirst
        """
        #k = bias = ** 1/N  # close to this-ish
        k = 1.
        ipoints = arange(ncells + 1) # ipoint
        # xmax = d0 * k**n
        kn = k**ipoints
        L0 = 1.0 / sum(kn)
        x = kn * L0
        return ipoints, x

    def make_hex_bars(self, unused_hex_id):
        #self.hexas = hexas
        #self.npoints = npoints
        #self.grading = grading
        points = []
        line_pairs = [
            # a set is 2 values
            #set1  set2 idir,face
            (0, 1, 4, 5, 0, 2),
            (3, 2, 7, 6, 0, 2),
            (0, 3, 1, 2, 1, 1),
            (4, 7, 5, 6, 1, 1),
            (0, 4, 1, 5, 2, 0),
            (3, 7, 2, 6, 2, 0),
        ]

        points = []
        for grading, unused_hexa in zip(self.grading, self.hexas):
            # x = 0-1 to 4-5
            #     3-2 to 7-6
            # y = 0-3 to 1-2
            #     4-7 to 5-6
            # z = 0-4 to 1-5
            #     3-7 to 2-6
            nx, ny, nz, unused_method = grading
            bias = (nx, ny, nz)
            for line_pair in line_pairs:
                i1, i2, i3, i4, idir, iface = line_pair
                n1 = self.nodes[i1, :]
                n2 = self.nodes[i2, :]
                n3 = self.nodes[i3, :]
                n4 = self.nodes[i4, :]
                unused_ncells_x = grading[idir]
                ncells_y = grading[iface]
                bias_x = bias[idir]
                bias_y = bias[iface]

                unused_npx, x = self.make_hex_bar(bias_x, ncells_y)
                unused_npy, y = self.make_hex_bar(bias_y, ncells_y)

                da = n1 - n2
                db = n3 - n4
                La = norm(da)
                Lb = norm(db)
                xa = La * x
                xb = Lb * x

                for i in range(1, ncells):
                    p1 = n2 + xa[i]
                    p2 = n4 + xb[i]
                    dp = p2 - p1
                    L = norm(dp)
                    unused_yout = L * y

                    unused_new_points = p1 + dp * y
                    points.append(1)
                    i / L
                    #p1 =

    def read_openfoam(self, block_mesh_name='blockMeshDict'):
        """reads the BlockMesh file"""
        self.log.info('block_mesh_name = %r' % block_mesh_name)
        foam_file = FoamFile(block_mesh_name)
        foam_lines = foam_file.read_foam_file()

        foam_file_dict = convert_to_dict(self, foam_lines, debug=True)
        #print(write_dict(foam_file_dict))
        unused_keys = foam_file_dict.keys()

        vertices = foam_file_dict['vertices']
        blocks = foam_file_dict['blocks']
        boundaries = foam_file_dict['boundary']
        #self.log.info(boundaries)

        nodes = []
        for unused_ivertex, vertex in iteritems(vertices):
            x, y, z = vertex.strip('() ').split()
            nodes.append([x, y, z])

        hexas = []
        npoints = []
        grading = []
        for key, block in iteritems(blocks):
            hexa, npointsi, gradingi, blank = block.split(')')
            assert blank == '', '%r' % blank

            hexa = hexa.replace('hex (', '').split()
            npointsi = npointsi.strip('( ').split()
            gradingi = gradingi.split('(')[1].split()
            #self.log.info('%s %s' % (npointsi, gradingi))
            #self.log.info(hexa)
            hexa = [int(i) for i in hexa]
            #print('hexa', key, hexa)
            hexas.append(hexa)
            npoints.append(npointsi)
            grading.append(gradingi)

        bc_map = {
            'wall' : 0,
            'patch' : 1,
            'symmetry' : 2,
        }
        iname_to_quads = defaultdict(list)
        inames = []
        bcs = []
        iname = 1
        iname_to_name = {}
        iname_to_type = {}
        quads = []
        for key, boundary in iteritems(boundaries):
            #print key, boundary.keys()
            Type = boundary['type']
            bc = bc_map[Type]
            #bcs.append(bc)
            faces = boundary['faces']
            iname_to_name[iname] = key
            iname_to_type[iname] = Type
            for unused_iface, face in iteritems(faces):
                quad = face.strip('() ').split()
                quad = [int(i) for i in quad]
                quads.append(quad)
                iname_to_quads[iname].append(quad)
                inames.append(iname)
                bcs.append(bc)
            self.log.info('iname=%s -> %s; bc=%s -> %s' % (iname, key, bc, Type))
            iname += 1

        nodes = array(nodes, dtype='float32')
        self.nodes = nodes

        hexas = array(hexas, dtype='int32')
        npoints = array(npoints, dtype='int32')
        grading = array(grading, dtype='float32')
        self.hexas = hexas
        self.npoints = npoints
        self.grading = grading

        quads = array(quads, dtype='int32')
        inames = array(inames, dtype='int32')
        bcs = array(bcs, dtype='int32')

        self.quads = quads
        self.iname_to_quads = iname_to_quads
        self.inames = inames
        self.bcs = bcs
        self.iname_to_name = iname_to_name
        self.iname_to_type = iname_to_type

        #bdf_filename = 'blockMesh.bdf'
        #self.write_bdf(bdf_filename, nodes, hexas)
        return nodes, hexas, quads, inames, bcs

    def write_bdf(self, bdf_filename, nodes, hexas):
        """writes the BlockMesh as a Nastran BDF"""
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            for inode, node in enumerate(nodes):
                (x, y, z) = node
                fields = ['GRID', inode + 1, None, float(x), float(y), float(z)]
                bdf_file.write(print_card_8(fields))

            pid = 1
            for ielement, hexa in enumerate(hexas):
                bdf_file.write(print_card_8(['CHEXA', ielement + 1, pid,] + list(hexa)))

            bdf_file.write('PSOLID, 1, 1\n')
            bdf_file.write('MAT1, 1, 1.0,,0.3\n')
            bdf_file.write('ENDDATA\n')
        #print(write_dict(d, baseword='BlockMesh'))

    def adjust_nodes_to_symmetry(self):
        # find where nodes have -y value
        unused_y_locations = self.nodes[:, 1]

        #i = where(y_locations < -0.0508)[0]
        #self.nodes[i, 1] = -0.0508

        # set -y to 0
        #i = where(y_locations < 0.0)[0]
        #self.nodes[i, 1] = 0.0
        neq = 1
        ieq = 0
        while neq > 0:
            neq = self.equivalence_nodes(0.0001)
            #self.write_block_mesh('blockMeshDict_n' + str(neq), make_symmetry=True)
            ieq += 1
            if ieq == 10:
                asdf
                break
            self.log.debug('neq = %s' % neq)
        self.log.debug('ieq = %s' % ieq)

    def equivalence_nodes(self, Rtol=0.1):
        """equivalences the nodes of a BlockMeshDict"""
        #same_location = []
        inode_map = {}
        neq = 0
        nodes = self.nodes
        nnodes = self.nodes.shape[0]
        self.log.info('nnodes = %s' % nnodes)
        for inode, nodei in enumerate(self.nodes):
            for jnode, nodej in enumerate(self.nodes):
                if inode < jnode:
                    delta = norm(nodei-nodej)
                    if delta < Rtol:
                        #print('delta=%6g %2s %2s' % (delta, inode, jnode))
                        #same_location.append([inode, jnode])
                        inode_map[jnode] = inode
                        neq += 1
        #same_location = array(same_location, dtype='int32')
        #print('same_location = \n%s' % same_location)
        #nleft = nnodes - neq
        self.log.info('neq = %s' % neq)
        if neq > 15:
            asdf

        # get the node mapping as a dict
        new_ids_map = {}
        new_ids_map2 = {}
        for i in range(nnodes):
            if i in inode_map:  # nodes to adjust
                new_ids_map[i] = inode_map[i]
            else:
                new_ids_map[i] = i

        stack_nodes = set(self.hexas.flatten())
        self.log.info(stack_nodes)
        for iname, faces in iteritems(self.iname_to_quads):
            for face in faces:
                stack_nodes.update(set(list(face)))

        #if 0:
            #print(stack_nodes, len(stack_nodes))
            #for key, value in sorted(iteritems(new_ids_map)):
                #if key not in stack_nodes:
                    #print(' #k=%s v=%s' % (key, value))
                #if key != value:
                    #print(' *k=%s v=%s' % (key, value))
                #else:
                    #print('  k=%s v=%s' % (key, value))

        # map them to 0...n
        i = 0
        j = 0
        j0 = 0
        self.log.info('-------------------')
        for key, value in sorted(iteritems(new_ids_map)):  # the dict of collapsed nodes
            if key not in stack_nodes:
                j += 1
                continue

            if key == value:
                new_ids_map2[i+j] = value - j
                i += 1
                j0 = j
            else:
                new_ids_map2[i+j] = value - j0
                j += 1

        #for key, value in sorted(iteritems(new_ids_map2)):
            #self.log.info('  k=%s v=%s' % (key, value))
        new_ids_map = new_ids_map2

        # get new array of nodes
        #not_indexes = get_not_indexes(self.nodes, same_location[:, 1])
        #nodes2 = self.nodes[not_indexes, :]

        #a[not_indices] = 888
        #print('nodes2 =', nodes2)
        nodes2 = []
        nodes2_written = []
        for inode, jnode in sorted(iteritems(new_ids_map2)):
            if jnode not in nodes2_written:
                node2 = nodes[inode, :]
                nodes2.append(node2)
                nodes2_written.append(jnode)

        nodes2 = array(nodes2, dtype='float32')
        nnodes2 = nodes2.shape[0]
        assert nnodes2 <= nnodes, 'nnodes=%s nnodes2=%s' % (nnodes, nnodes2)

        # update the hexas
        hexas2 = []
        npoints2 = []  # don't change
        grading2 = []  # don't change
        for unused_ihexa, hexa, npointsi, gradingi in zip(count(), self.hexas,
                                                          self.npoints, self.grading):
            hexa2 = []
            for j in hexa:
                i = new_ids_map[j]
                hexa2.append(i)
            #print hexa2
            #asdf
            if len(unique(hexa2)) == 8:

                hexas2.append(hexa2)
                npoints2.append(npointsi)
                grading2.append(gradingi)
        hexas2 = array(hexas2, dtype='int32')
        npoints2 = array(npoints2, dtype='int32')
        grading2 = array(grading2, dtype='float32')


        # update the faces
        iname_to_quads = {}
        for iname, faces in iteritems(self.iname_to_quads):
            faces2 = []
            for face in faces:
                face2 = []
                for j in face:
                    i = new_ids_map[j]
                    face2.append(i)
                faces2.append(face2)
            iname_to_quads[iname] = array(faces2, dtype='int32')

        #self.log.info(nodes.shape)
        #self.log.info(nodes2.shape)
        self.nodes = nodes2
        self.hexas = hexas2
        self.grading = grading2
        self.npoints = npoints2
        self.iname_to_quads = iname_to_quads

        return neq

    def write_block_mesh(self, block_mesh_name_out='blockMeshDict.out', make_symmetry=False):
        """filename interface to ``_write_block_mesh``"""
        with open(block_mesh_name_out, 'w') as block_mesh_file:
            self.log.info('writing %s' % block_mesh_name_out)
            self._write_block_mesh(block_mesh_file, make_symmetry)

    def _write_block_mesh(self, block_mesh_file, make_symmetry):
        """writes a BlockMeshDict with a file object"""
        nodes = self.nodes
        hexas = self.hexas

        header = (
            '/*--------------------------------*- C++ -*----------------------------------*\\\n'
            '| =========                 |                                                 |\n'
            '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n'
            '|  \\    /   O peration     | Version:  2.2.2                                 |\n'
            '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n'
            '|    \\/     M anipulation  |                                                 |\n'
            '\*---------------------------------------------------------------------------*/\n'
            'FoamFile\n'
            '{\n'
            '    version     0.0508;\n'
            '    format      ascii;\n'
            '    class       dictionary;\n'
            '    object      blockMeshDict;\n'
            '}\n'
            '\n'
            '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'
        )
        block_mesh_file.write(header)

        block_mesh_file.write('convertToMeters 1.0;\n')
        block_mesh_file.write('\n')
        block_mesh_file.write('vertices\n')
        block_mesh_file.write('(\n')

        unique_x = unique(nodes[:, 0])
        unique_y = unique(nodes[:, 1])
        unique_z = unique(nodes[:, 2])
        unique_x.sort()
        unique_y.sort()
        unique_z.sort()
        #print('unique_x = %s' % unique_x)
        self.log.info('unique_y = %s' % unique_y)
        #print('unique_z = %s' % unique_z)

        nnodes = nodes.shape[0]
        fmt_node = '%%%si' % len(str(nnodes))
        #print('fmt-node = %r' % fmt_node)
        for inode, (x, y, z) in enumerate(nodes):
            block_mesh_file.write('    (%7s %7s %7s) // %s\n' % (x, y, z, inode))
            if (inode + 1) % 5 == 0:
                block_mesh_file.write('\n')
        block_mesh_file.write(');\n\n')

        block_mesh_file.write('blocks\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write('      // hex                        npoints in each dir;    stretching\n')
        #print "nodes = ", nodes.shape
        hexai_fmt = 'hex (%s %s %s %s %s %s %s %s)' % tuple([fmt_node] * 8)

        m_to_inch = 39.3701
        for ihexa, (hexa, npointsi, gradingi) in enumerate(zip(hexas, self.npoints, self.grading)):
            save_element = False
            for inode in hexa:
                #print('inode = %s' % inode)
                node = nodes[inode, :]
                #print('nodes[%s] = %s' % (inode, str(node)))
                if node[1] >= 0.0:
                    save_element = True

            #npointsi = self.npoints[ihexa]
            #gradingi = self.grading[ihexa]
            save_element = True
            if save_element:
                voli = volume8(
                    nodes[hexa[0], :], nodes[hexa[1], :], nodes[hexa[2], :], nodes[hexa[3], :],
                    nodes[hexa[4], :], nodes[hexa[5], :], nodes[hexa[6], :], nodes[hexa[7], :],)
                if voli > 0:
                    shexai = hexai_fmt % tuple(hexa)
                    snpointsi = '(%s %s %s)' % tuple(npointsi)
                    sgradingi = 'simpleGrading (%g %g %g)' % tuple(gradingi)
                    svoli = 'vol = %g [in]' % (voli * m_to_inch**3)
                    block_mesh_file.write('      %s %s %s  // %s\n' % (
                        shexai, snpointsi, sgradingi, svoli))

            if (ihexa + 1) % 4 == 0:
                block_mesh_file.write('\n')
        block_mesh_file.write(');\n\n')

        block_mesh_file.write('edges\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write(');\n')

        block_mesh_file.write('\n')
        block_mesh_file.write('boundary\n')
        block_mesh_file.write('(\n')

        #iname_quads = {}

        symmetry_faces = []
        facei_fmt = '           (%s %s %s %s)' % tuple([fmt_node] * 4)
        for iname in sorted(self.iname_to_quads):
            name = self.iname_to_name[iname]
            Type = self.iname_to_type[iname]
            block_mesh_file.write('    %s\n' % name)
            block_mesh_file.write('    {\n')
            block_mesh_file.write('        type %s;\n' % Type)
            block_mesh_file.write('        faces\n')
            block_mesh_file.write('        (\n')
            faces = self.iname_to_quads[iname]
            for face in faces:
                n1 = nodes[face[0], :]
                n2 = nodes[face[1], :]
                n3 = nodes[face[2], :]
                n4 = nodes[face[3], :]
                area, centroid = area_centroid(n1, n2, n3, n4)
                centroid *= m_to_inch
                if area > 0.:
                    area *= m_to_inch ** 2
                    y_centroid = centroid[1]
                    if make_symmetry and allclose(y_centroid, 0.0):
                        symmetry_faces.append(face)
                    else:
                        if y_centroid <= 0.0:
                            block_mesh_file.write(
                                facei_fmt % tuple(face) +
                                ' // centroid=(%3g, %3g, %3g) [in]  Area=%.2f [in^2]\n' % (
                                    centroid[0], centroid[1], centroid[2], area))
                        else:
                            block_mesh_file.write(facei_fmt % tuple(face) + '\n') #+ ' // c=(%7s, %7s, %7s)\n' % tuple(centroid)
            block_mesh_file.write('        );\n')
            block_mesh_file.write('    }\n')
            #save_face = False

        if make_symmetry:
            #name = self.iname_to_name[iname]
            name = 'symmetry'
            foam_type = 'symmetryPlane'
            #foam_type = self.iname_to_type[iname]

            block_mesh_file.write('    %s\n' % name)
            block_mesh_file.write('    {\n')
            block_mesh_file.write('        type %s;\n' % foam_type)
            block_mesh_file.write('        faces\n')
            block_mesh_file.write('        (\n')
            for face in symmetry_faces:
                n1 = nodes[face[0], :]
                n2 = nodes[face[1], :]
                n3 = nodes[face[2], :]
                n4 = nodes[face[3], :]
                area, centroid = area_centroid(n1, n2, n3, n4)
                centroid *= m_to_inch
                area *= m_to_inch ** 2

                y_centroid = centroid[1]
                if y_centroid <= 0.0:
                    block_mesh_file.write(
                        facei_fmt % tuple(face) +
                        ' // centroid=(%3g, %3g, %3g) [in]  Area=%.2f [in^2]\n' % (
                            centroid[0], centroid[1], centroid[2], area))
                else:
                    block_mesh_file.write(facei_fmt % tuple(face) + '\n')
                    #+ ' // c=(%7s, %7s, %7s)\n' % tuple(centroid)

            block_mesh_file.write('        );\n')
            block_mesh_file.write('    }\n')

        block_mesh_file.write(');\n\n')
        block_mesh_file.write('mergeMatchPairs\n')
        block_mesh_file.write('(\n')
        block_mesh_file.write(');\n\n')
        block_mesh_file.write('// ************************************************************************* //\n')

def mirror_block_mesh(block_mesh_name, block_mesh_name_out, log=None, debug=True):
    """mirrors a blockMeshDict"""
    make_symmetry = True
    block_mesh_model = read_block_mesh(block_mesh_name, log=log, debug=debug)
    #out = block_mesh_model.read_openfoam(block_mesh_name)
    #unused_nodes, unused_hexas, unused_quads, unused_names, unused_bcs = out
    if make_symmetry:
        block_mesh_model.adjust_nodes_to_symmetry()
    block_mesh_model.write_block_mesh(block_mesh_name_out, make_symmetry=make_symmetry)

def main():  # pragma: no cover
    import sys
    block_mesh_name = sys.argv[1]
    block_mesh_name_out = sys.argv[2]
    mirror_block_mesh(block_mesh_name, block_mesh_name_out)

if __name__ == '__main__':  # pragma: no cover
    main()
