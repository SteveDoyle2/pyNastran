"""
defines:
 - model = Boundary(log=None, debug=False)
 - model = BoundaryFile(log=None, debug=False)

"""
import os
from collections import OrderedDict

import numpy as np
from cpylog import get_logger2

from pyNastran.converters.openfoam.points_file import PointFile
from pyNastran.converters.openfoam.face_file import FaceFile
from pyNastran.utils import check_path

def read_boundary(point_filename, face_filename, boundary_filename,
                  log=None, debug=False):
    boundary = Boundary(log=log, debug=debug)
    boundary.read_openfoam(point_filename, face_filename, boundary_filename)
    return boundary

def read_boundary_file(boundary_filename, log=None, debug=False):
    boundary_file = BoundaryFile(log=log, debug=debug)
    boundary_file.read_boundary_file(boundary_filename)
    return boundary_file

class BoundaryFile:
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_boundary_file(self, boundary_filename):
        #p = FoamFile(face_filename)
        #lines = p.read_foam_file()
        #self.log.info('converting')
        i = 0
        with open(boundary_filename, 'r') as boundary_file:
            lines = boundary_file.readlines()

        nboundaries = 0
        while nboundaries == 0:
            lineb = lines[i]
            i += 1
            try:
                nboundaries = int(lineb)
            except ValueError:
                print(lineb)
                raise
        lineb = lines[i]
        i += 1

        self.log.info('nboundaries = %s' % nboundaries)
        #self.log.info('lineA = %r' % line)

        self.log.info('building boundaries')
        boundaries = OrderedDict()
        boundaries = self._read_boundaries(lines, i, nboundaries, boundaries)

        return boundaries

    def _read_boundaries(self, lines, i, nboundaries, boundaries, basename=''):
        assert nboundaries > 0, nboundaries
        #read_next = 0
        for unused_j in range(nboundaries):
            #print('iboundary', unused_j)
            # 3(a b c) to [a, b, c]
            # 4(a b c d) to [a, b, c, d]
            nameline = lines[i]
            i += 1
            name = nameline.strip()
            self.log.info('name = %r' % (basename + name))

            unused_openline = lines[i]
            i += 1
            typeline = lines[i]
            i += 1
            word, boundary_type = typeline.strip('\n\r\t ;').split()

            nfacesline = lines[i]
            i += 1
            sline = nfacesline.strip('\n\r\t ;').split()
            word = sline[0]
            if word == 'inGroups' and len(sline) == 1:
                #nboundary_line = lines[i]; i+=1
                #nboundaries2 = int(nboundary_line)
                #openline = boundary_file.readline(); i+=1
                for unused_ii in range(7):
                    groupline = lines[i]
                    i += 1
                    #self.log.info(ii, groupline)
                sline = groupline.strip('\n\r\t ;').split()

                word, nfaces = sline
                nfaces = int(nfaces)
                self.log.info('nfaces = %r' % nfaces)

                startfacesline = lines[i]
                i += 1
                self.log.info('startfacesline = %r' % startfacesline)
                word, startfaces = startfacesline.strip('\n\r\t ;').split()
                startfaces = int(startfaces)
                unused_closeline = lines[i]
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
                    nfacesline = lines[i]
                    i += 1
                    word, nfaces = nfacesline.strip('\n\r\t ;').split()
                else:
                    word, nfaces = sline

                nfaces = int(nfaces)
                self.log.info('nfaces = %r' % (nfaces))

                startfacesline = lines[i]
                i += 1
                self.log.info('startfacesline = %r' % startfacesline)
                word, startfaces = startfacesline.strip('\n\r\t ;').split()
                startfaces = int(startfaces)

                unused_closeline = lines[i]
                i += 1

                boundary_name = basename + name
                if boundary_name in boundaries:
                    raise KeyError('boundary_name=%r is already defined...'
                                   'boundaries must have unique names' % boundary_name)
                boundaries[boundary_name] = [boundary_type, nfaces, startfaces]

        for name, boundary in boundaries.items():
            self.log.info('name=%s boundary=%s' % (name, boundary))
        return boundaries


class Boundary:
    """defines Boundary"""
    def __init__(self, log=None, debug=False):
        """creates a Boundary object"""
        self.debug = False
        #log = None
        self.log = get_logger2(log, debug=debug)

    def read_openfoam(self, point_filename, face_filename, boundary_filename):
        """reads a Boundary file"""
        check_path(face_filename, 'face_filename')
        check_path(point_filename, 'point_filename')
        check_path(boundary_filename, 'boundary_filename')

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
        for name, boundary in boundaries.items():
            # type            patch;  # 0
            # nFaces          nFaces; # 1
            # startFace       777700; # 2
            self.log.info('boundary[%s] = %s' % (name, boundary))
            nfacesi = boundary[1]
            startface = int(boundary[2])
            nfaces2 += nfacesi
            new_faces = list(np.arange(nfacesi, dtype='int32') + startface)
            #f_boundary_faces.write('boundary_faces[%s, %s] = %s\n' % (
                #name, len(new_faces), new_faces))
            ifaces_to_read += new_faces

        self.log.info('nfaces2 = %s' % nfaces2)
        ifaces_to_read = np.ravel(ifaces_to_read)
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
            ipoints_to_read = np.unique(faces.ravel())
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
        names = np.zeros(nfaces2, dtype='int32')
        iname = 1
        snames = [None] * (len(boundaries) + 1)
        self.log.info('')
        for name, boundary in boundaries.items():
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
            #pid_save = set()
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
