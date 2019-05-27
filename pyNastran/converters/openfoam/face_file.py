"""
defines:
 - model = FaceFile(log=None, debug=False)

"""
from copy import deepcopy

import numpy as np

from cpylog import get_logger2


class FaceFile:
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
            faces = self._read_all_faces(face_file, nfaces)
        else:
            faces = self._read_subset_faces(face_file, nfaces, ifaces_to_read)
        self.log.info('faces.shape = %s' % str(faces.shape))
        return faces

    def _read_all_faces(self, face_file, nfaces):
        """reads all the faces"""
        self.log.info('nfaces = %s' % nfaces)
        faces = np.ones((nfaces, 4), dtype='int32') * -1

        i = 0
        for j in range(nfaces):
            # 3(a b c) to [a, b, c]
            # 4(a b c d) to [a, b, c, d]
            face_line = face_file.readline()
            i += 1
            sline = face_line[1:].strip('( )\n\r').split()
            #self.log.info(str(sline))
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
        return faces

    def _read_subset_faces(self, face_file, nfaces, ifaces_to_read):
        """
        we need to sort our faces as we read the faces sequentially
        however the output needs to be in unsorted order so we save the sorting
        order and back apply it at the end
        """
        if isinstance(ifaces_to_read, list):
            ifaces_to_read = np.array(ifaces_to_read)

        isort = np.argsort(ifaces_to_read)
        self.log.info('isort = %s' % isort)
        ifaces_to_read_sorted = ifaces_to_read[isort]

        nfaces2 = len(ifaces_to_read)  # must be sorted
        assert len(isort) == nfaces2
        faces = np.ones((nfaces2, 4), dtype='int32') * -1
        ni = 0
        i = 0
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
        return faces
