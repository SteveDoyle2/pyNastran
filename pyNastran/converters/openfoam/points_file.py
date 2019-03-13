"""
defines:
 - model = PointFile(log=None, debug=False)
"""
from __future__ import print_function
from codecs import open

import numpy as np

from cpylog import get_logger2


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
                points = np.zeros((npoints2, 3), dtype='float32')

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
                points = np.zeros((npoints, 3), dtype='float32')
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
