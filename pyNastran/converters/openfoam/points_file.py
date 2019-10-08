"""
defines:
 - model = PointFile(log=None, debug=False)

"""
import numpy as np

from cpylog import get_logger2


def read_points_file(point_filename, ipoints_to_read=None,
                    log=None, debug=False):
    """functional interface to PointFile"""
    points = PointFile(log=log, debug=debug)
    xyz = points.read_point_file(point_filename, ipoints_to_read=ipoints_to_read)
    return xyz


class PointFile:
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_point_file(self, point_filename, ipoints_to_read=None):
        """rads the points file"""
        #p = FoamFile(face_filename)
        #lines = p.read_foam_file()
        #self.log.info('converting')
        i = 0
        npoints = 0
        iempty = 0
        with open(point_filename, 'r') as points_file:
            while npoints == 0:
                line = points_file.readline().rstrip()
                if len(line) == 0:
                    iempty += 1
                    assert iempty < 50, 'no lines found in the file'
                i += 1
                npoints = int(line)
            line = points_file.readline()
            i += 1

            self.log.info('npoints = %s' % npoints)
            #print('lineA = %r' % line)

            self.log.debug('building points')
            assert npoints > 0, npoints
            if ipoints_to_read is not None:
                #print('ipoints =', ipoints_to_read)
                ipoints_to_read.sort()
                npoints_to_read = len(ipoints_to_read)
                points = np.zeros((npoints_to_read, 3), dtype='float32')

                self.log.info('npoints_to_read = %s' % npoints_to_read)
                j = -1
                ni = 0

                for ipoint in ipoints_to_read:
                    # get the next point id
                    point = points_file.readline()
                    sline = point.strip('( )\n\r').split()
                    i += 1
                    j += 1
                    #print(f'j={j} ni={ni} ipoint={ipoint}')

                    while j < ipoint:
                        #print('skipping %s' % sline)
                        point = points_file.readline()
                        sline = point.strip('( )\n\r').split()
                        assert len(sline) == 3, f'iline={i} sline={sline}'
                        i += 1
                        j += 1
                    #print(point)

                    sline = point.strip('( )\n\r').split()
                    #print('sline =', sline)
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
                    #print(j, sline)
                    try:
                        points[j, :] = sline
                    except:
                        print('point = %r' % point)
                        print(sline, i)
                        raise
        self.log.info('points.shape = %s' % str(points.shape))
        return points
