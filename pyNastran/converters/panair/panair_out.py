"""reads the panair.out file"""
from pyNastran.utils.log import get_logger2
import numpy as np

def read_panair_out(panair_out_filename='panair.out'):
    """reads the panair.out file"""
    model = PanairOut()
    model.read_panair_out(panair_out_filename=panair_out_filename)
    return model

class Network(object):
    """Stores network/patch info"""
    def __init__(self, line):
        self.line = line
        self.data = []
        self.jc = 0
        self.ip = 0
        self.headers = ['x', 'y', 'z', 'wx', 'wy', 'wz', 'cp2ndu', 'cpisnu', 'lmachu', 'source', 'doublet']

    def add(self, jc, ip, x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet):
        """adds the data"""
        self.jc = jc
        self.ip = ip
        self.data.append([x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet])

    def to_numpy(self):
        """converts to a numpy array in point order"""
        self.data = np.array(self.data, dtype='float32')#.reshape(jc)

class PanairOut(object):
    """reads the panair.out file"""
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log=log, debug=debug, encoding='utf-8')
        self.headers = ['x', 'y', 'z', 'wx', 'wy', 'wz', 'cp2ndu', 'cpisnu', 'lmachu', 'source', 'doublet']
        self.networks = {}

    def read_panair_out(self, panair_out_filename='panair.out'):
        with open(panair_out_filename, 'r') as out_file:
            lines = out_file.readlines()

        iline = -1
        keywords = ['simultaneous solution number', 'network id:', 'force / moment data for network']
        isolution = 1
        iline_old = -1
        nlines = len(lines) - 1
        while iline < nlines:
            print('------------------------------------')
            print(iline, nlines)
            iline = self.get_solution(lines, iline, nlines, isolution=isolution)
            if not iline > iline_old:  # pragma: no cover
                raise RuntimeError('did not increment counter...')
            iline_old = iline
            isolution += 1
        self.log.debug('done with panair.out')

        for isolution, networks in self.networks.items():
            for inetwork, network in networks.items():
                network.to_numpy()

    def get_solution(self, lines, iline, nlines, isolution=1):
        networks = {}
        self.networks[isolution] = networks
        if iline == nlines:
            print('%s/%s: %s' % (iline, nlines, line))
            return iline

        line = lines[iline].rstrip()
        while 'simultaneous solution number' not in line and iline < nlines:
            line = lines[iline].rstrip()
            iline += 1
        print('%s: %s' % (iline, line))

        if iline == nlines:
            print('%s/%s: %s' % (iline, nlines, line))
            return iline

        inetwork = 1
        print('get to start of network %i' % inetwork)
        # get to start of first network
        while not line.startswith('1    network id:') and iline < nlines:
            line = lines[iline].rstrip()
            iline += 1

        if iline == nlines:
            print('%s/%s: %s' % (iline, nlines, line))
            return iline

        if line.startswith('1    network id:'):
            if inetwork not in networks:
                network = Network(line)
                networks[inetwork] = network
                #print(network)
                #print('%s: network!' % iline, line)

        #print('reading network---')
        #print(line)
        while 'force / moment data for network' not in line:
            line = lines[iline].rstrip()
            #print(line)
            if 'force / moment data for network' in line:
                iline += 1
                break
            if line.startswith('1    network id:'):
                if inetwork not in networks:
                    network = Network(line)
                    networks[inetwork] = network
                    #print(network)
                    #print('%s: network!' % iline, line)
                iline += 1
                continue
            elif 'jc   ip         x          y          z         wx         wy         wz      cp2ndu     cpisnu     lmachu     source    doublet' in line:
                iline += 1
                continue

            if not line:
                iline += 1
                continue

            if line[0] == '1' and len(line) == 1 or '0*b*for-mom-net#-' in line:
                #print("****************************")
                iline += 1
                continue

            sline = line.split()
            #print(iline, sline)
            try:
                jc, ip, x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet = sline
            except:
                print(iline, line)
                raise
            network.add(jc, ip, x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet)
            iline += 1
        #print('*%s: %s' % (iline, line))
        line = lines[iline].rstrip()
        inetwork += 1
        return iline


if __name__ == '__main__':  # pragma: no cover
    read_panair_out(panair_out_filename='panair.out')


