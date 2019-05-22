"""reads the panair.out file"""
import numpy as np
from cpylog import get_logger2

def read_panair_out(panair_out_filename='panair.out', log=None, debug=False):
    """reads the panair.out file"""
    model = PanairOut(log=log, debug=debug)
    model.read_panair_out(panair_out_filename=panair_out_filename)
    return model


class Ft13Network:
    """Stores network/patch info"""
    def __init__(self, inetwork, data):
        #print(data[0])
        self.inetwork = inetwork
        self.data = np.vstack(data)
        #print(self.data.shape)
        #print(self.data)
        #print('ft13: %i: %s' % (inetwork, len(data)))  # 720

    def __repr__(self):
        return '<Ft13Network>; i=%s len=%s' % (self.inetwork, len(self.data))

class Network:
    """Stores network/patch info"""
    def __init__(self, line):
        self.line = line
        self.data = []
        self.jc = 0
        self.ip = 0
        self.headers = [
            'x', 'y', 'z', 'wx', 'wy', 'wz', 'cp2ndu', 'cpisnu', 'lmachu', 'source', 'doublet',
        ]

    def add(self, jc, ip, x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet):
        """adds the data"""
        self.jc = jc
        self.ip = ip
        self.data.append([x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet])

    def to_numpy(self):
        """converts to a numpy array in point order"""
        self.data = np.array(self.data, dtype='float32')#.reshape(jc)

class PanairOut:
    """reads the panair.out file"""
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log=log, debug=debug, encoding='utf-8')
        self.headers = [
            'x', 'y', 'z', 'wx', 'wy', 'wz', 'cp2ndu', 'cpisnu', 'lmachu', 'source', 'doublet'
        ]
        self.headers_ft13 = [
            'x', 'y', 'z', 'd0', 'dx', 'dy', 'dz', 's0', 'anx', 'any', 'anz',
            'lmachu', 'wxu', 'wyu', 'wzu', 'pheu', 'vxu', 'vyu', 'vzu', 'cplinu', 'cpslnu', 'cp2ndu', 'cpisnu',
            'lmachl', 'wxl', 'wyl', 'wzl', 'phel', 'vxl', 'vyl', 'vzl', 'cplinl', 'cpslnl', 'cp2ndl', 'cpisnl',
            'wnu', 'wnl', 'pwnu', 'pwnl', 'vtu', 'vtl', 'pvtu', 'pvtl', 'cplind', 'cpslnd', 'cp2ndd', 'cpisnd']
        self.networks = {}
        self.networks_ft13 = {}

    def read_ft13(self, ft13_filename='ft13'):
        with open(ft13_filename, 'r') as out_file:
            lines = out_file.readlines()

        iline = 0
        nlines = len(lines) - 1
        while iline < nlines:
            line = lines[iline]
            #print(line)
            if 'jc   ip        x          y          z          d0         dx         dy         dz         s0        anx        any        anz' in line:
                #print('break')
                break
            iline += 1

        key_lines = [
            'jc   ip        x          y          z          d0         dx         dy         dz         s0        anx        any        anz',
            'lmachu        wxu        wyu        wzu       pheu       vxu        vyu        vzu       cplinu     cpslnu     cp2ndu     cpisnu',
            'lmachl        wxl        wyl        wzl       phel       vxl        vyl        vzl       cplinl     cpslnl     cp2ndl     cpisnl',
            'wnu        wnl        pwnu       pwnl        vtu        vtl        pvtu       pvtl      cplind     cpslnd     cp2ndd     cpisnd',
        ]

        inetwork = 0
        jcounter = 0
        networks = {}
        data = []
        datai = []
        #print(iline, nlines)
        isolution = 1
        all_networks = {isolution : networks}
        while iline < nlines:
            line = lines[iline].rstrip()
            if not line:
                iline += 1
                continue
            if len(line) == 1:
                if len(datai) and datai[0] != 'network':
                    #print('datai = ', datai)
                    data.append(datai)
                datai = []
                iline += 1
                continue
            anyi = any([key_line in line for key_line in key_lines])
            #print(anyi)
            if anyi:
                #print(line)
                #print(datai)
                assert jcounter == 0, jcounter
                assert datai == [], data
                iline += 1
                continue

            if 'network id:' in line:
                if datai and datai[0] != 'network':
                    #print('datai =', len(datai), datai)
                    data.append(datai)
                if len(data) < 4:
                    self.log.info(str(data))
                inetwork += 1
                if len(data):
                    assert len(data) > 0, data
                    networks[inetwork] = Ft13Network(inetwork, data)
                data = []
                datai = []

            assert 'jc' not in line, line
            assert 'lmachu' not in line, line
            #print(iline, jcounter, line)
            iline += 1
            if 'mach number =' in line:
                iline += 2
                isolution += 1
                inetwork = 0
                networks = {}
                all_networks[isolution] = networks
                assert len(datai) == 0, datai
                assert len(data) == 0
                continue

            sline = line.split()
            if sline[0] == 'network':
                continue

            # 12
            if jcounter == 0:
                sline2 = [
                    #line[:11],
                    line[11:22], line[22:33], line[33:44], line[44:55], line[55:66], line[66:77],
                    line[77:88], line[88:99], line[99:110], line[110:121], line[121:132], #line[132:143],
                ]
            else:
                sline0 = line[:11].strip()
                sline2 = [
                    line[:11],
                    line[11:22], line[22:33], line[33:44], line[44:55], line[55:66], line[66:77],
                    line[77:88], line[88:99], line[99:110], line[110:121], line[121:132], #line[132:143],
                ]
                assert ' ' not in sline0, sline2
            assert ' ' not in sline2[0].strip(), sline2
            #print(sline)
            #assert len(sline) == len(sline2), sline2
            datai += sline2

            #print(jcounter, datai)
            jcounter += 1
            if jcounter == 4:
                data.append(datai)
                assert len(datai) == len(self.headers_ft13)
                datai = []
                jcounter = 0
        self.log.debug('all_networks=%s' % str(all_networks.keys()))
        self.log.debug('  networks=%s' % str(networks.keys()))
        self.networks_ft13 = all_networks
        return all_networks

    def read_panair_out(self, panair_out_filename='panair.out'):
        with open(panair_out_filename, 'r') as out_file:
            lines = out_file.readlines()

        iline = -1
        unused_keywords = [
            'simultaneous solution number',
            'network id:',
            'force / moment data for network']
        isolution = 1
        iline_old = -1
        nlines = len(lines) - 1
        while iline < nlines:
            #print('------------------------------------')
            #print(iline, nlines)
            iline = self.get_solution(lines, iline, nlines, isolution=isolution)
            if not iline > iline_old:  # pragma: no cover
                raise RuntimeError('did not increment counter...')
            iline_old = iline
            isolution += 1
        self.log.debug('done with panair.out')

        for isolution, networks in sorted(self.networks.items()):
            self.log.debug('out: isolution=%s' % isolution)
            for inetwork, network in networks.items():
                network.to_numpy()
                self.log.debug('  %s %s' % (inetwork, str(network.data.shape)))

    def get_solution(self, lines, iline, nlines, isolution=1):
        networks = {}
        self.networks[isolution] = networks
        if iline == nlines:
            #print('%s/%s: %s' % (iline, nlines, line))
            return iline

        line = lines[iline].rstrip()
        while 'simultaneous solution number' not in line and iline < nlines:
            line = lines[iline].rstrip()
            iline += 1
        #print('%s: %s' % (iline, line))

        if iline == nlines:
            self.log.debug('%s/%s: %s' % (iline, nlines, line))
            return iline

        inetwork = 1
        self.log.debug('get to start of network %i' % inetwork)
        # get to start of first network
        while not line.startswith('1    network id:') and iline < nlines:
            line = lines[iline].rstrip()
            iline += 1

        if iline == nlines:
            self.log.debug('%s/%s: %s' % (iline, nlines, line))
            return iline

        if line.startswith('1    network id:') and inetwork not in networks:
            network = Network(line)
            networks[inetwork] = network
            #print(network)
            self.log.debug('%s: network!; %s' % (iline, line))

        #print('reading network---')
        #print(line)
        #self.log.info('while force/moment')
        while 'force / moment data for network' not in line:
            line = lines[iline].rstrip()
            if 'simultaneous solution number' in line:
                return iline
            #print(line)
            if 'force / moment data for network' in line:
                #self.log.info('  while network')
                while '1    network id:' not in line and 'simultaneous solution number' not in line and iline < nlines:
                    iline += 1
                    line = lines[iline].rstrip()
                    if 'simultaneous solution number' in line:
                        return iline
                    #self.log.debug("  breaking on 'force / moment data for network'")
                    #break

                #print(line)
                if line.startswith('1    network id:'):
                    inetwork += 1
                    if inetwork not in networks:
                        #self.log.info('***inetwork=%s' % inetwork)
                        network = Network(line)
                        networks[inetwork] = network

                #continue

            if iline == nlines:
                self.log.debug('%s/%s: %s' % (iline, nlines, line))
                break
                #return iline

            if line.startswith('1    network id:'):
                if inetwork not in networks:
                    network = Network(line)
                    networks[inetwork] = network
                    #print(network)
                    self.log.info('%s: network!; %s' % (iline, line))
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
                #print(iline, line)
                raise
            network.add(jc, ip, x, y, z, wx, wy, wz, cp2ndu, cpisnu, lmachu, source, doublet)
            iline += 1
        #print('return %s: %s' % (iline, line))
        line = lines[iline].rstrip()
        inetwork += 1
        return iline


if __name__ == '__main__':  # pragma: no cover
    read_panair_out(panair_out_filename='panair.out')
