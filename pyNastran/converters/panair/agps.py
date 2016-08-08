from __future__ import print_function
#from six import iteritems
import numpy as np

class AGPS(object):
    def __init__(self, log=None, debug=False):
        self.infilename = None
        self.pressures = {}
        self.log = log
        self.debug = False

    def read_agps(self, infilename):
        """Loads an AGPS file"""
        self.infilename = infilename
        with open(self.infilename, 'r') as agps_file:
            lines = agps_file.readlines()

        i = 0
        while i < len(lines):
            if 'irow' in lines[i] and 'cp1' in lines[i]:
                break
            i += 1
        i += 1

        #print("lines[%i] = %r" % (i, lines[i]))
        col = []
        patches = []
        patch = []
        ipatch = 1
        while i < len(lines):
            while '*eof' not in lines[i]:
                col.append(lines[i].strip())
                #print(lines[i].strip())
                i += 1
            i += 1
            if i < len(lines):
                n, c = lines[i].strip('\n\r n').split('c')
                n, c = int(n), int(c)
                if n == ipatch:
                    #print('saving column')
                    # still on same patch
                    patch.append(col)
                else:
                    #print('saving patch')
                    patch.append(col)
                    patches.append(patch)
                    # new patch
                    patch = []
                    ipatch += 1
                    #print('new patch', ipatch)
                #print('lines[%i] = %r' % (i, (n, c)))
                i += 1
            col = []
        patches.append(patch)

        # time to parse the patches
        for ipatch, patch in enumerate(patches):
            if self.debug:
                print('ipatch=%s' % ipatch)

            nrows = len(patch)
            ncols = len(patch[0])
            XYZ = np.zeros((nrows, ncols, 3), dtype='float32')
            Cp = np.zeros((nrows, ncols), dtype='float32')
            for icol, col in enumerate(patch):
                for inode, node in enumerate(col):
                    # dropping the counter with [1:]
                    xi, yi, zi, cpi = node.strip().split()[1:]
                    XYZ[icol, inode, :] = [xi, yi, zi]
                    Cp[icol, inode] = cpi

            self.pressures[ipatch] = Cp
            if self.debug:
                print('')
        #for ipatch, Cp in sorted(iteritems(self.pressures)):
            #print(Cp)

def main():
    agps = AGPS()
    agps.read_agps('agps')

if __name__ == '__main__':  # pragma: no cover
    main()
