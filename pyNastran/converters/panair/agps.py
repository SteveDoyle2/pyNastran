from __future__ import print_function
#from six import iteritems
from numpy import array

class AGPS(object):
    def __init__(self, log=None, debug=False):
        self.infilename = None
        self.pressures = {}
        self.log = log
        self.debug = False

    def read_agps(self, infilename):
        self.infilename = infilename
        with open(self.infilename, 'r') as f:
            lines = f.readlines()

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
                    #print("saving column")
                    # still on same patch
                    patch.append(col)
                else:
                    #print("saving patch")
                    patch.append(col)
                    patches.append(patch)
                    # new patch
                    patch = []
                    ipatch += 1
                    #print("new patch", ipatch)
                #print("lines[%i] = %r" % (i, (n, c)))
                i += 1
            col = []
        #print("*saving patch")
        patches.append(patch)

        #print("len(patches) =", len(patches))
        #assert len(patches) < 20
        # time to parse the patches
        #for ipatch, patch in enumerate(patches):

        for ipatch, patch in enumerate(patches):
            if self.debug:
                print("ipatch=%s" % ipatch)
            #print("ipatch =", ipatch)
            #print(patch)
            X = []
            Y = []
            Z = []
            Cp = []
            for icol, col in enumerate(patch):
                x = []
                y = []
                z = []
                cp = []

                #print("next col")
                for inode, node in enumerate(col):
                    #print(node)
                    # dropping the counter with [1:]
                    xi, yi, zi, cpi = node.strip().split()[1:]
                    x.append(float(xi))
                    y.append(float(yi))
                    z.append(float(zi))
                    cp.append(float(cpi))
                X.append(x)
                Y.append(y)
                Z.append(z)
                Cp.append(cp)
                #print("len(cp) = %s" % len(cp))
            #print(X)
            if self.debug:
                print(len(X))
                for x in X:
                    print('%s %s' % (x, len(x)))

            X = array(X, dtype='float32')
            Y = array(Y, dtype='float32')
            Z = array(Z, dtype='float32')
            Cp = array(Cp, dtype='float32')
            self.pressures[ipatch] = Cp
            #print(Cp.shape)
            if self.debug:
                print("")
        #for ipatch, Cp in sorted(iteritems(self.pressures)):
            #print(Cp)

if __name__ == '__main__':  # pragma: no cover
    a = AGPS()
    a.read_agps('agps')
