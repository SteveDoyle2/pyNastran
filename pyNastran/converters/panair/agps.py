"""
defines:
 - AGPS(log=None, debug=False)

"""
import numpy as np

class AGPS:
    """Interface to the AGPS file"""
    def __init__(self, log=None, debug=False):
        """
        Initializes the AGPS object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has

        """
        self.infilename = None
        self.pressures = {}
        self.log = log
        self.debug = debug

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
            self.log.debug('ipatch=%s' % ipatch)

            nrows = len(patch)
            ncols = len(patch[0])
            XYZ = np.zeros((nrows, ncols, 3), dtype='float32')

            cps = []
            node0 = patch[0][0].split()
            cp0 = node0[4:]
            ncp = len(cp0)
            cp = np.zeros((ncp, nrows, ncols), dtype='float32')
            for icol, col in enumerate(patch):
                for inode, node in enumerate(col):
                    sline = node.strip().split()

                    # dropping the counter with [1:]
                    # counter, xi, yi, zi, cp1, cp2, cp3, cp4
                    xyz = sline[1:4]
                    cpi = sline[4:]
                    assert len(xyz) == 3, 'sline=%s xyz=%s' % (sline, xyz)
                    assert 1 <= len(cp) <= 3, 'sline=%s cpi=%s' % (sline, cpi)
                    XYZ[icol, inode, :] = xyz # xi, yi, zi

                    # we don't support multiple Cps, but you can hack it...
                    cp[:, icol, inode] = cpi

            self.pressures[ipatch] = cp
        #for ipatch, Cp in sorted(self.pressures.items()):
            #print(Cp)
