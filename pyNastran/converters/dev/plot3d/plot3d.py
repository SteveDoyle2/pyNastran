from numpy import zeros
from cpylog import __version__ as CPYLOG_VERSION
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger


class Plot3d:  # pragma: no cover
    def __init__(self, log=None, debug=False):
        self.x = {}
        self.y = {}
        self.z = {}
        self.block_shapes = {}

        self.log = get_logger(log, debug)

    def read_plot3d(self, p3d_name):  # pragma: no cover
        self.read_plot3d_ascii(p3d_name)

    def read_plot3d_ascii(self, p3d_name):  # pragma: no cover
        p3d_file = open(p3d_name, 'r')
        sline = p3d_file.readline().strip().split()
        assert len(sline) == 1, sline
        nblocks = int(sline[0])

        npoints = 0
        for i in range(nblocks):
            nx, ny, nz = p3d_file.readline().strip().split()
            nx = int(nx)
            ny = int(ny)
            nz = int(nz)
            self.block_shapes[i] = (nx, ny, nz)
            self.x[i] = zeros((nx * ny * nz), 'float32')
            self.y[i] = zeros((nx * ny * nz), 'float32')
            self.z[i] = zeros((nx * ny * nz), 'float32')
            npoints += nx * ny * nz

        nleft = npoints * 3
        iblock = 0
        nxyzi = 0
        ixyz = 0
        block = self.x[iblock]
        nxyz = len(block)
        nxyzi2 = None
        while nleft > 0:
            sline = p3d_file.readline().strip().split()
            floats = [float(s) for s in sline]
            nxyzi2 = nxyzi + len(floats)
            block[nxyzi : nxyzi2] = floats

            #print "sline = ", sline, nxyzi, nxyzi2, nxyz
            nxyzi = nxyzi2
            if nxyzi2 == nxyz:
                print("finished with block %i ixyz=%s" % (iblock, ixyz))
                #block = self.blocks[iblock]
                print("reshaping...", self.block_shapes[iblock])
                nleft -= nxyz
                blocki = block.reshape(self.block_shapes[iblock])
                if ixyz == 0:
                    self.x[iblock] = blocki
                elif ixyz == 1:
                    self.y[iblock] = blocki
                elif ixyz == 2:
                    self.z[iblock] = blocki
                else:
                    raise RuntimeError()

                # next block
                nxyzi = 0
                if ixyz == 0:
                    block = self.y[iblock]
                    nxyz = len(block)
                    ixyz = 1
                elif ixyz == 1:
                    block = self.z[iblock]
                    nxyz = len(block)
                    ixyz = 2
                elif ixyz == 2:
                    iblock += 1
                    if iblock == nblocks:
                        break
                    block = self.x[iblock]
                    nxyz = len(block)
                    ixyz = 0
                else:
                    raise NotImplementedError()
                print("iblock=%s icoeff=%s nleft=%s" %(iblock, ixyz, nleft))

            elif nxyzi2 > nxyz:
                asdf2

        print("finished with all blocks")

if __name__ == '__main__':  # pragma: no cover
    fname = 'HSCT-1.p3d'
    p3d = Plot3d()
    p3d.read_plot3d(fname)
