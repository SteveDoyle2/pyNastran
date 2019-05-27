import numpy as np
from cpylog import get_logger2
from pyNastran.utils import check_path


def read_fgrid(fgrid_filename, unused_dimension_flag, log=None, debug=False):
    """loads a *.fgrid file"""
    model = FGridReader(log=log, debug=debug)
    model.read_fgrid(fgrid_filename, unused_dimension_flag=3)
    return model


class FGridReader:
    """FGRID interface class"""
    def __init__(self, log=None, debug=False):
        """
        Initializes the FGridReader object

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
        self.log = get_logger2(log, debug=debug)
        self.debug = debug
        self.nodes = None
        self.bcs = None
        self.tris = None
        self.tets = None

    def read_fgrid(self, fgrid_filename, unused_dimension_flag=3):
        """extracts the nodes, tris, bcs, tets"""
        check_path(fgrid_filename, 'fgrid_filename')
        with open(fgrid_filename, 'r') as fgrid:
            nnodes, ntris, ntets = fgrid.readline().split()
            nnodes = int(nnodes)
            ntris = int(ntris)
            ntets = int(ntets)

            self.log.info('nnodes=%s ntris=%s ntets=%s' % (nnodes, ntris, ntets))
            assert nnodes > 0, nnodes
            #inode = 0
            # I think this goes xxx, yyy, zzz
            # instead of x, y, z
            #            x, y, z
            xyz = []
            nfloats = 0
            while nfloats < nnodes * 3:
                sline = fgrid.readline().split()
                nfloatsi = len(sline)
                nfloats += nfloatsi
                xyz.extend(sline)
            nodes = np.array(xyz, dtype='float32')
            assert nfloats == nnodes * 3, 'nfloats=%s nnodes*3=%s' % (nfloats, nnodes * 3)
            assert nodes.max() > 0, nodes.max()

            # we want a contiguous array
            self.nodes = nodes.reshape((3, nnodes)).T.ravel().reshape(nnodes, 3)
            #---------------------------------------------------------------------

            tris = []
            nints = 0
            while nints < ntris * 3:
                sline = fgrid.readline().split()
                nintsi = len(sline)
                nints += nintsi
                tris.extend(sline)
            if tris:
                self.tris = np.array(tris, dtype='int32').reshape(ntris, 3)

            #---------------------------------------------------------------------
            nints = 0
            bcs = []
            while nints < ntris:
                sline = fgrid.readline().split()
                nintsi = len(sline)
                nints += nintsi
                bcs.extend(sline)
            if bcs:
                self.bcs = np.array(bcs, dtype='int32')

            #---------------------------------------------------------------------
            nints = 0
            tets = []
            while nints < ntets * 4:
                sline = fgrid.readline().split()
                nintsi = len(sline)
                nints += nintsi
                tets.extend(sline)
            if bcs:
                self.tets = np.array(tets, dtype='int32').reshape((ntets, 4))
                #print(self.tets.shape)
