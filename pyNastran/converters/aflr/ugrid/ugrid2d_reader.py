"""
Defines the following classes:
    - UGRID2D_Reader
"""
from numpy import array
from cpylog import get_logger2


class UGRID2D_Reader:
    """
    Interface to the AFLR UGrid2D format.
    """
    def __init__(self, log=None, debug=None):
        self.log = get_logger2(log=log, debug=debug)
        self.nodes = None
        self.tris = None
        self.quads = None

    def read_ugrid(self, ugrid_filename):
        """
        Reads a ugrid2d file of the form::

          #(nnodes, ntrias, nquads), ntets, npyram5, npenta6, nhexas8s
          '5 1 1   0 0 0 0\n'
          # nodes
          '0. 0. 0.\n'
          '1. 0. 0.\n'
          '1. 1. 0.\n'
          '0. 1. 0.\n'
          '0. 2. 0.\n'
          # tris
          '3 4 5\n'
          # quads
          '1 2 3 4\n'

        .. note:: comment lines should not be included and exist for reference
        """
        with open(ugrid_filename, 'r') as ugrid_file:
            ugrid_str = ugrid_file.read()
            data = ugrid_str.split()

        try:
            nnodes, ntrias, nquads, ntets, npyram5, npenta6, nhexas8s = (
                int(val) for val in data[:7])
        except ValueError:
            self.log.error('data[:7] = %s' % data[:7])
            raise
        i = 7

        self.log.debug('nnodes=%s ntrias=%s nquads=%s ntets4=%s '
                       'npyram5=%s npenta6=%s nhexas8s=%s' % (
                           nnodes, ntrias, nquads, ntets,
                           npyram5, npenta6, nhexas8s))


        #nodes = zeros(nnodes * 3, dtype=ndarray_float)
        #tris  = zeros(ntris * 3, dtype='int32')
        #quads = zeros(nquads * 4, dtype='int32')
        #pids = zeros(npids, dtype='int32')

        #tets = zeros(ntets * 4, dtype='int32')
        #penta5s = zeros(npenta5s * 5, dtype='int32')
        #penta6s = zeros(npenta6s * 6, dtype='int32')
        #hexas = zeros(nhexas * 8, dtype='int32')

        # nodes
        iend = i + nnodes * 3
        nodes = array(data[i:iend], dtype='float64')
        if len(nodes) != (iend - i):
            raise RuntimeError('len(xyz)=%s for nnodes=%s (expected=%s)' % (
                len(nodes), (iend - i)//3, iend - i))
        nodes = nodes.reshape((nnodes, 3))
        self.log.debug(nodes[0, :])
        self.log.debug(nodes[nnodes-1, :])
        self.log.debug(nodes[-1, :])
        assert nodes[:, 2].max() == 0.0
        assert nodes[:, 2].min() == 0.0
        i = iend

        # tris
        iend = i + ntrias * 3
        tris = array(data[i:iend], dtype='int32')
        if len(tris) != (iend - i):
            raise RuntimeError('len(tri_nodes)=%s for ntris=%s (expected=%s)' % (
                len(tris), (iend - i)//3, iend - i))
        tris = tris.reshape((ntrias, 3))
        #print(nodes[0, :])
        #print(nodes[nnodes-1, :])
        #print(nodes[-1, :])
        #assert nodes[:, 2].max() == 0.0
        #assert nodes[:, 2].min() == 0.0
        i = iend

        iend = i + nquads * 4
        quads = array(data[i:iend], dtype='int32')
        if len(quads) != (iend - i):
            raise RuntimeError('len(quad_nodes)=%s for nquads=%s (expected=%s)' % (
                len(quads), (iend - i)//3, iend - i))
        quads = quads.reshape((nquads, 4))
        #print(nodes[0, :])
        #print(nodes[nnodes-1, :])
        #print(nodes[-1, :])
        #assert nodes[:, 2].max() == 0.0
        #assert nodes[:, 2].min() == 0.0
        i = iend

        self.nodes = nodes
        self.tris = tris - 1
        self.quads = quads - 1
