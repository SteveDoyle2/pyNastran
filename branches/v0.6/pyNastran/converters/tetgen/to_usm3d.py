from struct import pack, unpack
from numpy import array, transpose, zeros, unique

from pyNastran.utils.log import get_logger

from pyNastran.converters.tetgen.tetgen_reader import TetgenReader
from pyNastran.converters.usm3d.usm3d_reader import write_usm3d_volume


def main():
    m = TetgenReader()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

    m2 = TetgenReader()
    base = 'tetgen_test_flipped.1'
    m2.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=2)
    ntris, four = m2.tris.shape
    #boundary_nodes = unique(m2.tris)
    #nboundary_nodes, = boundary_nodes.shape
    nboundary_nodes, three = m2.nodes.shape
    assert isinstance(nboundary_nodes, int), nboundary_nodes

    ntets, four = m.tets.shape
    nnodes, three = m.nodes.shape
    m.header = {
        'inew': -1,
        'nElements' : ntets,
        'nPoints'   : nnodes,
        'nBoundPts' : nboundary_nodes,
        'nViscPts'  : nnodes,
        'nViscElem' : ntets,
        'tc'        : 0.0, # d
    }
    write_usm3d_volume(m, base)

    #model = Usm3dReader()
    #basename = 'HSCT_inviscid'
    #basename = 'box'
    #model.read_usm3d(basename)
    #model.write_usm3d(basename + '_2')


if __name__ == '__main__':
    main()