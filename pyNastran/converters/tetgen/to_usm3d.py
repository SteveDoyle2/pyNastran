from struct import pack, unpack

from pyNastran.utils.log import get_logger
from pyNastran.converters.dev.tetgen.tetgen import Tetgen
from pyNastran.converters.usm3d.usm3d_reader import write_usm3d_volume


def main():
    m = Tetgen()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

    m2 = Tetgen()
    base = 'tetgen_test_flipped.1'
    m2.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=2)
    ntris = m2.tris.shape[0]
    #boundary_nodes = unique(m2.tris)
    #nboundary_nodes, = boundary_nodes.shape
    nboundary_nodes = m2.nodes.shape[0]
    assert isinstance(nboundary_nodes, int), nboundary_nodes

    ntets = m.tets.shape[0]
    nnodes = m.nodes[0]
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

    #model = Usm3d()
    #basename = 'HSCT_inviscid'
    #basename = 'box'
    #model.read_usm3d(basename)
    #model.write_usm3d(basename + '_2')


if __name__ == '__main__':  # pragma: no cover
    main()