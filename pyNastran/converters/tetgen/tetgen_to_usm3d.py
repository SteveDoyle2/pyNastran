from pyNastran.converters.tetgen.tetgen import Tetgen
from pyNastran.converters.usm3d.usm3d_reader import write_usm3d_volume


def tetgen_to_usm3d(base):
    m = Tetgen()
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

    m2 = Tetgen()
    m2.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=2)
    ntris = m2.tris.shape[0]

    nboundary_nodes = m2.nodes.shape[0]
    if not isinstance(nboundary_nodes, int):
        msg = 'nboundary_nodes=%r; type=%s must be an integer' % (
            nboundary_nodes, type(nboundary_nodes))
        raise TypeError(msg)

    ntets = m.tets.shape[0]
    nnodes = m.nodes.shape[0]
    #assert isinstance(nnodes, int), type(nnodes)
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

def main():
    #base = 'tetgen_test_flipped.1'
    base = 'tetgen_test.1'
    tetgen_to_usm3d(base)

if __name__ == '__main__':  # pragma: no cover
    main()
