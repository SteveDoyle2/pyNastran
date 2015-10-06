from six import iteritems
import os

from numpy import zeros, array
from pyNastran.converters.ugrid.ugrid_reader import UGRID
from pyNastran.converters.tecplot.tecplot import Tecplot


def ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename, log=None):
    """
    Converts a UGRID to a Tecplot ASCII file.

    Parameters
    ----------
    ugrid_filename : str
        the input UGRID filename
    bdf_filename : str
        the output Tecplot filename
    """
    model = UGRID(log=log)
    assert os.path.exists(ugrid_filename), '%r doesnt exist' % ugrid_filename
    model.read_ugrid(ugrid_filename)
    model.write_tecplot(tecplot_filename)



#def write_tecplot(self, tecplot_filename):
    #self.check_hanging_nodes()
    #tecplot = ugrid_to_tecplot(self)
    #tecplot.write_tecplot(tecplot_filename, adjust_nids=True)  # is adjust correct???
    #tecplot.results = array([], dtype='float32')

def ugrid_to_tecplot(ugrid_model):
    nnodes = len(ugrid_model.nodes)
    nodes = zeros((nnodes, 3), dtype='float64')
    ugrid_model.check_hanging_nodes()
    elements = []

    ntets = len(ugrid_model.tets)
    non_tets = len(ugrid_model.penta5s) + len(ugrid_model.penta6s) + len(ugrid_model.hexas)
    assert ntets + non_tets > 0, 'nsolids=%s' % (ntets + non_tets)

    tecplot = Tecplot()
    tecplot.xyz = ugrid_model.nodes

    if ntets and non_tets == 0:
        elements = ugrid_model.tets
        tecplot.tet_elements = elements - 1
    elif non_tets:
        for element in ugrid_model.tets:
            n1, n2, n3, n4 = element
            elements.append([n1, n2, n3, n4,
                             n4, n4, n4, n4])
        for element in ugrid_model.penta5s:
            n1, n2, n3, n4, n5 = element
            elements.append([n1, n2, n3, n4,
                             n5, n5, n5, n5])
        for element in ugrid_model.penta6s:
            n1, n2, n3, n4, n5, n6 = element
            elements.append([n1, n2, n3, n4,
                             n5, n6, n6, n6])
        for element in ugrid_model.hexas:
            n1, n2, n3, n4, n5, n6, n7, n8 = element
            elements.append([n1, n2, n3, n4,
                             n5, n6, n7, n8])
        elements = array(elements, dtype='int32') - 1
        tecplot.hexa_elements = elements
    else:
        raise RuntimeError()
    return tecplot


def main():
    import sys
    assert len(sys.argv) == 3, 'number of arguments must be 2; ugrid_filename, tecplot_filename; nargs=%s; args=%s' % (len(sys.argv[1:]), sys.argv[1:])
    ugrid_filename = sys.argv[1]
    tecplot_filename = sys.argv[2]
    ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
