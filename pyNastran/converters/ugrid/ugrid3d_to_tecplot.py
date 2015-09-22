from six import iteritems
import os

from numpy import zeros, array
from pyNastran.converters.ugrid.ugrid_reader import UGRID


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

def main():
    import sys
    assert len(sys.argv) == 3, 'number of arguments must be 2; ugrid_filename, tecplot_filename; nargs=%s; args=%s' % (len(sys.argv[1:]), sys.argv[1:])
    ugrid_filename = sys.argv[1]
    tecplot_filename = sys.argv[2]
    ugrid3d_to_tecplot_filename(ugrid_filename, tecplot_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
