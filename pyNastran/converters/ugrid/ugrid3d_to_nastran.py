from __future__ import print_function
import os

from pyNastran.converters.ugrid.ugrid_reader import UGRID


def ugrid3d_to_nastran(ugrid_filename, bdf_filename, include_shells=True, include_solids=True,
                       convert_pyram_to_penta=False, encoding=None, size=16, is_double=False):
    """
    Converts a UGRID to a BDF.

    Parameters
    ----------
    ugrid_filename : str
        the input UGRID filename
    bdf_filename : str
        the output BDF filename
    size : int; {8, 16}; default=16
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    """
    model = UGRID(log=None, debug=False)
    assert os.path.exists(ugrid_filename), '%r doesnt exist' % ugrid_filename
    model.read_ugrid(ugrid_filename)
    model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids,
                    convert_pyram_to_penta=convert_pyram_to_penta,
                    encoding=encoding, size=size, is_double=is_double)


def main():
    import sys
    assert len(sys.argv) == 3, 'number of arguments must be 2; ugrid_filename, bdf_filename; nargs=%s; args=%s' % (len(sys.argv[1:]), sys.argv[1:])
    ugrid_filename = sys.argv[1]
    bdf_filename = sys.argv[2]
    ugrid3d_to_nastran(ugrid_filename, bdf_filename,
                       size=size, is_double=is_double)


if __name__ == '__main__':
    main()
