import os

from pyNastran.utils import check_path
from pyNastran.converters.aflr.ugrid.ugrid_reader import UGRID


def ugrid3d_to_nastran(ugrid_filename, bdf_filename, include_shells=True, include_solids=True,
                       convert_pyram_to_penta=False, encoding=None,
                       size=16, is_double=False, log=None):
    """
    Converts a UGRID to a BDF.

    Parameters
    ----------
    ugrid_filename : str
        the input UGRID filename
    bdf_filename : str
        the output BDF filename
    include_shells : bool; default=True
        should the shells be written
    include_solids : bool; default=True
        should the solids be written
    convert_pyram_to_penta : bool; default=False
        False : NX Nastran
        True : MSC Nastran
    size : int; {8, 16}; default=16
        the bdf write precision
    is_double : bool; default=False
        the field precision to write
    log : logger; default=None
        a logger object

    Returns
    -------
    ugrid_model : UGRID()
        the ugrid model
    """
    model = UGRID(log=log, debug=False)
    check_path(ugrid_filename, 'ugrid_filename')
    model.read_ugrid(ugrid_filename)
    model.write_bdf(bdf_filename, include_shells=include_shells, include_solids=include_solids,
                    convert_pyram_to_penta=convert_pyram_to_penta,
                    encoding=encoding, size=size, is_double=is_double)
    return model

#def main():
    #import sys
    #if len(sys.argv) != 3:
        #msg = 'number of arguments must be 2; ugrid_filename, bdf_filename; nargs=%s; args=%s' % (
            #len(sys.argv[1:]), sys.argv[1:])
        #raise SyntaxError(msg)
    #ugrid_filename = sys.argv[1]
    #bdf_filename = sys.argv[2]
    #ugrid3d_to_nastran(ugrid_filename, bdf_filename,
                       #size=size, is_double=is_double)


#if __name__ == '__main__':   # pragma: no cover
    #main()
