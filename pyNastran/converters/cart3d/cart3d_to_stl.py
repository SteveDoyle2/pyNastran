"""
defines:
 - stl = cart3d_to_stl(cart3d, stl_filename=None, is_binary=False, log=None, debug=False)
 - stl = cart3d_to_stl_filename(cart3d_filename, stl_filename=None, is_binary=False, log=None, debug=False)
"""
from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
from pyNastran.converters.stl.stl import STL

def cart3d_to_stl(cart3d, stl_filename=None, is_binary=False, log=None, debug=False):
    """
    Converts a Cart3D object to STL format.

    Parameters
    ----------
    cart3d : Cart3D()
        a Cart3D object
    stl_filename : str; default=None
        path to the output STL file (or None to skip)
    log : log
        a logger object (or None)
    debug : bool; default=False
        True/False (used if log is not defined)

    Returns
    -------
    stl : STL()
        an STL object

    .. todo:: this seems to be broken...
    """
    #normals = cart3d.get_normals()
    stl = STL(log=log, debug=debug)
    stl.nodes = cart3d.nodes #- 1
    stl.elements = cart3d.elements #- 1
    assert stl.elements.min() == 0, stl.elements.min()
    if stl_filename:
        stl.write_stl(stl_filename, is_binary=is_binary)
    return stl

def cart3d_to_stl_filename(cart3d_filename, stl_filename=None, is_binary=False,
                           log=None, debug=False):
    """
    Converts a Cart3D file to STL format.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    stl_filename : str; default=None
        path to the output STL file (or None to skip
    is_binary : bool; default=False
        writes the stl in binary
    log : log
        a logger object (or None)
    debug : bool; default=False
        True/False (used if log is not defined)

    Returns
    -------
    stl : STL()
        an STL object
    """
    if isinstance(cart3d_filename, Cart3D):
        cart3d = cart3d_filename
    else:
        cart3d = read_cart3d(cart3d_filename, log=log, debug=debug)
    return cart3d_to_stl(cart3d, stl_filename, is_binary=is_binary,
                         log=cart3d.log, debug=debug)
