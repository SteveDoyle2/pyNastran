import numpy as np

from pyNastran.converters.stl.stl import STL, read_stl
from pyNastran.converters.cart3d.cart3d import Cart3D


def stl_to_cart3d(stl_filename, cart3d_filename=None, log=None, debug=False,
                  is_binary=False, float_fmt='%6.7f'):
    """
    Converts a Cart3D object to STL format.

    Parameters
    ----------
    stl_filename : str / STL()
        str : path to the input STL file
        STL() : an STL object
    cart3d_filename : str; default=None
        str : path to the output Cart3D file (or None to skip)
    log : log
        a logger object (or None)
    debug : bool; default=False
        True/False (used if log is not defined)
    is_binary : bool; default=False
        should the cart3d file be binary
    float_fmt : str; default='6.7f'
        the cart3d node precision

    Returns
    -------
    stl : STL()
        an STL object

    """
    if isinstance(stl_filename, str):
        stl = read_stl(stl_filename, log=log, debug=debug)
    elif isinstance(stl_filename, STL):
        stl = stl_filename
    else:
        raise TypeError('stl_filename must be a string or STL; type=%s' % type(stl_filename))

    cart3d = Cart3D(log=log, debug=debug)
    cart3d.nodes = stl.nodes
    cart3d.elements = stl.elements
    nelements = len(stl.elements)
    cart3d.regions = np.zeros(nelements, dtype='int32')
    if cart3d_filename:
        cart3d.write_cart3d(cart3d_filename, is_binary=is_binary, float_fmt=float_fmt)
    return cart3d
