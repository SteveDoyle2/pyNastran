from pyNastran.converters.cart3d.cart3d import Cart3D
from pyNastran.converters.stl.stl import STL

# def cart3d_to_stl(cart3d, log=None, debug=False):
    # """
    # Converts a Cart3DReader object to STL format.

    # :param cart3d: a Cart3DReader object
    # :param log:    a logger object (or None)
    # :param debug:  True/False (used if log is not defined)

    # :returns stl: an STLReader object
    # """
    # normals = cart3d.get_normals()
    # stl = STLReader(log=log, debug=debug)
    # stl.nodes = cart3d.nodes
    # stl.elements = cart3d.elements
    # stl.write_stl(stl_filename)
    # return stl

def cart3d_to_stl_filename(cart3d_filename, stl_filename, is_binary=False, log=None, debug=False):
    """
    Converts a Cart3D file to STL format.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    stl_filename : str
        path to the output STL file
    is_binary : bool; default=False
        writes the stl in binary
    log : log
        a logger object (or None)
    debug : bool; default=False
        True/False (used if log is not defined)
    """
    cart3d = Cart3D(log=log, debug=debug)
    cart3d.read_cart3d(cart3d_filename)

    stl = STL()
    stl.nodes = cart3d.nodes - 1
    stl.elements = cart3d.elements - 1
    stl.write_stl(stl_filename, is_binary=is_binary)


if __name__ == '__main__':  # pragma: no cover
    bdf_filename = 'g278.bdf'
    cart3d_filename = 'g278.tri'
    nastran_to_cart3d(bdf_filename, cart3d_filename)
