from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
from pyNastran.converters.stl.stl_reader import STLReader

def cart3d_to_stl(cart3d, log=None, debug=False):
    """
    Converts a Cart3DReader object to STL format.
    
    :param cart3d: a Cart3DReader object
    :param log:    a logger object (or None)
    :param debug:  True/False (used if log is not defined)

    :returns stl: an STLReader object
    """
    normals = cart3d.normals()
    stl = STLReader(log=log, debug=debug)
    stl.nodes = nodes
    stl.elements = elements
    stl.write_stl(stl_filename)
    return stl

def cart3d_to_stl_filename(cart3d_filename, stl_filename, log=None, debug=False):
    """
    Converts a Cart3D file to STL format.
    
    :param cart3d_filename: path to the input Cart3D file
    :param stl_filename:    path to the output STL file
    :param log:             a logger object (or None)
    :param debug:           True/False (used if log is not defined)
    """
    cart3d = Cart3DReader(log=log, debug=debug)
    cart3d.read_cart3d(cart3d_filename)

    stl = STLReader()
    stl.nodes = cart3d.nodes
    stl.elements = cart3d.elements
    stl.write_stl(stl_filename)
    return
    if 0:
        normals = model.normals()

        f = open(stl_filename, 'wb')
        f.write('solid cart3d_model\n')

        i = 0
        for (n1, n2, n3) in model.elements:
                #f.write('solid cart3d_model\n')
                n = normals[i]
                f.write('loop\n')
                f.write('  facet normal %f %f %f\n' % (n[0], n[1], n[2]))
                f.write('    vertex %f %f %f\n' % (n1[0], n1[1], n1[2]))
                f.write('    vertex %f %f %f\n' % (n1[0], n1[1], n1[2]))
                f.write('    vertex %f %f %f\n' % (n1[0], n1[1], n1[2]))
                f.write('  endfacet')
                f.write('endloop\n')
                i += 1
        f.close()

if __name__ == '__main__':
    bdf_filename = 'g278.bdf'
    cart3d_filename = 'g278.tri'
    nastran_to_cart3d(bdf_filename, cart3d_filename)