from six.moves import zip
from numpy import unique

from pyNastran.converters.cart3d.cart3d_reader import Cart3DReader
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


def cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=None, debug=False):
    """
    Converts a Cart3D file to STL format.

    :param cart3d_filename: path to the input Cart3D file
    :param bdf_filename:    path to the output BDF file
    :param log:             a logger object (or None)
    :param debug:           True/False (used if log is not defined)
    """
    cart3d = Cart3DReader(log=log, debug=debug)
    (nodes, elements, regions, loads) = cart3d.read_cart3d(cart3d_filename)

    #bdf = BDF()
    #bdf.nodes = cart3d.nodes
    #bdf.elements = cart3d.elements
    #bdf.write_bdf(bdf_filename)
    #return
    f = open(bdf_filename, 'wb')
    f.write('CEND\n')
    f.write('BEGIN BULK\n')
    f.write('$Nodes\n')

    i = 0
    nid = 1
    cid = 0
    for node in nodes:
        card = print_card_16(['GRID', nid, cid] + list(node))
        f.write(card)
        nid += 1

    eid = 1
    f.write('$Elements\n')
    for (n1, n2, n3), pid in zip(elements, regions):
        card = print_card_8(['CTRIA3', eid, pid, n1, n2, n3])
        f.write(card)
        eid += 1

    t = 0.1
    E = 1e7
    nu = 0.3
    f.write('$Properties\n')
    for pid in unique(regions):
        mid = pid
        card = print_card_8(['PSHELL', pid, mid, t])
        f.write(card)
        card = print_card_8(['MAT1', mid, E, None, nu])
        f.write(card)
    f.write('ENDDATA\n')
    f.close()

def main():
    cart3d_filename = 'threePlugs.tri'
    bdf_filename = 'threePlugs.bdf'
    cart3d_to_nastran_filename(cart3d_filename, bdf_filename)

    if __name__ == '__main__':  # pragma: no cover
        main()
