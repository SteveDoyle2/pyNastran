from numpy import unique

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d import Cart3D, read_cart3d
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


def cart3d_to_nastran_model(cart3d_filename, log=None, debug=False):
    """
    Converts a Cart3D file to Nastran format and returns a BDF() object.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    log : log / None
        log : a logger object
        None : a log will be defined
    debug : bool
        True/False (used if log is not defined)

    Returns
    -------
    bdf_model : BDF
        BDF() model object

    """
    if isinstance(cart3d_filename, Cart3D):
        cart3d = cart3d_filename
    else:
        cart3d = read_cart3d(cart3d_filename, log=log, debug=debug, result_names=None)
    nodes = cart3d.nodes
    elements = cart3d.elements + 1
    regions = cart3d.regions

    if regions.min() == 0:
        # bit of a hack to take an invalid cart3d model and make it
        # work in Nastran, which requires property_ids > 0
        regions += 1

    i = 0
    nid = 1
    cid = 0
    model = BDF(log=log, debug=debug)
    for node in nodes:
        card = ['GRID', nid, cid] + list(node)
        model.add_card(card, 'GRID', is_list=True)
        nid += 1

    eid = 1
    for (n1, n2, n3), pid in zip(elements, regions):
        card = ['CTRIA3', eid, pid, n1, n2, n3]
        model.add_card(card, 'CTRIA3', is_list=True)
        #print(model.elements[eid])
        eid += 1

    t = 0.1
    E = 1e7
    nu = 0.3
    for pid in unique(regions):
        mid = pid
        card = ['PSHELL', pid, mid, t]
        model.add_card(card, 'PSHELL', is_list=True)
        card = ['MAT1', mid, E, None, nu]
        model.add_card(card, 'MAT1', is_list=True)
    model.pop_parse_errors()
    return model


def cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=None, debug=False):
    """
    Converts a Cart3D file to Nastran format.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    bdf_filename : str
        path to the output BDF file
    log : log / None
        log : a logger object
        None : a log will be defined
    debug : bool
        True/False (used if log is not defined)

    Examples
    --------
    >>> cart3d_filename = 'threePlugs.tri'
    >>> bdf_filename = 'threePlugs.bdf'
    >>> cart3d_to_nastran_filename(cart3d_filename, bdf_filename)
    """
    if isinstance(cart3d_filename, Cart3D):
        cart3d = cart3d_filename
    else:
        cart3d = read_cart3d(cart3d_filename, log=log, debug=debug, result_names=None)

    nodes = cart3d.nodes
    elements = cart3d.elements + 1
    regions = cart3d.regions

    if regions.min() == 0:
        # bit of a hack to take an invalid cart3d model and make it
        # work in Nastran, which requires property_ids > 0
        regions += 1

    #bdf = BDF()
    #bdf.nodes = cart3d.nodes
    #bdf.elements = cart3d.elements
    #bdf.write_bdf(bdf_filename)
    #return
    with open(bdf_filename, 'w') as bdf_file:
        bdf_file.write('CEND\n')
        bdf_file.write('BEGIN BULK\n')
        bdf_file.write('$Nodes\n')

        i = 0
        nid = 1
        cid = 0
        for node in nodes:
            card = print_card_16(['GRID', nid, cid] + list(node))
            bdf_file.write(card)
            nid += 1

        eid = 1
        bdf_file.write('$Elements\n')
        assert 0 not in elements
        for (n1, n2, n3), pid in zip(elements, regions):
            card = print_card_8(['CTRIA3', eid, pid, n1, n2, n3])
            bdf_file.write(card)
            eid += 1

        t = 0.1
        E = 1e7
        nu = 0.3
        bdf_file.write('$Properties\n')
        for pid in unique(regions):
            mid = pid
            card = print_card_8(['PSHELL', pid, mid, t])
            bdf_file.write(card)
            card = print_card_8(['MAT1', mid, E, None, nu])
            bdf_file.write(card)
        bdf_file.write('ENDDATA\n')
