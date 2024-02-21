from __future__ import annotations
from typing import TYPE_CHECKING
from numpy import unique

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.cart3d.cart3d import Cart3D, read_cart3d
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    CASE = dict[int, float]
    LOADS_DICT = dict[int, CASE]

def cart3d_to_nastran_model(cart3d_filename: str,
                            loads_map: Optional[LOADS_DICT]=None,
                            log: Optional[SimpleLogger]=None,
                            debug: bool=False) -> BDF:
    """
    Converts a Cart3D file to Nastran format and returns a BDF() object.

    Parameters
    ----------
    cart3d_filename : str
        path to the input Cart3D file
    loads_dict: dict[load_id, case]
        case : dict[eid, pressure]
        write the loads
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
    model = BDF(log=log, debug=debug)
    for xyz in nodes:
        model.add_grid(nid, xyz)
        nid += 1

    eid = 1
    for nids, pid in zip(elements, regions):
        model.add_ctria3(eid, pid, nids)
        #print(model.elements[eid])
        eid += 1

    t = 0.1
    E = 1e7
    G = None
    nu = 0.3
    for pid in unique(regions):
        mid = pid
        model.add_pshell(pid, mid1=mid, t=t)
        model.add_mat1(mid, E, G, nu)

    if loads_map:
        for name, load_id in loads_map.items():
            loads = cart3d.loads[name]
            comment = name
            for inid, Cp in enumerate(loads):
                loadi = model.add_sload(load_id, [inid+1], [Cp], comment=comment)
                comment = ''

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


        load_id = 1
        if 'Cp' in cart3d.loads:
            Cp = cart3d.loads['Cp']
            #+--------+-----+------+------+------+------+------+------+------+
            #|    1   |   2 |  3   |  4   |   5  |   6  |   7  |   8  |   9  |
            #+========+=====+======+======+======+=============+======+======+
            #| PLOAD2 | SID |  P   | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 |
            #+--------+-----+------+------+------+------+------+------+------+

            #from collections import defaultdict
            #nid_to_eids = defaultdict(list)

            assert len(Cp) == len(nodes)
            inode1 = elements[:, 0] - 1
            inode2 = elements[:, 1] - 1
            inode3 = elements[:, 2] - 1
            cp1 = Cp[inode1]
            cp2 = Cp[inode2]
            cp3 = Cp[inode3]
            cp_avg = (cp1 + cp2 + cp3) / 3.

            for eid, cp_avgi in enumerate(cp_avg):
                card = print_card_8(['PLOAD2', load_id, cp_avgi, eid + 1])  # +1 b/c it's 0-based
                bdf_file.write(card)

        bdf_file.write('ENDDATA\n')

if __name__ == '__main__':  # pragma: no cover
    import sys
    cart3d_filename = sys.argv[1]
    print(cart3d_filename)
    bdf_filename = 'spike.bdf'
    cart3d_to_nastran_filename(cart3d_filename, bdf_filename, log=None, debug=False)
