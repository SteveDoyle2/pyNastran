from six import iteritems
from pyNastran.bdf.bdf import BDF
from pyNastran.converters.ugrid.ugrid_reader import UGRID

from numpy import array, hstack

def nastran_to_ugrid(bdf_model, ugrid_filename_out=None, properties=None,
                     check_shells=True, check_solids=True):
    """
    set xref=False
    """
    # pids_to_inlcude = []
    # for pid, prop in iteritems(model.properties):
        # if prop.type == 'PSHELL':
            # pids_to_include.append(pid)
    if properties is not None:
        for pid, pid_new in iteritems(properties):
            bdf_model.properties[pid].pid = pid_new

    card_types = ['CQUAD4', 'CTRIA3', 'CTETRA', 'CHEXA', 'GRID', 'CPENTA', 'CPYRAM']
    out = bdf_model.get_card_ids_by_card_types(card_types)
    node_ids = out['GRID']
    ctria3 = out['CTRIA3']
    cquad4 = out['CQUAD4']

    ctetra = out['CTETRA']
    cpyram = out['CPYRAM']
    cpenta = out['CPENTA']
    chexa = out['CHEXA']

    nnodes = len(node_ids)
    ntris = len(ctria3)
    nquads = len(cquad4)
    nshells = ntris + nquads
    ntetra = len(ctetra)
    npyram = len(cpyram)
    npenta = len(cpenta)
    nhexa = len(chexa)
    nsolids = ntetra + npyram + npenta + nhexa
    if nnodes == 0:
        raise RuntimeError('nnodes=%i nshells=%i nsolids=%i' % (nnodes, nshells, nsolids))
    if nshells == 0 and check_shells:
        raise RuntimeError('nnodes=%i nshells=%i nsolids=%i' % (nnodes, nshells, nsolids))
    if nsolids == 0 and check_solids:
        raise RuntimeError('nnodes=%i nshells=%i nsolids=%i' % (nnodes, nshells, nsolids))

    nodes = bdf_model.nodes
    elements = bdf_model.elements
    xyz_cid0 = array([nodes[nid].xyz for nid in node_ids], dtype='float64')

    pids = []
    model = UGRID()
    model.nodes = xyz_cid0
    if ntris:
        model.tris = array([elements[eid].node_ids for eid in ctria3], dtype='int32')
        ptris = array([elements[eid].Pid() for eid in ctria3], dtype='int32')
        pids.append(ptris)
    if nquads:
        model.quads = array([elements[eid].node_ids for eid in cquad4], dtype='int32')
        pquads = array([elements[eid].Pid() for eid in cquad4], dtype='int32')
        pids.append(pquads)

    if check_shells:
        if len(pids) == 1:
            model.pids = pids[0]
        elif len(pids) == 2:
            model.pids = hstack(pids)
        else:
            raise RuntimeError(pids)

    if ntetra:
        model.tets = array([elements[eid].node_ids for eid in ctetra], dtype='int32')
    if npyram:
        model.penta5s = array([elements[eid].node_ids for eid in cpyram], dtype='int32')
    if nhexa:
        model.penta6s = array([elements[eid].node_ids for eid in cpenta], dtype='int32')
    if nhexa:
        model.hexas = array([elements[eid].node_ids for eid in chexa], dtype='int32')

    if ugrid_filename_out is not None:
        model.write_ugrid(ugrid_filename_out, check_shells=check_shells)
    return model

def merge_ugrids(a_model, b_model):
    """
    Merges two UGrid models

    TODO: not implemented
    """
    a_model

def main():
    """
    Converts a Nastran model to UGRID model and renumbers the properties.
    Also creates a fun3d.mapbc file.
    """
    properties_orig = {
        'bay' : ('viscous', [13, 15, 17, 17]),
        'inlet' : ('freestream', [1, 2, 3,
                                  18, 22, 25]),
        'left_wall' : ('reimann', [3, 26,
                                   6, 42,
                                   9, 54,
                                   12, 70]),
        'outlet' : ('reimann', [10, 11, 12, 62, 66]),
        'right_wall' : ('reimann', [1, 4, 7, 10,
                                    19, 37, 50, 63]),
        'top_wall' : ('viscous', [21, 24, 28,
                                  39, 41, 44,
                                  52, 56,
                                  65, 58, 72]),
    }
    bcname_to_bckey = {
        'viscous' : 4000,
        'reimann' : 5025,
        'freestream' : 5050,
    }
    pid = 1
    properties = {}
    with open('fun3d.mapbc', 'wb') as mapbc:
        for name, (bcname, pids) in sorted(iteritems(properties_orig)):
            for pidi in pids:
                properties[pidi] = pid
            bc = bcname_to_bckey[bcname]

            mapbc.write('%s %s # name=%s bcname=%s\n' % (bc, pid, name, bcname))
            pid += 1

    import sys
    from pyNastran.bdf.bdf_interface.dev_utils import bdf_renumber
    bdf_filename = sys.argv[1]
    ugrid_filename_out = sys.argv[2]
    # bdf_model = BDF(debug=False)
    # bdf_model.read_bdf(bdf_filename)

    bdf_filename_out = bdf_filename[:-4] + '.re.bdf'
    if 0:
        starting_id_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False,
                     starting_id_dict=starting_id_dict)
        bdf_model = BDF(debug=False)
        bdf_model.read_bdf(bdf_filename_out, xref=False)
    else:
        bdf_model = BDF(debug=False)
        bdf_model.read_bdf(bdf_filename_out, xref=False)


    # bdf_model = BDF(debug=False)
    # bdf_model.read_bdf(bdf_filename)
    nastran_to_ugrid(bdf_model, ugrid_filename_out) #, properties=properties


if __name__ == '__main__':
    main()
