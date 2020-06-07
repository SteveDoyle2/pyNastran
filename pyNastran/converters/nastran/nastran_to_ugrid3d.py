"""
defines:
 - merge_ugrid3d_and_bdf_to_ugrid3d_filename(ugrid_filename, bdf_filename, ugrid_filename_out,
                                             pshell_pids_to_remove,
                                             update_equivalence=True, tol=0.01)
 - equivalence_ugrid3d_and_bdf_to_bdf(ugrid_filename, bdf_filename,
                                       pshell_pids_to_remove,
                                       tol=0.01, renumber=True)

"""
import os
from struct import Struct
from numpy import array, unique #, hstack

from cpylog import get_logger2
from pyNastran.utils import check_path
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber

from pyNastran.converters.aflr.ugrid.ugrid_reader import (
    UGRID, determine_dytpe_nfloat_endian_from_ugrid_filename)


def merge_ugrid3d_and_bdf_to_ugrid3d_filename(ugrid_filename, bdf_filename, ugrid_filename_out,
                                              pshell_pids_to_remove,
                                              update_equivalence=True, tol=0.01, log=None):
    """
    assumes cid=0

        Parameters
    ----------
    ugrid_filename : str
        the AFLR3/UGrid3d filename
    bdf_filename : str
        the BDF filename
    ugrid_filename_out : str
        the output AFLR3/UGrid3d filename
    pshell_pids_to_remove : List[int, ...]
        ???
    tol : float; default=0.01
        the equivalence tolerance
    update_equivalence : bool; default=True
        calls ``equivalence_ugrid3d_and_bdf_to_bdf`` to equivalence nodes
    """
    #base, ext = os.path.splitext(ugrid_filename_out)
    #bdf_filename = base + '.bdf'
    log = get_logger2(log, debug=True)
    log.debug(f'merge_ugrid3d_and_bdf_to_ugrid3d_filename - bdf_filename = {bdf_filename}')
    log.debug(f'merge_ugrid3d_and_bdf_to_ugrid3d_filename - ugrid_filename = {ugrid_filename}')

    if update_equivalence:
        bdf_filename2 = equivalence_ugrid3d_and_bdf_to_bdf(ugrid_filename, bdf_filename,
                                                           pshell_pids_to_remove,
                                                           tol, renumber=True, log=log)
    else:
        base = os.path.splitext(bdf_filename)[0]
        #bdf_merged_filename = base + '_merged.bdf'
        #bdf_equivalence_filename = base + '_equivalence.bdf'
        bdf_renumber_filename = base + '_renumber.bdf'
        bdf_filename2 = bdf_renumber_filename

    log.debug(f'**** bdf_filename2 = {bdf_filename2}')
    model = read_bdf(bdf_filename2, xref=False, log=log)

    outi = determine_dytpe_nfloat_endian_from_ugrid_filename(ugrid_filename)
    ndarray_float, float_fmt, nfloat, endian, ugrid_filename = outi
    #ndarray_float, float_fmt, nfloat, endian, ugrid_filename

    # more for documentation than anything else
    assert ndarray_float in ['float32', 'float64'], ndarray_float
    assert float_fmt in ['f', 'd'], float_fmt
    assert nfloat in [4, 8], nfloat
    assert endian in ['<', '>'], endian

    cards_to_get = ['GRID', 'CTRIA3', 'CQUAD4', 'PSHELL', 'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA']
    out = model.get_card_ids_by_card_types(cards_to_get)
    nids = out['GRID']
    nnodes = len(nids)

    ntris = len(out['CTRIA3'])
    nquads = len(out['CQUAD4'])

    ntets = len(out['CTETRA'])
    npyramids = len(out['CPYRAM'])
    npentas = len(out['CPENTA'])
    nhexas = len(out['CHEXA'])
    nshells = ntris + nquads
    nsolids = ntets + npyramids + npentas + nhexas
    if nshells == 0:
        raise RuntimeError('the UGRID model does not have any boundary condition surfaces...')
    assert nsolids > 0, nsolids

    #nodes = zeros((nnodes, 3), dtype=ndarray_float)

    #tris = zeros((ntris, 3), dtype='int32')
    #quads = zeros((nquads, 4), dtype='int32')

    #tets = zeros((ntetras, 4), dtype='int32')
    #pyramids = zeros((npyramids, 5), dtype='int32')
    #pentas = zeros((npyramids, 6), dtype='int32')
    #hexas = zeros((nhexas, 6), dtype='int32')

    xyz = array([model.nodes[nid].xyz for nid in sorted(nids)],
                dtype=ndarray_float)

    # get the pshells
    #pshells = out['PSHELL']
    # TODO: need to think about property IDs
    tris = out['CTRIA3']
    quads = out['CQUAD4']
    eids = tris + quads
    pids = [model.elements[eid].pid for eid in eids]

    #tris_shrink = [eid for eid, pid in zip(eids[:ntris], pids[:ntris])
                   #if pid in pshell_pids_to_save]
    #quads_shrink = [eid for eid, pid in zip(eids[ntris:], pids[ntris:])
                    #if pid in pshell_pids_to_save]
    #ntris = len(tris_shrink)
    #nquads = len(quads_shrink)
    nshells = nquads + ntris
    npids = len(pids)
    if not nshells == npids:
        raise RuntimeError('nshells=%s npids=%s; must be the same' % (nshells, npids))

    #pids_shrink = [pidi for pidi in pids
                   #if pidi in pshell_pids_to_save]
    pids_shrink = pids

    with open(ugrid_filename_out, 'wb') as f_ugrid:
        #element_ids = hstack([
            #out['CTRIA3'], out['CQUAD4'],
            #out['CTETRA'] + out['CPYRAM'], out['CPENTA'], out['CHEXA']
        #])
        #fmt = '%ii' % (nsolids * nnodes)
        structi = Struct(endian + '7i')
        f_ugrid.write(structi.pack(nnodes, ntris, nquads, ntets, npyramids, npentas, nhexas))

        # %3f or %3d
        fmt = endian + '%i%s' % (nnodes * 3, float_fmt) # len(x,y,z) = 3
        structi = Struct(fmt)
        f_ugrid.write(structi.pack(*xyz.ravel()))

        for card_type in cards_to_get[1:]:  # drop the GRIDs & PSHELLs
            if card_type == 'PSHELL':
                assert len(pids) > 0, 'pids=%s' % pids
                #print('writing %s' % card_type)

                # %10i
                fmt = endian + '%ii' % (nshells)
                structi = Struct(fmt)
                pids = pids_shrink
                f_ugrid.write(structi.pack(*pids))
            elif card_type in ['CTRIA3', 'CQUAD4'] and 0:
                if card_type == 'CTRIA3':
                    eids = tris
                elif card_type == 'CQUAD4':
                    eids = quads

                # if there are cards
                if len(eids):
                    #print('writing %s' % card_type)
                    nelements = len(eids)
                    eid0 = eids[0]

                    # get the 0th card so we can size the formatter
                    element0 = model.elements[eid0]
                    nnodes_per_element = len(element0.nodes)

                    node_ids = array([model.elements[eid].node_ids for eid in sorted(eids)],
                                     dtype='int32')

                    # '%8i'
                    fmt = endian + '%ii' % (nelements * nnodes_per_element)
                    structi = Struct(fmt)
                    f_ugrid.write(structi.pack(*node_ids.ravel()))
            else:
                eids = out[card_type]

                # if there are cards
                if len(eids):
                    #print('writing %s' % card_type)
                    nelements = len(eids)
                    eid0 = eids[0]

                    # get the 0th card so we can size the formatter
                    element0 = model.elements[eid0]
                    nnodes_per_element = len(element0.nodes)

                    node_ids = array([model.elements[eid].node_ids for eid in sorted(eids)],
                                     dtype='int32')
                    # '%8i'
                    fmt = endian + '%ii' % (nelements * nnodes_per_element)
                    structi = Struct(fmt)
                    f_ugrid.write(structi.pack(*node_ids.ravel()))



def equivalence_ugrid3d_and_bdf_to_bdf(ugrid_filename, bdf_filename,
                                       pshell_pids_to_remove,
                                       tol=0.01, renumber=True, log=None):
    """
    Merges a UGRID3D (*.ugrid) with a BDF and exports a BDF that is
    equivalenced and renumbered.

    Parameters
    ----------
    ugrid_filename : str
        the AFLR3/UGrid3d filename
    bdf_filename : str
        the BDF filename
    pshell_pids_to_remove : List[int, ...]
    tol : float; default=0.01
        the equivalence tolerance
    renumber : bool; default=True
        calls ``bdf_renumber`` to renumber the output BDF model

    Returns
    -------
    out_bdf_filename : str
        the output BDF filename
    """
    log = get_logger2(log, debug=True)
    log.info(f'equivalence_ugrid3d_and_bdf_to_bdf - bdf_filename={bdf_filename}')
    log.info(f'equivalence_ugrid3d_and_bdf_to_bdf - ugrid_filename={ugrid_filename}')
    check_path(ugrid_filename, 'ugrid_filename')

    base = os.path.splitext(bdf_filename)[0]
    #bdf_merged_filename = base + '_merged.bdf'
    bdf_equivalence_filename = base + '_equivalence.bdf'
    bdf_renumber_filename = base + '_renumber.bdf'

    update_merge = True
    if update_merge:
        ugrid_model = UGRID(log=log, debug=False)
        ugrid_model.read_ugrid(ugrid_filename)

        bdf_model = read_bdf(bdf_filename, xref=False, log=log)
        #bdf_model.write_bdf(bdf_merged_filename, interspersed=False, enddata=False)

        tol = 0.01
        nid0 = max(bdf_model.nodes) + 1  # new node ids start at max+1
        nid_offset = nid0 - 1            # node_ids are 1-based, so we must offset them
        eid = max(bdf_model.elements) + 1

        cp = None
        for nid, node in enumerate(ugrid_model.nodes):
            #assert len(node) == 3, node
            card = ['GRID', nid + nid0, cp] + list(node)
            bdf_model.add_card(card, 'GRID', is_list=True)
            #f.write(print_card_double(card))

        pid_solid = 100
        mid = 1

        pids = unique(ugrid_model.pids)
        for pidi in pids:
            if pidi not in pshell_pids_to_remove:
                card = ['PSHELL', pidi, mid, 0.1]
                bdf_model.add_card(card, 'PSHELL', is_list=True)

        card = ['PSOLID', pid_solid, mid]
        bdf_model.add_card(card, 'PSOLID', is_list=True)

        card = ['MAT1', mid, 3.0e7, None, 0.3]
        bdf_model.add_card(card, 'MAT1', is_list=True)

        shells = [
            ('CQUAD4', ugrid_model.quads),
            ('CTRIA3', ugrid_model.tris),
        ]
        for card_type, card_nodes in shells:
            if card_nodes.shape[0]:
                for pid, nodes in zip(ugrid_model.pids, card_nodes + nid_offset):
                    if pid not in pshell_pids_to_remove:
                        card = [card_type, eid, pid, ] + list(nodes)
                        bdf_model.add_card(card, card_type, is_list=True)
                        eid += 1

        solids = [
            ('CTETRA', ugrid_model.tets),
            ('CPYRAM', ugrid_model.penta5s),
            ('CPENTA', ugrid_model.penta6s),
            ('CHEXA', ugrid_model.hexas),
        ]
        for card_type, card_nodes in solids:
            if card_nodes.shape[0]:
                for nodes in card_nodes + nid_offset:
                    card = [card_type, eid, pid_solid, ] + list(nodes)
                    bdf_model.add_card(card, card_type, is_list=True)
                    eid += 1

        # tol = min_edge_length / 2.0
        # TODO:  remove this...
        bdf_model.write_bdf('model_join.bdf', interspersed=False)
        bdf_model.cross_reference()
        bdf_equivalence_nodes(bdf_model, bdf_equivalence_filename, tol,
                              renumber_nodes=False, neq_max=10, xref=False, log=log)

    if renumber:
        starting_ids_dict = {
            'cid' : 1,
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(bdf_equivalence_filename, bdf_renumber_filename, size=16, is_double=False,
                     starting_id_dict=starting_ids_dict, log=log)
        #os.remove(bdf_equivalence_filename)
        out_bdf_filename = bdf_renumber_filename
    else:
        out_bdf_filename = bdf_equivalence_filename

    #os.remove(bdf_merged_filename)
    #os.remove(bdf_renumber_filename)
    os.remove('model_join.bdf')
    return out_bdf_filename

    #bdf_model.write_bdf(bdf_renumber_filename, interspersed=False)


def main():  # pragma: no cover
    """demo problem"""
    PKG_PATH = ''
    model_dir = os.path.join(PKG_PATH, 'aflr_work', 'bay')
    ugrid_filename = os.path.join(model_dir, 'model.b8.ugrid')
    bdf_filename = os.path.join(model_dir, 'solid.bdf')
    ugrid_filename_out = os.path.join(model_dir, 'model_merged.b8.ugrid')
    pshell_pids_to_remove = []
    tol = 0.01
    merge_ugrid3d_and_bdf_to_ugrid3d_filename(ugrid_filename, bdf_filename, ugrid_filename_out,
                                              pshell_pids_to_remove, tol=tol)


if __name__ == '__main__':   # pragma: no cover
    main()
