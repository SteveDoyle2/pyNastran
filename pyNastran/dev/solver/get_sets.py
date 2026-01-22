import numpy as np
from pyNastran.bdf.bdf import BDF


def get_aset(model: BDF) -> set[tuple[int, int]]:
    aset_map = set()
    for aset in model.asets:
        if aset.type == 'ASET1':
            comp = aset.components
            for nid in aset.ids:
                for compi in comp:
                    aset_map.add((nid, int(compi)))
        elif aset.type == 'ASET':
            for nid, comp in zip(aset.ids, aset.components):
                for compi in comp:
                    aset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(aset)
    return aset_map


def get_bset(model: BDF) -> set[tuple[int, int]]:
    """creates the b-set"""
    bset_map = set()
    for bset in model.bsets:
        if bset.type == 'BSET1':
            comp = bset.components
            for nid in bset.ids:
                for compi in comp:
                    bset_map.add((nid, int(compi)))
        elif bset.type == 'BSET':
            for nid, comp in zip(bset.ids, bset.components):
                for compi in comp:
                    bset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(bset)
    return bset_map


def get_cset(model: BDF) -> set[tuple[int, int]]:
    """creates the c-set"""
    cset_map = set()
    for cset in model.csets:
        if cset.type == 'CSET1':
            comp = cset.components
            for nid in cset.ids:
                for compi in comp:
                    cset_map.add((nid, int(compi)))
        elif cset.type == 'CSET':
            for nid, comp in zip(cset.ids, cset.components):
                for compi in comp:
                    cset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(cset)
    return cset_map


def get_omit_set(model: BDF) -> set[tuple[int, int]]:
    """creates the o-set"""
    omit_set_map = set()
    for omit in model.omits:
        if omit.type == 'OMIT1':
            comp = omit.components
            for nid in omit.ids:
                for compi in comp:
                    omit_set_map.add((nid, int(compi)))
        elif omit.type == 'OMIT':
            for nid, comp in zip(omit.ids, omit.components):
                for compi in comp:
                    omit_set_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(omit)
    return omit_set_map


def get_rset(model: BDF) -> set[tuple[int, int]]:
    """creates the r-set"""
    rset_map = set()
    for rset in model.suport:
        for nid, comp in zip(rset.ids, rset.components):
            for compi in comp:
                rset_map.add((nid, int(compi)))

    for suport in model.suport1:
        comp = suport.components
        for nid in suport.ids:
            for compi in comp:
                rset_map.add((nid, int(compi)))
    return rset_map


def get_qset(model: BDF) -> set[tuple[int, int]]:
    """creates the q-set"""
    qset_map = set()
    for qset in model.qsets:
        if qset.type == 'QSET1':
            comp = qset.components
            for nid in qset.ids:
                for compi in comp:
                    qset_map.add((nid, int(compi)))
        elif qset.type == 'QSET':
            for nid, comp in zip(qset.ids, qset.components):
                for compi in comp:
                    qset_map.add((nid, int(compi)))
        else:
            raise NotImplementedError(qset)
    return qset_map
