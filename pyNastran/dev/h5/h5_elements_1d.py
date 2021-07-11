from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import h5py
from .h5_utils import read_basic_element
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF

def read_ctube(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G', 'DOMAIN_ID')
    read_basic_element(group, geom_model, geom_model.add_ctube)
    #EID = group['EID']
    #PID = group['PID']
    #NIDS = group['G']
    #DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids in zip(EID, PID, NIDS):
        #obj = geom_model.add_ctube(eid, pid, nids, comment='')
        #obj.validate()

def read_cbend(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    pass

def read_crod(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    read_basic_element(group, geom_model, geom_model.add_crod)
    #EID = group['EID']
    #PID = group['PID']
    #NIDS = group['G']
    #DOMAIN_ID = group['DOMAIN_ID']
    #for eid, pid, nids in zip(EID, PID, NIDS):
        #obj = geom_model.add_crod(eid, pid, nids, comment='')
        #obj.validate()

def read_conrod(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'G1', 'G2', 'MID', 'A', 'J', 'C', 'NSM', 'DOMAIN_ID')
    assert len(group.dtype.names) == 9, group.dtype.names
    EID = group['EID']
    G1 = group['G1']
    G2 = group['G2']
    MID = group['MID']
    A = group['A']
    J = group['J']
    C = group['C']
    NSM = group['NSM']

    NIDS = np.stack([G1, G2], axis=1)
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, mid, nids, a, j, c, nsm in zip(EID, MID, NIDS, A, J, C, NSM):
        obj = geom_model.add_conrod(eid, mid, nids, A=a, j=j, c=c, nsm=nsm)
        obj.validate()
    geom_model.card_count[name] = len(EID)

def read_cbar(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #setattr(geom_model, name, group)
    #return
    assert len(group.dtype.names) == 18, group.dtype.names
    EID = group['EID']
    PID = group['PID']
    GA = group['GA']
    GB = group['GB']
    G0 = group['GO']

    PA = group['PA']
    PB = group['PB']

    FLAG = group['FLAG']
    X1 = group['X1']
    X2 = group['X2']
    X3 = group['X3']

    W1A = group['W1A']
    W2A = group['W2A']
    W3A = group['W3A']

    W1B = group['W1B']
    W2B = group['W2B']
    W3B = group['W3B']
    NIDS = np.stack([GA, GB], axis=1)
    X = np.stack([X1, X2, X3], axis=1)
    WA = np.stack([W1A, W2A, W3A], axis=1)
    WB = np.stack([W1B, W2B, W3B], axis=1)
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, flag, x, g0, pa, pb, wa, wb in zip(EID, PID, NIDS,
                                                           FLAG, X, G0, PA, PB, WA, WB):
        if g0 > 0:
            x = None
        else:
            g0 = None
        assert flag in [1, 2], 'flag=(1) GGG; expected'
        obj = geom_model.add_cbar(eid, pid, nids, x, g0, offt='GGG',
                                  pa=pa, pb=pb, wa=wa, wb=wb, comment='')
        obj.validate()
        str(obj)
    geom_model.card_count[name] = len(EID)

def read_cbeam(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'GA', 'GB', 'SA', 'SB', 'X', 'G0', 'F', 'PA', 'PB', 'WA', 'WB', 'DOMAIN_ID')
    assert len(group.dtype.names) == 14, group.dtype.names
    EID = group['EID']
    PID = group['PID']
    GA = group['GA']
    GB = group['GB']
    G0 = group['G0']

    SA = group['SA']
    SB = group['SB']

    PA = group['PA']
    PB = group['PB']

    FLAG = group['F']
    X = group['X']
    #X1 = group['X1']
    #X2 = group['X2']
    #X3 = group['X3']

    WA = group['WA']
    WB = group['WB']

    #W1A = group['W1A']
    #W2A = group['W2A']
    #W3A = group['W3A']

    #W1B = group['W1B']
    #W2B = group['W2B']
    #W3B = group['W3B']
    NIDS = np.stack([GA, GB], axis=1)
    #X = np.stack([X1, X2, X3], axis=1)
    #WA = np.stack([W1A, W2A, W3A], axis=1)
    #WB = np.stack([W1B, W2B, W3B], axis=1)
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, flag, x, g0, pa, pb, wa, wb, sa, sb in zip(
        EID, PID, NIDS, FLAG, X, G0, PA, PB, WA, WB, SA, SB):
        if g0 > 0:
            x = None
        else:
            g0 = None
        assert flag == 1, 'flag=(1) GGG; expected'
        obj = geom_model.add_cbeam(eid, pid, nids, x, g0,
                                   pa=pa, pb=pb, wa=wa, wb=wb, sa=sa, sb=sb)
        obj.validate()
        str(obj)
    geom_model.card_count[name] = len(EID)
