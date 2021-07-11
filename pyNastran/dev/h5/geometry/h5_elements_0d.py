from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
import h5py
from ..h5_utils import read_basic_element
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF

def read_conm2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'G', 'CID', 'M', 'X1', 'X2', 'X3', 'I1', 'I2', 'I3', 'DOMAIN_ID')
    EID = group['EID']
    NID = group['G']
    CID = group['CID']
    MASS = group['M']

    X1 = group['X1']
    X2 = group['X2']
    X3 = group['X3']

    I1 = group['I1']
    I2 = group['I2']
    I3 = group['I3']

    X = np.stack([X1, X2, X3], axis=1)
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, nid, cid, mass, x, i1, i2, i3 in zip(EID, NID, CID, MASS, X, I1, I2, I3):
        i11 = i1
        i21, i22 = i2
        i31, i32, i33 = i3
        i = [i11, i21, i22, i31, i32, i33]
        obj = geom_model.add_conm2(eid, nid, mass, cid=cid, X=x, I=i, comment='')
        obj.validate()

def read_celas1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    EID = group['EID']
    PID = group['PID']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, c1, c2 in zip(EID, PID, NIDS, C1, C2):
        obj = geom_model.add_celas1(eid, pid, nids, c1=c1, c2=c2, comment='')
        obj.validate()

def read_celas2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'K', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    EID = group['EID']
    K = group['K']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, k, nids, c1, c2 in zip(EID, K, NIDS, C1, C2):
        obj = geom_model.add_celas2(eid, k, nids, c1=c1, c2=c2, comment='')
        obj.validate()

def read_celas3(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    PID = group['PID']
    G1 = group['S1']
    G2 = group['S2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids in zip(EID, PID, NIDS):
        obj = geom_model.add_celas3(eid, pid, nids, comment='')
        obj.validate()

def read_celas4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    K = group['K']
    G1 = group['S1']
    G2 = group['S2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, k, nids in zip(EID, K, NIDS):
        obj = geom_model.add_celas4(eid, k, nids)
        obj.validate()

def read_cdamp1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    EID = group['EID']
    PID = group['PID']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids, c1, c2 in zip(EID, PID, NIDS, C1, C2):
        obj = geom_model.add_cdamp1(eid, pid, nids, c1=c1, c2=c2, comment='')
        obj.validate()

def read_cdamp2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'B', 'G1', 'G2', 'C1', 'C2', 'DOMAIN_ID')
    EID = group['EID']
    B = group['B']
    G1 = group['G1']
    G2 = group['G2']
    C1 = group['C1']
    C2 = group['C2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, b, nids, c1, c2 in zip(EID, B, NIDS, C1, C2):
        obj = geom_model.add_cdamp2(eid, b, nids, c1=c1, c2=c2, comment='')
        obj.validate()

def read_cdamp3(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    PID = group['PID']
    G1 = group['S1']
    G2 = group['S2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, pid, nids in zip(EID, PID, NIDS):
        obj = geom_model.add_cdamp3(eid, pid, nids, comment='')
        obj.validate()

def read_cdamp4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    EID = group['EID']
    B = group['B']
    G1 = group['S1']
    G2 = group['S2']
    NIDS = np.stack([G1, G2], axis=1)
    assert NIDS.shape[1] == 2, NIDS.shape
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, b, nids in zip(EID, B, NIDS):
        obj = geom_model.add_cdamp4(eid, b, nids)
        obj.validate()

def read_cvisc(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'PID', 'G', 'DOMAIN_ID')
    read_basic_element(group, geom_model, geom_model.add_cvisc)

def read_cbush(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    asdf
def read_cbush1d(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    asdf
def read_cbush2d(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    asdf
