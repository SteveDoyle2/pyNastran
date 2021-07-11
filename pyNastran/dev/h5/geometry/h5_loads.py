from __future__ import annotations
from typing import TYPE_CHECKING
#import numpy as np
import h5py
from ..h5_utils import get_tree, passer
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.loads.static_loads import PLOAD1

def read_dload(*args):
    pass
def read_load(*args):
    pass

def read_pload(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('SID', 'P', 'G', 'DOMAIN_ID')
    SID = group['SID']
    P = group['P']
    G = group['G']
    DOMAIN_ID = group['DOMAIN_ID']
    for sid, pressure, nodes in zip(SID, P, G):
        if nodes[-1] == 0:
            obj = geom_model.add_pload(sid, pressure, nodes[:-1], comment='')
        else:
            obj = geom_model.add_pload(sid, pressure, nodes, comment='')
        obj.validate()
        str(obj)

def read_pload1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('SID', 'EID', 'TYPE', 'SCALE', 'X1', 'P1', 'X2', 'P2', 'DOMAIN_ID')
    SID = group['SID']
    EID = group['EID']
    TYPE = group['TYPE']
    SCALE = group['SCALE']
    X1 = group['X1']
    P1 = group['P1']
    X2 = group['X2']
    P2 = group['P2']
    DOMAIN_ID = group['DOMAIN_ID']
    for sid, eid, load_type, scale, x1, x2, p1, p2 in zip(SID, EID, TYPE, SCALE, X1, X2, P1, P2):
        load_type_str = PLOAD1.valid_types[load_type - 1]
        scale_str = PLOAD1.valid_scales[scale - 1]

        obj = geom_model.add_pload1(sid, eid, load_type_str, scale_str, x1, p1, x2=x2, p2=p2, comment='')
        obj.validate()
        str(obj)

def read_pload2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('SID', 'P', 'EID', 'DOMAIN_ID')
    SID = group['SID']
    P = group['P']
    EID = group['EID']
    DOMAIN_ID = group['DOMAIN_ID']
    for sid, pressure, eid in zip(SID, P, EID):
        eids = [eid]
        obj = geom_model.add_pload2(sid, pressure, eids, comment='')
        obj.validate()
        str(obj)

def read_pload4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('SID', 'EID', 'P', 'G1', 'G34', 'CID', 'N', 'SORL', 'LDIR', 'DOMAIN_ID')
    SID = group['SID']
    P = group['P']
    EID = group['EID']
    G1 = group['G1']
    G34 = group['G34']
    CID = group['CID']
    N = group['N']
    SORL = group['SORL']
    LDIR = group['LDIR']
    DOMAIN_ID = group['DOMAIN_ID']
    for sid, pressures, eid, g1, g34, cid, nvector, sorl, ldir in zip(
        SID, P, EID, G1, G34, CID, N, SORL, LDIR):
        g1 = g1 if g1 != 0 else None
        g34 = g34 if g34 != 0 else None
        sorl_str = sorl.strip().decode('latin1')
        ldir_str = ldir.strip().decode('latin1')
        eids = [eid]
        obj = geom_model.add_pload4(
            sid, eids, pressures, g1=g1, g34=g34, cid=cid,
            nvector=nvector, surf_or_line=sorl_str, line_load_dir=ldir_str, comment='')
        obj.validate()
        str(obj)

def read_tload1(*args):
    pass
def read_tload2(*args):
    pass
def read_rload1(*args):
    pass
def read_rload2(*args):
    pass

def read_force(*args):
    pass
def read_force1(*args):
    pass
def read_force2(*args):
    pass
def read_moment(*args):
    pass
def read_moment1(*args):
    pass
def read_moment2(*args):
    pass

def read_sload(*args):
    pass

load_map = {
    'DLOAD': read_dload,
    'LOAD': read_load,

    'PLOAD': read_pload,
    'PLOAD1': read_pload1,
    'PLOAD2': read_pload2,
    'PLOAD4': read_pload4,

    'TLOAD1': read_tload1,
    'TLOAD2': read_tload2,
    'RLOAD1': read_rload1,
    'RLOAD2': read_rload2,

    'FORCE': read_force,
    'FORCE1': read_force1,
    'FORCE2': read_force2,

    'MOMENT': read_moment,
    'MOMENT1': read_moment1,
    'MOMENT2': read_moment2,

    'SLOAD': read_sload,
    'DAREA': passer,
    'PLOADX1': passer,

    # --------------------------------
    # thermal
    'CONV': passer,
}
