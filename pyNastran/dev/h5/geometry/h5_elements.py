from __future__ import annotations
from typing import List, Callable, Any # , TYPE_CHECKING
import numpy as np
import h5py
from ..h5_utils import get_tree, passer
from .h5_elements_0d import (
    read_cbush, read_cbush1d, read_cbush2d,
    read_cdamp1, read_cdamp2, read_cdamp3, read_cdamp4,
    read_celas1, read_celas2, read_celas3, read_celas4,
    read_cvisc,
    read_conm2,
)
from .h5_elements_1d import (read_crod, read_conrod, read_ctube,
                             read_cbar, read_cbeam, read_cbend)
from .h5_elements_2d import (read_cshear, read_ctria3, read_cquad4,
                             read_cquad8, read_ctria6,
                             read_cquadr, read_ctriar)
#if TYPE_CHECKING:  # pragma: no cover
from pyNastran.bdf.bdf import BDF

GeomCallable = Callable[[h5py._hl.dataset.Dataset, BDF],  # inputs
                        None]  # output


def read_rbe2(name: str, group: h5py._hl.group.Group, geom_model: BDF) -> None:
    gm_group = group.get('GM')
    rb_group = group.get('RB')

    assert gm_group is not None, gm_group
    assert rb_group is not None, rb_group
    assert isinstance(group, h5py._hl.group.Group)
    #names = group.dtype.names

    EID = rb_group['EID']
    GN = rb_group['GN']
    CM = rb_group['CM']
    GM_POS = rb_group['GM_POS']
    GM_LEN = rb_group['GM_LEN']
    ALPHA = rb_group['ALPHA']
    DOMAIN = rb_group['DOMAIN_ID']

    GM = gm_group['ID']
    for eid, gn, cm, gm_pos, gm_len, alpha, domain in zip(EID, GN, CM, GM_POS, GM_LEN, ALPHA, DOMAIN):
        Gmi = GM[gm_pos:gm_pos+gm_len]
        assert len(Gmi) == gm_len
        rbe2 = geom_model.add_rbe2(eid, gn, cm, Gmi, alpha=alpha, comment='')
        rbe2.validate()

def read_plotel(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    #('EID', 'G', 'DOMAIN_ID')
    EID = group['EID']
    NIDS = group['G']
    DOMAIN_ID = group['DOMAIN_ID']
    for eid, nodes in zip(EID, NIDS):
        obj = geom_model.add_plotel(eid, nodes, comment='')
        obj.validate()


SolidElementCallable = Callable[
    [int, int, list[int], str],
    Any,
]
def _read_solid(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF,
                nnodes: int, add_card: SolidElementCallable) -> None:
    assert len(group.dtype.names) == 4, group.dtype.names
    geom_model.card_count[name] = len(group['EID'])
    setattr(geom_model, name, group)
    return
    #geom_model.CTETRA = group

    EID = group['EID']
    PID = group['PID']
    NIDS = group['G']
    DOMAIN_ID = group['DOMAIN_ID']
    quadratic_nids = NIDS[:, nnodes:]
    is_large = quadratic_nids.max(axis=1) == 1

    for eid, pid, nids in zip(EID[is_large], PID[is_large], NIDS[~is_large, nnodes:]):
        obj = geom_model.add_ctetra(eid, pid, nids, comment='')
        obj.validate()
    for eid, pid, nids in zip(EID[is_large], PID[is_large], NIDS[is_large, :]):
        nids2 = [nid if nid != 0 else nid for nid in nids]
        obj = geom_model.add_ctetra(eid, pid, nids2, comment='')
        obj.validate()

def read_ctetra(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    nnodes = 4
    _read_solid(name, group, geom_model,
                nnodes, geom_model.add_ctetra)

def read_cpenta(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    nnodes = 6
    _read_solid(name, group, geom_model,
                nnodes, geom_model.add_cpenta)

def read_chexa(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    nnodes = 8
    _read_solid(name, group, geom_model,
                nnodes, geom_model.add_chexa)

#------------------------------------------------------------
def read_caero1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    if geom_model.flags['aero'] is False:
        return
    #('EID', 'PID', 'CP', 'NSPAN', 'NCHORD', 'LSPAN', 'LCHORD', 'IGID', 'X1', 'Y1', 'Z1', 'X12', 'X4', 'Y4', 'Z4', 'X43', 'DOMAIN_ID')
    EID = group['EID']
    PID = group['PID']
    CP = group['CP']
    NSPAN = group['NSPAN']
    NCHORD = group['NCHORD']
    LSPAN = group['LSPAN']
    LCHORD = group['LCHORD']
    IGID = group['IGID']
    X1 = group['X1']
    Y1 = group['Y1']
    Z1 = group['Z1']
    X12 = group['X12']
    X4 = group['X4']
    Y4 = group['Y4']
    Z4 = group['Z4']

    X43 = group['X43']
    DOMAIN_ID = group['DOMAIN_ID']
    XYZ1 = np.stack([X1, Y1, Z1], axis=1)
    XYZ4 = np.stack([X4, Y4, Z4], axis=1)
    assert XYZ1.shape[1] == 3
    for eid, pid, cp, nspan, nchord, lspan, lchord, igroup, xyz1, x12, xyz4, x43 in zip(
        EID, PID, CP, NSPAN, NCHORD, LSPAN, LCHORD, IGID, XYZ1, X12, XYZ4, X43):
        obj = geom_model.add_caero1(
            eid, pid, igroup, xyz1, x12, xyz4, x43,
            cp=cp, nspan=nspan, lspan=lspan,
            nchord=nchord, lchord=lchord, comment='')
        obj.validate()
        str(obj)

element_map = {
    # mass
    'CMASS1': passer,
    'CMASS2': passer,
    'CMASS3': passer,
    'CMASS4': passer,
    'CONM2': read_conm2,

    # 0D
    'CELAS1': read_celas1,
    'CELAS2': read_celas2,
    'CELAS3': read_celas3,
    'CELAS4': read_celas4,

    'CDAMP1': read_cdamp1,
    'CDAMP2': read_cdamp2,
    'CDAMP3': read_cdamp3,
    'CDAMP4': read_cdamp4,
    'CVISC': read_cvisc,
    'CBUSH': read_cbush,
    'CBUSH1D': read_cbush1d,
    'CBUSH2D': read_cbush2d,

    # 1D
    'CROD': read_crod,
    'CONROD': read_conrod,
    'CBAR': read_cbar,
    'CBEAM': read_cbeam,
    'CTUBE': read_ctube,
    'CBEND': read_cbend,

    # 2D
    'CTRIA3': read_ctria3,
    'CTRIA6': read_ctria6,
    'CQUAD4': read_cquad4,
    'CQUAD8': read_cquad8,
    'CTRIAR': read_ctriar,
    'CQUADR': read_cquadr,
    'CSHEAR': read_cshear,
    'PLOTEL': read_plotel,

    # 3d
    'CTETRA': read_ctetra,
    'CPENTA': read_cpenta,
    'CHEXA': read_chexa,

    # axisymmetric
    'CTRIAX6': passer,

    # rigid
    'RBAR': passer,
    'RROD': passer,
    'RTRPLT': passer,

    'RBE1': passer,
    'RBE2': read_rbe2,
    'RBE3': passer,

    # other
    'GENEL': passer,
    # --------------------------
    # thermal
    'CHBDYE': passer,
    # --------------------------
    #non-elements
    'CAERO1' : read_caero1,
}
