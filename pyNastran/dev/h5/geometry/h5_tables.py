from __future__ import annotations
from typing import Callable, TYPE_CHECKING
#import numpy as np
import h5py
from ..h5_utils import get_tree, passer
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def read_tabled1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table1(name, group, geom_model, geom_model.add_tabled1)

def read_tabled2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table2(name, group, geom_model, geom_model.add_tabled2)

def read_tabled3(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table3(name, group, geom_model, geom_model.add_tabled3)

def read_tabled4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table4(name, group, geom_model, geom_model.add_tabled4)

def read_tablem1(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table1(name, group, geom_model, geom_model.add_tablem1)

def read_tablem2(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table2(name, group, geom_model, geom_model.add_tablem2)

def read_tablem3(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table3(name, group, geom_model, geom_model.add_tablem3)

def read_tablem4(name: str, group: h5py._hl.dataset.Dataset, geom_model: BDF) -> None:
    _read_table4(name, group, geom_model, geom_model.add_tablem4)

def _read_table1(name: str, group: h5py._hl.group.Group, geom_model: BDF,
                 add_table1: Callable) -> None:
    identity = group.get('IDENTITY')
    TID = identity['ID']
    CODEX = identity['CODEX']
    CODEY = identity['CODEY']
    POS = identity['POS']
    LEN = identity['LEN']
    DOMAIN_ID = identity['DOMAIN_ID']

    axis_map = {
        0: 'LINEAR',
    }
    xy = group.get('XY')
    X = xy['X']
    Y = xy['Y']
    for tid, codex, codey, pos, leni in zip(TID, CODEX, CODEY, POS, LEN):
        if tid >= 100_000_000:
            continue
        x = X[pos:pos+leni]
        y = Y[pos:pos+leni]
        obj = add_table1(
            tid, x, y,
            xaxis=axis_map[codex], yaxis=axis_map[codey],
            extrap=0, comment='')
        obj.validate()
        str(obj)

def _read_table2(name: str,
                 group: h5py._hl.group.Group,
                 geom_model: BDF,
                 add_table2: Callable):
    #identity = ('ID', 'X1', 'POS', 'LEN', 'DOMAIN_ID')
    identity = group.get('IDENTITY')
    TID = identity['ID']
    X1 = identity['X1']
    POS = identity['POS']
    LEN = identity['LEN']
    DOMAIN_ID = identity['DOMAIN_ID']

    xy = group.get('XY')
    X = xy['X']
    Y = xy['Y']
    for tid, x1, pos, leni in zip(TID, X1, POS, LEN):
        if tid >= 100_000_000:
            continue
        x = X[pos:pos+leni]
        y = Y[pos:pos+leni]
        obj = add_table2(tid, x1, x, y, comment='')
        obj.validate()
        str(obj)

def _read_table3(name: str,
                 group: h5py._hl.group.Group,
                 geom_model: BDF,
                 add_table3: Callable):
    #('ID', 'X1', 'X2', 'POS', 'LEN', 'DOMAIN_ID')
    identity = group.get('IDENTITY')
    xy = group.get('XY')

    TID = identity['ID']
    X1 = identity['X1']
    X2 = identity['X2']
    POS = identity['POS']
    LEN = identity['LEN']
    DOMAIN_ID = identity['DOMAIN_ID']

    xy = group.get('XY')
    X = xy['X']
    Y = xy['Y']
    for tid, x1, x2, pos, leni in zip(TID, X1, X2, POS, LEN):
        if tid >= 100_000_000:
            continue
        x = X[pos:pos+leni]
        y = Y[pos:pos+leni]
        obj = add_table3(tid, x1, x2, x, y)
        obj.validate()
        str(obj)

def _read_table4(name: str,
                 group: h5py._hl.group.Group,
                 geom_model: BDF,
                 add_table4: Callable):
    #('ID', 'X1', 'X2', 'X3', 'X4', 'POS', 'LEN', 'DOMAIN_ID')
    identity = group.get('IDENTITY')
    coef = group.get('COEF')

    TID = identity['ID']
    X1 = identity['X1']
    X2 = identity['X2']
    X3 = identity['X3']
    X4 = identity['X4']
    POS = identity['POS']
    LEN = identity['LEN']
    DOMAIN_ID = identity['DOMAIN_ID']

    A = coef['A']
    for tid, x1, x2, x3, x4, pos, leni in zip(TID, X1, X2, X3, X4, POS, LEN):
        if tid >= 100_000_000:
            continue
        a = A[pos:pos+leni]
        obj = add_table4(tid, x1, x2, x3, x4, a)
        obj.validate()
        str(obj)

table_map = {
    'MKAERO1': passer,
    'TABLED1': read_tabled1,
    'TABLED2': read_tabled2,
    'TABLED3': read_tabled3,
    'TABLED4': read_tabled4,

    'TABLEM1': read_tablem1,
    'TABLEM2': read_tablem2,
    'TABLEM3': read_tablem3,
    'TABLEM4': read_tablem4,
}
