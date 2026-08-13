from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.bdf.cards.aero.zaero import ZAERO
    from pyNastran.bdf.cards.aero.zaero_cards.atm import ATMOS, FIXMATM, FIXHATM, FIXMACH, FIXMDEN
    from pyNastran.bdf.cards.aero.zaero_cards.geometry import (
        PAFOIL7,
        PAFOIL8,
    )
    from pyNastran.bdf.cards.aero.zaero_cards.flutter import MKAEROZ


def get_flutter_table(model: ZAERO, fix_id: int,
                      msg: str = "") -> FIXMATM | FIXHATM | FIXMACH | FIXMDEN:
    """gets a flutter table (FIXMATM, FIXHATM, FIXMACH, FIXMDEN)"""
    try:
        return model.flutter_table[fix_id]
    except KeyError:
        flutter_tables = np.unique(list(model.flutter_table.keys()))
        raise KeyError(f"fix_id={fix_id} not found{msg}.  Allowed flutter_tables={flutter_tables}")

def get_pafoil(model: ZAERO, pid: int, msg: str = "") -> PAFOIL7 | PAFOIL8:
    """gets a pafoil profile (PAFOIL7, PAFOIL8)"""
    try:
        return model.pafoil[pid]
    except KeyError:
        pafoils = np.unique(list(model.pafoil.keys()))
        raise KeyError(f"pid={pid} not found{msg}.  Allowed pafoils={pafoils}")


def get_mkaeroz(model: ZAERO, mkaeroz_id: int, msg: str = "") -> MKAEROZ:
    """gets an MKAEROZ"""
    try:
        return model.mkaeroz[mkaeroz_id]
    except KeyError:
        mkaerozs = np.unique(list(model.mkaeroz.keys()))
        raise KeyError(f"mkaeroz_id={mkaeroz_id} not found{msg}.  Allowed mkaerozs={mkaerozs}")

def get_atmos(model: ZAERO, atmos_id: int, msg: str = "") -> ATMOS:
    """gets an ATMOS"""
    try:
        return model.atmos[atmos_id]
    except KeyError:
        atmos_ids = np.unique(list(model.atmos.keys()))
        raise KeyError(f"atmos_id={atmos_id} not found{msg}.  Allowed atmos_ids={atmos_ids}")


def get_extfile(model: ZAERO, filename: str | int) -> str | None:
    if isinstance(filename, str):
        filename_ref = None
    elif isinstance(filename, int):
        filename_ref = model.extfile[filename]
    else:
        raise TypeError(f"filename={filename!r} type={str(filename)}")
    return filename_ref
