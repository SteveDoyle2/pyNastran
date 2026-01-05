from __future__ import annotations
from typing import Any, TYPE_CHECKING
from pyNastran.utils import PathLike
from pyNastran.bdf.cards.aero.zaero import ZAERO
from pyNastran.bdf.bdf import BDF, read_bdf


def zaero_to_nastran(zaero_filename: PathLike | ZAERO | BDF,
                     save: bool=True) -> tuple[dict, list, bool]:
    """Converts a ZAERO model to Nastran"""
    model = get_zaero_bdf_model(zaero_filename)
    if model.nastran_format not in ['zona', 'zaero']:
        caeros = {}
        caero2s = []
        make_paero1 = False
        return caeros, caero2s, make_paero1

    zaero = model.zaero
    caeros, caero2s, make_paero1 = _convert_caeros(zaero)
    splines = _convert_splines(zaero)
    aesurf, aelists = _convert_aesurf_aelist(zaero)

    trims = _convert_trim(zaero)
    aeros, aero = zaero.model.aeros.convert_to_zona(zaero.model)

    aelinks = _convert_trimlnk(zaero)

    if save:
        zaero.clear()
        zaero.model.splines = splines
        zaero.model.aesurf = aesurf
        zaero.model.aelists = aelists
        zaero.model.aelinks = aelinks
        zaero.model.trims = trims
        zaero.model.aeros = aeros
        zaero.model.aero = aero
    return caeros, caero2s, make_paero1


def get_zaero_bdf_model(zaero_filename: PathLike | ZAERO | BDF) -> BDF:
    if isinstance(zaero_filename, ZAERO):
        model = zaero_filename.model
        return model

    if isinstance(zaero_filename, PathLike):
        model = read_bdf(zaero_filename, mode='zaero')
        return model

    if isinstance(zaero_filename, BDF):
        model = zaero_filename
        return model
    # else:  # pragma: no cover
    raise TypeError(type(zaero_filename))
    # return model


def _convert_caeros(zaero: ZAERO) -> dict[int, Any]:
    """Converts ZONA CAERO7/BODY7 to CAERO1/CAERO2"""
    model = zaero.model
    caeros = {}
    caero2s = []
    paero1_id = 1
    make_paero1 = False
    for caero_id, caero in sorted(model.caeros.items()):
        if caero.type == 'CAERO7':
            caero_new = caero.convert_to_nastran(paero1_id)
            make_paero1 = True
        elif caero.type == 'BODY7':
            caero2s.append(caero)
            continue
        else:
            raise NotImplementedError(caero)
        caeros[caero_id] = caero_new

    # if make_paero1:
    #     paero = model.add_paero1(paero_id)
    zaero.add_caero2s(caero2s, add=False)
    return caeros, caero2s, make_paero1


def _convert_aesurf_aelist(zaero: ZAERO) -> tuple[dict[int, Any], list]:
    """
    Converts ZONA AESURFZ to AESURF/AELIST

    +---------+--------+-------+-------+-------+--------+--------+
    |    1    |   2    |   3   |   4   |   5   |   6    |    7   |
    +=========+========+=======+=======+=======+========+========+
    | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |
    +---------+--------+-------+-------+-------+--------+--------+
    | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |
    +---------+--------+-------+-------+-------+--------+--------+
    """
    model = zaero.model
    aelist_id = max(model.aelists) + 1 if model.aelists else 1
    aesurf_id = aelist_id
    aesurf = {}
    aelists = {}
    for unused_aesurf_name, aesurfi in sorted(model.aesurf.items()):
        aelist, aesurfi2 = aesurfi.convert_to_nastran(model, aesurf_id, aelist_id)
        aelists[aelist.sid] = aelist
        aesurf[aesurfi2.aesurf_id] = aesurfi2
        aesurf_id += 1
        aelist_id += 1
    return aesurf, aelists


def _convert_splines(zaero: ZAERO) -> dict[int, Any]:
    """Converts ZONA splines to splines"""
    splines = {}
    model = zaero.model
    for unused_spline_id, spline in model.splines.items():
        # print(spline)
        if spline.type == 'SPLINE1_ZAERO':
            splines_new = spline.convert_to_nastran(model)
        elif spline.type == 'SPLINE3_ZAERO':
            splines_new = spline.convert_to_nastran(model)
        else:
            raise NotImplementedError(spline)
        for spline_new in splines_new:
            splines[spline.eid] = spline_new
    return splines


def _convert_trim(zaero: ZAERO) -> dict[int, Any]:
    """Converts ZAERO TRIM to TRIM"""
    trims = {}
    model = zaero.model
    for trim_id, trim in sorted(model.trims.items()):
        trim_new = trim.convert_to_nastran(model)
        trims[trim_id] = trim_new
    return trims

def _convert_trimlnk(zaero: ZAERO) -> dict[int, Any]:
    """Converts ZAERO TRIMLNK to AELINK"""
    model = zaero.model
    assert isinstance(model.aelinks, dict), model.aelinks
    aelinks = {}
    for trim_id, trimlnk in sorted(zaero.trimlnk.items()):
        aelink = trimlnk.convert_to_nastran(model)
        aelinks[trim_id] = aelink
    return aelinks
