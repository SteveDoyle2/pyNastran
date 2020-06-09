"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import List, Union, TYPE_CHECKING

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header
if TYPE_CHECKING:
    from pyNastran.bdf.cards.aero.static_loads import AEROS # , AESTAT, CSSCHD, DIVERG, TRIM, TRIM2
    from pyNastran.bdf.cards.aero.dynamic_loads import AERO # , FLUTTER, FLFACT, MKAERO1, MKAERO2
    from pyNastran.op2.op2_geom import OP2Geom, BDF

def write_edt(op2, op2_ascii, model, endian=b'<'):
    """writes the EDT/EDTS table"""
    card_types = [
        'MKAERO1',
        'AERO', 'AEROS',
        'CAERO1', 'PAERO1',
    ]

    cards_to_skip = [
    ]
    #card_types = [
        #'CAERO1', 'CAERO2',
        #'SPLINE1', 'SPLINE2',
        #'AEFACT', 'AELIST', 'AESTAT',
    #]
    out = defaultdict(list)
    # geometry
    for eid, caero in sorted(model.caeros.items()):
        out[caero.type].append(eid)
    for pid, paero in sorted(model.paeros.items()):
        out[paero.type].append(pid)
    for spline_id, spline in sorted(model.splines.items()):
        out[spline.type].append(spline_id)
    for aelink_id, aelink in sorted(model.aelinks.items()):
        out[aelink.type].append(aelink_id)
    for name, aecomp in sorted(model.aecomps.items()):
        out[aecomp.type].append(name)
    #for name, aecompl in sorted(model.aecompl.items()):
        #out[aecompl.type].append(name)

    # loads
    for trim_id, trim in sorted(model.trims.items()):
        out[trim.type].append(trim_id)
    for aestat_id, aestat in sorted(model.aestats.items()):
        out[aestat.type].append(aestat_id)
    if model.aeros:
        out[model.aeros.type].append(model.aeros)

    # flutter
    for flutter_id, flutter in sorted(model.flutters.items()):
        out[flutter.type].append(flutter_id)
    for flfact_id, flfact in sorted(model.flfacts.items()):
        out[flfact.type].append(flfact_id)
    for mkaero in model.mkaeros:
        out[mkaero.type].append(mkaero)
    if model.aero:
        out[model.aero.type].append(model.aero)

    for card_type in list(out):
        if card_type not in card_types:
            del out[card_type]
            model.log.warning(f'removing {card_type} in OP2 writer')
    # other
    if len(out) == 0:
        return

    write_geom_header(b'EDT', op2, op2_ascii, endian=endian)
    itable = -3

    for name, ids in sorted(out.items()):
        model.log.debug('EDT %s %s' % (name, ids))
        ncards = len(ids)
        assert ncards > 0, ncards
        if name in cards_to_skip:
            model.log.warning('skipping EDT-%s' % name)
            continue

        #if nmaterials == 0:
            #continue
        if name == 'CAERO1':
            nbytes = _write_caero1(model, name, ids, ncards, op2, op2_ascii, endian)
        elif name == 'PAERO1':
            nbytes = _write_paero1(model, name, ids, ncards, op2, op2_ascii, endian)
        elif name == 'MAKERO1':
            nbytes = _write_makero1(model, name, ids, ncards, op2, op2_ascii, endian)
        elif name == 'AERO':
            nbytes = _write_aero(model, name, ids, ncards, op2, op2_ascii, endian)
        elif name == 'AEROS':
            nbytes = _write_aeros(model, name, ids, ncards, op2, op2_ascii, endian)
        else:  # pragma: no cover
            raise NotImplementedError(name)
        op2.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2, op2_ascii, itable)
    #-------------------------------------

def _write_caero1(model: Union[BDF, OP2Geom], name: str,
                  caero_ids: List[int], ncards: int,
                  op2, op2_ascii, endian: bytes) -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 EID    I
    2 PID    I
    3 CP     I
    4 NSPAN  I
    5 NCHORD I
    6 LSPAN  I
    7 LCHORD I
    8 IGID   I
    9 X1    RS
    10 Y1   RS
    11 Z1   RS
    12 X12  RS
    13 X4   RS
    14 Y4   RS
    15 Z4   RS
    16 X43  RS

    """
    key = (3002, 30, 263)
    nfields = 16
    structi = Struct(endian + b'8i 8f')
    nbytes = write_header(name, nfields, ncards, key, op2, op2_ascii)

    for caero_id in caero_ids:
        caero = model.caeros[caero_id]
        x1, y1, z1 = caero.p1
        x4, y4, z4 = caero.p4
        #print(caero.get_stats())
        data = [caero.eid, caero.pid, caero.cp,
                caero.nspan, caero.nchord,
                caero.lspan, caero.lchord,
                caero.igroup, x1, y1, z1, caero.x12,
                x4, y4, z4, caero.x43]

        assert None not in data, data
        op2_ascii.write(f'  CAERO1 data={data}\n')
        op2.write(structi.pack(*data))
    return nbytes

def _write_paero1(model: Union[BDF, OP2Geom], name: str,
                  paero_ids: List[int], ncards: int,
                  op2, op2_ascii, endian: bytes) -> int:
    """
    (3102, 31, 264)
    MSC 2018.2

    Word Name Type Description
    1 PID      I
    2 B1       I
    3 B2       I
    4 B3       I
    5 B4       I
    6 B5       I
    7 B6       I
    8 UNDEF none

    PAERO1  100001

    """
    key = (3102, 31, 264)
    nfields = 8
    structi = Struct(endian + b'8i')
    nbytes = write_header(name, nfields, ncards, key, op2, op2_ascii)

    for paero_id in paero_ids:
        paero = model.paeros[paero_id]
        n0s = 7 - len(paero.caero_body_ids)
        data = [paero.pid, ] + paero.caero_body_ids
        if n0s:
            data += [0] * n0s
        assert None not in data, data
        op2_ascii.write(f'  PAERO1 data={data}\n')
        op2.write(structi.pack(*data))
    return nbytes

def _write_makero1(model: Union[BDF, OP2Geom], name: str,
                   mkaero1s, ncards: int,
                   op2, op2_ascii, endian: bytes) -> int:
    """writes the MKAERO1

    data = (1.3, -1, -1, -1, -1, -1, -1, -1,
           0.03, 0.04, 0.05, -1, -1, -1, -1, -1)
    """
    key = (3802, 38, 271)
    nfields = 16 * ncards
    #spack = Struct(endian + b'i10fi')
    nbytes = write_header(name, nfields, ncards, key, op2, op2_ascii)

    for mkaero in mkaero1s:
        data = []
        nmachs = len(mkaero.machs)
        nkfreqs = len(mkaero.reduced_freqs)
        nint_mach = 8 - nmachs
        nint_kfreq = 8 - nkfreqs
        fmt1 = b'%if' % nmachs + b'i' * nint_mach
        fmt2 = b'%if' % nkfreqs + b'i' * nint_kfreq
        spack = Struct(endian + fmt1 + fmt2)
        data.extend(mkaero.machs.tolist())
        if nint_mach:
            data.extend([-1]*nint_mach)
        data.extend(mkaero.reduced_freqs.tolist())
        if nint_mach:
            data.extend([-1]*nint_mach)

        op2_ascii.write(f'  mkaero1 data={data}\n')
        op2.write(spack.pack(*data))
    return nbytes

def _write_aero(model: Union[BDF, OP2Geom], name: str,
                aero: List[AERO], ncards: int,
                op2, op2_ascii, endian: bytes) -> int:
    """
    Word Name Type Description
    1 ACSID     I
    2 VELOCITY RS
    3 REFC     RS
    4 RHOREF   RS
    5 SYMXZ     I
    6 SYMXY     I

    """
    aeroi = aero[0]
    #print(aeroi.get_stats())
    spack = Struct(endian + b'i 3f 2i')

    key = (3202, 32, 265)
    nfields = 6
    nbytes = write_header(name, nfields, ncards, key, op2, op2_ascii)

    # cref   : 0.3048
    # is_anti_symmetric_xy : False
    # is_anti_symmetric_xz : False
    # is_symmetric_xy : False
    # is_symmetric_xz : False
    # rho_ref : 1.225
    # sym_xy : 0
    # sym_xz : 0
    # type   : 'AERO'
    # velocity : 0.0
    velocity = 0.0
    if aeroi.velocity is not None:
        velocity = aeroi.velocity
    data = [aeroi.acsid, velocity, aeroi.cref, aeroi.rho_ref, aeroi.sym_xz, aeroi.sym_xy]
    op2_ascii.write(f'  AERO data={data}\n')
    op2.write(spack.pack(*data))
    return nbytes

def _write_aeros(model: Union[BDF, OP2Geom], name: str,
                 aeros: List[AEROS], ncards: int,
                 op2, op2_ascii, endian: bytes) -> int:
    """
    AEROS(2202, 22, 340)

    AEROS   0       100     36.     360.    12960.
    data = (0, 100, 36.0, 360.0, 12960.0, 0, 0)

    """
    aeroi = aeros[0]
    #print(aeroi.get_stats())
    spack = Struct(endian + b'2i 3f 2i')

    key = (2202, 22, 340)
    nfields = 7
    nbytes = write_header(name, nfields, ncards, key, op2, op2_ascii)

    data = [aeroi.acsid, aeroi.rcsid,
            aeroi.cref, aeroi.bref, aeroi.sref,
            aeroi.sym_xz, aeroi.sym_xy]
    print(data)
    assert None not in data, data
    op2_ascii.write(f'  AEROS data={data}\n')
    op2.write(spack.pack(*data))
    return nbytes
