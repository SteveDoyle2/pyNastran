"""writes the MPT/MPTS table"""
from __future__ import annotations
from copy import deepcopy
from collections import defaultdict
from struct import pack, Struct
from typing import Union, Any, TYPE_CHECKING

import numpy as np

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues

integer_types = (int, np.int32, np.int64)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.cards.aero.aero import CAERO1, CAERO2, PAERO1, PAERO2, AESURF, AESURFS
    from pyNastran.bdf.cards.aero.static_loads import AEROS, AESTAT # , CSSCHD, DIVERG, TRIM, TRIM2
    from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, FLUTTER # , FLFACT, MKAERO2
    from pyNastran.op2.op2_geom import OP2Geom, BDF

def write_edt(op2_file, op2_ascii, model: Union[BDF, OP2Geom],
              endian: bytes=b'<', nastran_format: str='nx') -> None:
    """writes the EDT/EDTS table"""
    if not hasattr(model, 'loads'):  # OP2
        return
    card_types = [
        'MKAERO1', # 'MKAERO2',
        'AERO', 'AEROS',
        'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        'PAERO1', 'PAERO2', 'PAERO5',
        'SPLINE1', 'SPLINE2', #'SPLINE3',
        'AELIST',
        'AEFACT',
        'AESURF', 'AESURFS',
        'AESTAT',
        'TRIM', 'DIVERG', 'FLUTTER',
        'DEFORM',
        'FLFACT',
        'SET1',
        'AELINK',
        'MONPNT1', 'MONPNT2', 'MONPNT3',
    ]

    cards_to_skip = [
        #'GUST',  # part of DIT
    ]
    out = defaultdict(list)
    # geometry
    for unused_load_id, loads in model.loads.items():
        for load in loads:
            if load.type in ['DEFORM', 'CLOAD']:
                out[load.type].append(load)

    for set_id, seti in sorted(model.sets.items()):
        out[seti.type].append(set_id)
    for eid, caero in sorted(model.caeros.items()):
        out[caero.type].append(eid)
    for pid, paero in sorted(model.paeros.items()):
        out[paero.type].append(pid)
    for spline_id, spline in sorted(model.splines.items()):
        out[spline.type].append(spline_id)
    for aelist_id, aelist in sorted(model.aelists.items()):
        out[aelist.type].append(aelist_id)
    for aesurf_id, aesurf in sorted(model.aesurf.items()):
        out[aesurf.type].append(aesurf_id)
    for aesurfs_id, aesurfs in sorted(model.aesurfs.items()):
        out[aesurfs.type].append(aesurfs_id)
    for aelink_id, aelinks in sorted(model.aelinks.items()):
        for aelink in aelinks:
            out[aelink.type].append(aelink_id)
    for name, aecomp in sorted(model.aecomps.items()):
        out[aecomp.type].append(name)
    #for name, aecompl in sorted(model.aecompl.items()):
        #out[aecompl.type].append(name)

    # loads
    for diverg_id, diverg in sorted(model.divergs.items()):
        out[diverg.type].append(diverg_id)
    for trim_id, trim in sorted(model.trims.items()):
        out[trim.type].append(trim_id)
    #for gust_id, gust in sorted(model.gusts.items()):  # part of DIT
        #out[gust.type].append(gust_id)
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

    for monitor_id, monitor in enumerate(model.monitor_points):
        out[monitor.type].append(monitor_id)
    if model.aero:
        out[model.aero.type].append(model.aero)

    remove_unsupported_cards(out, card_types, model.log)
    # other
    if len(out) == 0:
        return

    write_geom_header(b'EDT', op2_file, op2_ascii, endian=endian)
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
        #model.log.info(f'EDT: {name}; n={len(ids)}')
        try:
            func = EDT_MAP[name]
        except KeyError:  # pragma: no cover
            #continue
            raise NotImplementedError(name)

        nbytes = func(model, name, ids, ncards, op2_file, op2_ascii,
                      endian, nastran_format)
        op2_file.write(pack('i', nbytes))
        itable -= 1
        data = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4]
        op2_file.write(pack('9i', *data))
        op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    #print('itable', itable)
    close_geom_table(op2_file, op2_ascii, itable)
    #-------------------------------------

def remove_unsupported_cards(card_dict: dict[str, Any],
                             card_types: list[str],
                             log: SimpleLogger):

    for card_type in list(card_dict):
        if card_type not in card_types:
            del card_dict[card_type]
            log.warning(f"removing {card_type} in OP2 writer because it's unsupported")

def write_trim(model: Union[BDF, OP2Geom], name: str,
               trim_ids: list[int], ncards: int,
               op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    (2402, 24, 342)
    MSC 2018.2

    Word Name Type Description
    1 ID           I
    2 MACH        RS
    3 Q           RS
    4 AEQRATIO    RS

    5 LABEL(2) CHAR4
    7 UX          RS
    Words 5 through 7 repeat until (-1,-1,-1) occurs

    """
    #struct1 = Struct(endian + b'i 3f')
    #struct2 = Struct(endian + b'8sf')
    #struct_end = Struct(endian + b'3i')
    key = (2402, 24, 342)
    #nfields = 16

    nlabels_total = 0
    for trim_id in trim_ids:
        trim = model.trims[trim_id]
        assert len(trim.labels) == len(trim.uxs), trim.get_stats()
        nlabels = len(trim.labels)
        #nuxs = len(trim.uxs)
        nlabels_total += nlabels # + nuxs

    # the *3 comes from:
    #  - label (2 fields of 4s)
    #  - ux
    nvalues = 7 * ncards + nlabels_total * 3
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for trim_id in trim_ids:
        trim = model.trims[trim_id]
        nlabels = len(trim.labels)
        data = [trim.sid, trim.mach, trim.q, trim.aeqr]
        for label, ux in zip(trim.labels, trim.uxs):
            label_bytes = b'%-8s' % label.encode('latin1')
            data.extend([label_bytes, ux])
        data += [-1, -1, -1]
        assert None not in data, data
        fmt = endian + b'ifff ' + b' 8sf' * nlabels + b' 3i'
        structi = Struct(fmt)
        op2_ascii.write(f'  TRIM data={data}\n')
        op2_file.write(structi.pack(*data))
        all_data += data
    return nbytes

def write_diverg(model: Union[BDF, OP2Geom], name: str,
                 diverg_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Record – DIVERG(2702,27,387)
    Divergence analysis data.
    Word Name Type Description
    1 SID     I    Unique set identification number
    2 NROOT   I    Number of divergence roots to output
    3 M       RS   Mach number
    Word 3 repeats until -1 occurs
    """
    key = (2702, 27, 387)

    nvalues = 0
    all_data = []
    fmt = b''
    for diverg_id in diverg_ids:
        diverg = model.divergs[diverg_id]
        #assert 1 == 0, trim.get_stats()

        # machs  : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        # nroots : 21
        # sid    : 100
        nmachs = len(diverg.machs)
        nvalues += 3 + nmachs


    #nvalues = 7 * ncards + nlabels_total * 3
    assert nvalues > 0, nvalues
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for diverg_id in diverg_ids:
        diverg = model.divergs[diverg_id]
        data = [diverg.sid, diverg.nroots] + list(diverg.machs) + [-1]
        fmt = endian + f'2i {nmachs}f i'.encode('ascii')

        structi = Struct(fmt)
        #print(f'  DIVERG data={data}\n')
        op2_ascii.write(f'  DIVERG data={data}\n')
        op2_file.write(structi.pack(*data))
        all_data += data
    assert len(data) == nvalues, f'ndata={len(data)}; nvalues={nvalues}'
    return nbytes

def write_caero1(model: Union[BDF, OP2Geom], name: str,
                 caero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
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
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for caero_id in caero_ids:
        caero = model.caeros[caero_id] # type: CAERO1
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
        op2_file.write(structi.pack(*data))
    return nbytes

def write_caero2(model: Union[BDF, OP2Geom], name: str,
                 caero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 EID          I
    2 PID          I
    3 CP           I
    4 NSB          I
    5 NINT         I
    6 LSB          I
    7 LINT         I
    8 IGID         I
    9 X1          RS
    10 Y1         RS
    11 Z1         RS
    12 X12        RS
    13 UNDEF(4) none

    data = (54000, 4020, 0, 8, 8, 0, 0, 1, -5.0, 0, 0, 40.0, 0, 0, 0, 0),
    """
    key = (4301, 43, 167)
    nfields = 16
    structi = Struct(endian + b'8i 8f')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for caero_id in caero_ids:
        caero = model.caeros[caero_id] # type: CAERO2
        x1, y1, z1 = caero.p1
        #print(caero.get_stats())
        data = [caero.eid, caero.pid, caero.cp,
                caero.nsb, caero.nint, caero.lsb, caero.lint,
                caero.igroup, x1, y1, z1, caero.x12,
                0.0, 0.0, 0.0, 0.0]

        assert None not in data, data
        op2_ascii.write(f'  CAERO2 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_caero3(model: Union[BDF, OP2Geom], name: str,
                 caero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Aerodynamic panel element configuration.

    Word Name Type Description

    1 EID    I Element identification number
    2 PID    I Property identification number of a PAERO3 entry
    3 CP     I Coordinate system for locating points 1 and 4
    4 LISTW  I Identification number of an AEFACT entry that lists
               coordinate pairs for structural interpolation of the wing
    5 LISTC1 I Identification number of an AEFACT entry that lists
               coordinate pairs for control surfaces
    6 LISTC2 I Identification number of an AEFACT entry that lists
               coordinate pairs for control surfaces
    7 UNDEF(2) None
    9  X1   RS X-coordinate of point 1 in coordinate system CP
    10 Y1   RS Y-coordinate of point 1 in coordinate system CP
    11 Z1   RS Z-coordinate of point 1 in coordinate system CP
    12 X12  RS Edge chord length in aerodynamic coordinate system
    13 X4   RS X-coordinate of point 4 in coordinate system CP
    14 Y4   RS Y-coordinate of point 4 in coordinate system CP
    15 Z4   RS Z-coordinate of point 4 in coordinate system CP
    16 X43  RS Edge chord length in aerodynamic coordinate system
    """
    #ntotal = 64 * self.factor # 4*16
    #ndatai = len(data) - n
    #ncards = ndatai // ntotal
    #assert ndatai % ntotal == 0
    structi = Struct(endian + b'6i 2i 8f')

    key = (4401, 44, 168)
    nfields = 16
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for caero_id in caero_ids:
        caero: CAERO3 = model.caeros[caero_id]
        x1, y1, z1 = caero.p1
        x4, y4, z4 = caero.p4
        #print(caero.get_stats())
        list_w = 0 if caero.list_w is None else caero.list_w
        list_c1 = 0 if caero.list_c1 is None else caero.list_c1
        list_c2 = 0 if caero.list_c2 is None else caero.list_c2
        data = [
            caero.eid, caero.pid, caero.cp,
            list_w, list_c1, list_c2,
            0, 0,
            x1, y1, z1, caero.x12,
            x4, y4, z4, caero.x43,]

        assert None not in data, data
        op2_ascii.write(f'  CAERO3 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_caero4(model: Union[BDF, OP2Geom], name: str,
                 caero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Word Name Type Description
    1 EID   I Element identification number
    2 PID   I Property identification number of a PAERO4 entry
    3 CP    I Coordinate system for locating points 1 and 4
    4 NSPAN I Number of strips
    5 LSPAN I Identification number of an AEFACT entry
              containing a list of division points for strips
    6 UNDEF(3) None
    9  X1  RS X-coordinate of point 1 in coordinate system CP
    10 Y1  RS Y-coordinate of point 1 in coordinate system CP
    11 Z1  RS Z-coordinate of point 1 in coordinate system CP
    12 X12 RS Edge chord length in aerodynamic coordinate system
    13 X4  RS X-coordinate of point 4 in coordinate system CP
    14 Y4  RS Y-coordinate of point 4 in coordinate system CP
    15 Z4  RS Z-coordinate of point 4 in coordinate system CP
    16 X43 RS Edge chord length in aerodynamic coordinate system
    """
    structi = Struct(endian + b'5i 3i 8f')

    key = (4501, 45, 169)
    nfields = 16
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for caero_id in caero_ids:
        caero = model.caeros[caero_id] # type: CAERO3
        x1, y1, z1 = caero.p1
        x4, y4, z4 = caero.p4
        #print(caero.get_stats())
        # 4 NSPAN I Number of strips
        # 5 LSPAN I Identification number of an AEFACT entry
        #           containing a list of division points for strips
        # 6 UNDEF(3) None

        data = [
            caero.eid, caero.pid, caero.cp,
            caero.nspan, caero.lspan, 0,
            0, 0,
            x1, y1, z1, caero.x12,
            x4, y4, z4, caero.x43,]

        assert None not in data, data
        op2_ascii.write(f'  CAERO4 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_caero5(model: Union[BDF, OP2Geom], name: str,
                 caero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 EID        I
    2 PID        I
    3 CP         I
    4 NSPAN      I
    5 LSPAN      I
    6 NTHRY      I
    7 NTHICK     I
    8 UNDEF   none
    9 X1        RS
    10 Y1       RS
    11 Z1       RS
    12 X12      RS
    13 X4       RS
    14 Y4       RS
    15 Z4       RS
    16 X43      RS

    """
    structi = Struct(endian + b'8i 8f')

    key = (5001, 50, 175)
    nfields = 16
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for caero_id in caero_ids:
        caero = model.caeros[caero_id] # type: CAERO5
        x1, y1, z1 = caero.p1
        x4, y4, z4 = caero.p4
        #print(caero.get_stats())
        # 4 NSPAN      I
        # 5 LSPAN      I
        # 6 NTHRY      I
        # 7 NTHICK     I
        # 8 UNDEF   none
        # 9 X1        RS

        data = [
            caero.eid, caero.pid, caero.cp,
            caero.nspan, caero.lspan, caero.ntheory,
            caero.nthick, 0,
            x1, y1, z1, caero.x12,
            x4, y4, z4, caero.x43,]

        assert None not in data, data
        op2_ascii.write(f'  CAERO5 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_paero1(model: Union[BDF, OP2Geom], name: str,
                 paero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
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
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for paero_id in paero_ids:
        paero = model.paeros[paero_id] # type: PAERO1
        n0s = 7 - len(paero.caero_body_ids)
        data = [paero.pid, ] + paero.caero_body_ids
        if n0s:
            data += [0] * n0s
        assert None not in data, data
        op2_ascii.write(f'  PAERO1 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_paero2(model: Union[BDF, OP2Geom], name: str,
                 paero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 PID        I
    2 ORIENT CHAR4
    3 UNDEF   none (orient carryover)
    4 WIDTH     RS
    5 AR        RS
    6 LRSB       I
    7 LRIB       I
    8 LTH1       I
    9 LTH2       I
    10 THI1      I
    11 THN1      I
    12 THI2      I
    13 THN2      I
    14 THI3      I
    15 THN3      I
    PAERO2  100001

    data = (100001, 100001, 0, 10, 0, 0, 24, 1,
            99.3, 21.45, -11.65, 42.86, 101.8387, 122.62, -2.69, 32.71)

    """
    key = (4601, 46, 170)
    nfields = 15
    structi = Struct(endian + b'i 4s 3f 10i')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for paero_id in paero_ids:
        paero = model.paeros[paero_id] # type: PAERO2
        orient = b'%-4s' % paero.orient.encode('ascii')
        lrsb = 0 if paero.lrsb is None else paero.lrsb
        lrib = 0 if paero.lrib is None else paero.lrib
        lth1 = 0 if paero.lth[0] is None else paero.lth[0]
        lth2 = 0 if paero.lth[1] is None else paero.lth[1]

        data = [
            paero.pid, orient, 0, paero.width, paero.AR,
            lrsb, lrib, lth1, lth2,
        ]
        # we need to copy this so we don't modify the original
        thi = deepcopy(paero.thi)
        nthi = len(thi)
        if nthi < 3:
            thi.extend([0]*(3-nthi))

        thn = deepcopy(paero.thn)
        nthn = len(paero.thn)
        if nthn < 3:
            thn.extend([0]*(3-nthn))

        #assert len(paero.thi) == 3, paero.thi
        #assert len(paero.ltn) == 3, paero.ltn
        for thii, thni in zip(thi, thn):
            assert isinstance(thii, integer_types), f'i={i} thi={thii}'
            assert isinstance(thni, integer_types), f'i={i} thn={thni}'
            data.extend([thii, thni])
        assert None not in data, data
        op2_ascii.write(f'  PAERO2 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
    return nbytes

def write_paero5(model: Union[BDF, OP2Geom], name: str,
                 paero_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2
    Word Name Type Description
    1 PID    I
    2 NALPHA I
    3 LALPHA I
    4 NXIS   I
    5 LXIS   I
    6 NTAUS  I
    7 LTAUS  I
    8 CAOCI RS
    Word 8 repeats until End of Record
    """
    key = (5101, 51, 176)

    ncaoci = 0
    for paero_id in paero_ids:
        paero = model.paeros[paero_id] # type: PAERO5
        ncaoci += len(paero.caoci)

    nvalues = 8 * ncards + ncaoci
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for paero_id in paero_ids:
        paero = model.paeros[paero_id] # type: PAERO5
        ncaocii = len(paero.caoci)

        fmt = endian + b'7i' + b'f' * ncaocii + b'i'
        structi = Struct(fmt)
        #lrsb = 0 if paero.lrsb is None else paero.lrsb
        #lrib = 0 if paero.lrib is None else paero.lrib
        #lth1 = 0 if paero.lth[0] is None else paero.lth[0]
        #lth2 = 0 if paero.lth[1] is None else paero.lth[1]

        data = [
            # 2 NALPHA I
            # 3 LALPHA I
            # 4 NXIS   I
            # 5 LXIS   I
            # 6 NTAUS  I
            # 7 LTAUS  I
            paero.pid,
            paero.nalpha, paero.lalpha,
            paero.nxis, paero.lxis,
            paero.ntaus, paero.ltaus] + list(paero.caoci) + [-1]
        #print(data)

        assert None not in data, data
        op2_ascii.write(f'  PAERO5 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
        all_data.extend(data)
    assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'
    return nbytes

def write_spline1(model: Union[BDF, OP2Geom], name: str,
                  spline_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Word Name Type Description
    1 EID   I
    2 CAERO I
    3 BOX1  I
    4 BOX2  I
    5 SETG  I
    6 DZ    RS
    7 METHOD(2) CHAR4 Method: IPS|TPS|FPS
    9 USAGE(2) CHAR4 Usage flag: FORCE|DISP|BOTH
    11 NELEM I Number of elements for FPS on x-axis
    12 MELEM I Number of elements for FPS on y-axis

    """
    key = (3302, 33, 266)
    nfields = 12
    structi = Struct(endian + b'5if 8s 8s 2i')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for spline_id in spline_ids:
        spline = model.splines[spline_id] # type: SPLINE1

        method = b'%-8s' % spline.method.encode('ascii')
        usage = b'%-8s' % spline.usage.encode('ascii')
        data = [spline.eid, spline.CAero(), spline.box1, spline.box2,
                spline.Set(), spline.dz, method, usage, spline.nelements, spline.melements]
        #print(spline.get_stats())
        assert None not in data, data
        op2_ascii.write(f'  SPLINE1 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_spline2(model: Union[BDF, OP2Geom], name: str,
                  spline_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Writes the SPLINE2 card

    Word Name Type Description
    1 EID   I
    2 CAERO I
    3 ID1   I
    4 ID2   I
    5 SETG  I
    6 DZ    RS
    7 DTOR  RS
    8 CID   I
    9 DTHX  RS
    10 DTHY RS
    11 USAGE(2) CHAR4 Usage flag: FORCE|DISP|BOTH

    """
    key = (3402, 34, 267)
    nfields = 12
    structi = Struct(endian + b'5i 2f i 2f 8s')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for spline_id in spline_ids:
        spline = model.splines[spline_id] # type: SPLINE1

        usage_bytes = b'%-8s' % spline.usage.encode('ascii')
        dthx = 0.0 if spline.dthx is None else spline.dthx
        dthy = 0.0 if spline.dthx is None else spline.dthx
        #print(spline.get_stats())
        data = [spline.eid, spline.caero, spline.box1, spline.box2,
                spline.setg, spline.dz, spline.dtor, spline.cid,
                dthx, dthy, usage_bytes]
        assert None not in data, data
        op2_ascii.write(f'  SPLINE2 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_spline3(model: Union[BDF, OP2Geom], name: str,
                  spline_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:  # pragma: no cover
    """
    Writes the SPLINE3 card

    Word Name Type Description
    1 EID          I Element identification number
    2 CAERO        I Identification number of the macro-element on which
                     the element to be interpolated lies
    3 BOXID        I Identification number of the aerodynamic element
    4 COMP         I The component of motion to be interpolated
    5 USAGE(2) CHAR4 Spline usage flag to determine whether this spline applies
                     to the force transformation, displacement transformation, or both:
                     FORCE, DISP, or BOTH
    7 G            I Identification number of the independent grid point
    8 C            I Component number in the displacement coordinate system
    9 A           RS Coefficient of the constraint relationship
    Words 7 thru 9 repeat until -1 occurs

    How do you interpret the -1. Is it -1 or (-1, -1, -1)?
    """
    key = (4901, 49, 173)
    nvalues = 0
    msg = b''
    for spline_id in spline_ids:
        spline = model.splines[spline_id]
        ngrid = len(spline.nodes)
        nvalues += 6 + ngrid * 3

    for spline_id in spline_ids:
        spline = model.splines[spline_id]
        usage = spline.usage.encode('ascii')
        datai = [spline.eid, spline.caero, spline.box_id, spline.components, usage]
        fmt = '4i 8s'
        for nid, comp, coeff in zip(spline.nodes, spline.displacement_components, spline.coeffs):
            fmt += '2if'
            datai.extend([nid, comp, coeff])
        #values.extend(datai)
        assert None not in datai, datai
        msgi = Struct(fmt).pack(*datai)
        msg += msgi
        op2_ascii.write(f'  SPLINE3 data={datai}\n')
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    op2_file.write(msg)

    # structi = Struct(endian + b'5i 2f i 2f 8s')
    # nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)
    #
    # for spline_id in spline_ids:
    #     spline = model.splines[spline_id]
    #
    #     usage_bytes = b'%-8s' % spline.usage.encode('ascii')
    #     dthx = 0.0 if spline.dthx is None else spline.dthx
    #     dthy = 0.0 if spline.dthx is None else spline.dthx
    #     #print(spline.get_stats())
    #     data = [spline.eid, spline.caero, spline.box1, spline.box2,
    #             spline.setg, spline.dz, spline.dtor, spline.cid,
    #             dthx, dthy, usage_bytes]
    #     assert None not in data, data
    #     op2_ascii.write(f'  SPLINE3 data={data}\n')
    #     op2_file.write(structi.pack(*data))
    return nbytes

def write_aesurf(model: Union[BDF, OP2Geom], name: str,
                 aesurf_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 ID           I Identification of an aerodynamic trim variable degree of freedom >0, no default
    2 LABEL(2) CHAR4 Control Surface (CS) name, no default
    4 CID1         I IDentification of a rectangular Coordinate system with y-axis that defines the hinge line of the CS component
    5 ALID1        I IDentification of an AELIST bulk data entry that identifies all aerodynamic elements that make up the CS comp
    6 CID2         I IDentification of a rectangular Coordinate system with y-axis that defines the hinge line of the CS component
    7 ALID2        I IDentification of an AELIST bulk data entry that identifies all aerodynamic elements that make up the CS comp
    8 EFF         RS Control surface EFFectiveness, default=1.0
    9 LDW          I =0 create a linear down wash, >0 no linear downwash
    10 CREFC      RS REFerence Chord Length for the CS > 0.0, default = 1.0
    11 CREFS      RS REFerence area for the CS > 0.0, default = 1.0
    12 PLLIM      RS Lower deflection   Limit for the control surface in radians, default=no limit
    13 PULIM      RS Upper deflection   Limit for the control surface in radians, default=no limit
    14 HMLLIM     RS Lower Hinge Moment Limit for the control surface in force-length units, default=no limit
    15 HMULIM     RS Upper Hinge Moment Limit for the control surface in force-length units, default=no limit
    16 TQLLIM     RS Lower deflection   Limit for the control surface as fct(q), >0, default=no limit
    17 TQULIM     RS Upper deflection   Limit for the control surface as fct(q), >0, default=no limit

    """
    key = (2002, 20, 338)
    nfields = 17
    structi = Struct(endian + b'i 8s 4i fi 2f')

    structf = Struct(endian + b'f')
    structs = Struct(endian + b'4s')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for aesurf_id in aesurf_ids:
        aesurf = model.aesurf[aesurf_id] # type: AESURF

        label_bytes = b'%-8s' % aesurf.label.encode('ascii')

        if aesurf.ldw == 'LDW':
            ldw_int = 0
        elif aesurf.ldw == 'NOLDW':
            ldw_int = 1
        else:  # pragma: no cover
            raise NotImplementedError(aesurf.ldw)
        cid2 = 0 if aesurf.cid2 is None else aesurf.cid2
        aelist_id2 = 0 if aesurf.aelist_id2 is None else aesurf.aelist_id2

        #print(aesurf.get_stats())
        data = [
            aesurf.aesid, label_bytes,
            aesurf.cid1, aesurf.aelist_id1,
            cid2, aelist_id2,
            aesurf.eff, ldw_int,
            aesurf.crefc, aesurf.crefs,
        ]
        assert None not in data, data
        op2_file.write(structi.pack(*data))

        float_strs = [
            aesurf.pllim, aesurf.pulim,
            aesurf.hmllim, aesurf.hmulim,
            aesurf.tqllim, aesurf.tqulim]
        for val in float_strs:
            if val is None:
                op2_file.write(structs.pack(b'    '))
            else:
                op2_file.write(structf.pack(val))
        op2_ascii.write(f'  AESURF data={data}; float_strs={float_strs}\n')
    return nbytes

def write_aesurfs(model: Union[BDF, OP2Geom], name: str,
                  aesurfs_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Word Name Type Description
    1 ID       I     Identification of an aerodynamic trim variable degree
                     of freedom >0, no default
    2 LABEL(2) CHAR4 Control Surface (CS) name, no default
    4 LIST1    I     Identification of a SET1 that contains the grids ids
                     associated with this control surface
    5 LIST2    I     Identification of a SET1 that contains the grids ids
                     associated with this control surface
    """
    key = (7701, 77, 581)
    nfields = 5
    structi = Struct(endian + b'i 8s 2i')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for aesurfs_id in aesurfs_ids:
        aesurfs = model.aesurfs[aesurfs_id] # type: AESURFS
        label_bytes = b'%-8s' % aesurfs.label.encode('ascii')

        #print(aesurfs.get_stats())
        data = [
            aesurfs.aesid, label_bytes,
            aesurfs.list1, aesurfs.list2]
        assert None not in data, data
        op2_ascii.write(f'  AESURFS data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_aestat(model: Union[BDF, OP2Geom], name: str,
                 aestat_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 ID I
    2 LABEL(2) CHAR4

    """
    key = (2102, 21, 339)
    nfields = 3
    structi = Struct(endian + b'i 8s')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for aestat_id in aestat_ids:
        aestat = model.aestats[aestat_id] # type: AESTAT
        label_bytes = b'%-8s' % aestat.label.encode('ascii')

        #print(aestat.get_stats())
        data = [aestat.aestat_id, label_bytes]
        assert None not in data, data
        op2_ascii.write(f'  AESTAT data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_flutter(model: Union[BDF, OP2Geom], name: str,
                  flutter_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    (3902, 39, 272)
    MSC 2018.2

    Word Name Type Description
    1 SID           I
    2 METHOD(2) CHAR4
    4 DENS          I
    5 MACH          I
    6 RFREQ         I
    7 IMETH(2)  CHAR4
    SFLG=0 (std)
      9 NEIGN  I  nvalue
      10 EPR  RS
      11 SFLG  I SWEEP FLAG
    SFLG=1 (sweep)
      9 FMAX  RS maximum frequency
      10 EPR  RS
      11 SFLG  I SWEEP FLAG
    End SFLG
    Words 1 through max repeat until End of Record

    NX:
      data = (30, PK, 1, 2, 3, L, 3, 0.001, -1)
    """
    key = (3902, 39, 272)
    nfields = 11
    structi = Struct(endian + b'i 8s 3i 8s ifi')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for flutter_id in flutter_ids:
        flutter = model.flutters[flutter_id]  # type: FLUTTER
        #print(flutter.get_stats())
        assert isinstance(flutter.method, str), f'flutter.method={flutter.method}; type={type(flutter.method)}'
        assert isinstance(flutter.imethod, str), f'flutter.imethod={flutter.imethod}; type={type(flutter.imethod)}'

        #if isinstance(flutter.method, str):
        method = b'%-8s' % flutter.method.encode('latin1')
        imethod = b'%-8s' % flutter.imethod.encode('latin1')
        #if 0:
            #try:
                #method = b'-%8s' % flutter.method.decode('latin1')
                #imethod = b'-%8s' % flutter.imethod.decode('latin1')
            #except AttributeError:
                #print(flutter.get_stats())
                #raise
        nvalue = flutter.nvalue
        if nvalue is None:
            nvalue = 0
        #assert isinstance(flutter.nvalue, int), f'flutter.nvalue={flutter.nvalue}; type={type(flutter.nvalue)}'

        data = [flutter.sid, method,
                flutter.density, flutter.mach, flutter.reduced_freq_velocity,
                imethod, nvalue, flutter.epsilon, -1]
        assert len(data) == 9, data
        assert None not in data, data
        op2_ascii.write(f'  FLUTTER data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def _makero_temp(data, i: int, nloops: int):
    if i == (nloops - 1):
        data_temp = data[i*8:]
    else:
        #machs_temp = [-1] * 8
        data_temp = data[i*8:i*8+8]
    return data_temp

def write_mkaero1(model: Union[BDF, OP2Geom], name: str,
                  mkaero1s: list[MKAERO1], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """writes the MKAERO1

    data = (1.3, -1, -1, -1, -1, -1, -1, -1,
           0.03, 0.04, 0.05, -1, -1, -1, -1, -1)
    """
    key = (3802, 38, 271)

    makero1s_temp = []
    makero1s_final = []
    for mkaero in mkaero1s:
        nmachs = len(mkaero.machs)
        nkfreqs = len(mkaero.reduced_freqs)
        assert nmachs > 0, mkaero
        assert nkfreqs > 0, mkaero

        if nmachs <= 8 and nkfreqs <= 8:
            # no splitting required
            makero1s_final.append((mkaero.machs, mkaero.reduced_freqs))
        elif nmachs <= 8 or nkfreqs <= 8:
            # one of machs or kfreqs < 8
            makero1s_temp.append((mkaero.machs, mkaero.reduced_freqs))
        else:
            # both machs and kfreqs > 8
            nloops_mach = int(np.ceil(nmachs/8))
            for i in range(nloops_mach):
                machs_temp = _makero_temp(mkaero.machs, i, nloops_mach)
                assert len(machs_temp) > 0, (i, nloops_mach, machs_temp)
                makero1s_temp.append((machs_temp, mkaero.reduced_freqs))

    for (machs, reduced_freqs) in makero1s_temp:
        nmachs = len(machs)
        nkfreqs = len(reduced_freqs)
        assert nmachs > 0, nmachs
        assert nkfreqs > 0, nkfreqs

        if nmachs <= 8 and nkfreqs <= 8:  # pragma: no cover
            raise RuntimeError(f'this should never happen...nmachs={nmachs} knfreqs={nkfreqs}')
        if nmachs <= 8:
            # nkfreqs > 8
            nloops = int(np.ceil(nkfreqs/8))
            for i in range(nloops):
                reduced_freqs_temp = _makero_temp(reduced_freqs, i, nloops)
                makero1s_final.append((machs, reduced_freqs_temp))
        elif nkfreqs <= 8:
            # nmachs > 8
            nloops = int(np.ceil(nmachs/8))
            for i in range(nloops):
                machs_temp = _makero_temp(machs, i, nloops)
                assert len(machs_temp) > 0, (i, nloops_mach, machs_temp)
                makero1s_final.append((machs_temp, reduced_freqs))
        else:  # pragma: no cover
            raise RuntimeError(f'this should never happen...nmachs={nmachs} knfreqs={nkfreqs}')
            #raise RuntimeError((nmachs, nkfreqs))

    ncards = len(makero1s_final)
    nfields = 16
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for machs, reduced_freqs in makero1s_final:
        data = []
        nmachs = len(machs)
        nkfreqs = len(reduced_freqs)
        assert nmachs > 0, machs
        assert nkfreqs > 0, reduced_freqs

        nint_mach = 8 - nmachs
        nint_kfreq = 8 - nkfreqs
        fmt1 = b'%if' % nmachs + b'i' * nint_mach
        fmt2 = b'%if' % nkfreqs + b'i' * nint_kfreq
        spack = Struct(endian + fmt1 + fmt2)

        data.extend(machs.tolist())
        assert nint_mach < 8, nint_mach
        if nint_mach:
            data.extend([-1]*nint_mach)
        data.extend(reduced_freqs.tolist())
        if nint_kfreq:
            data.extend([-1]*nint_kfreq)
        op2_ascii.write(f'  mkaero1 data={data}\n')
        op2_file.write(spack.pack(*data))
    return nbytes

def write_aero(model: Union[BDF, OP2Geom], name: str,
               aero: list[AERO], ncards: int,
               op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
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
    spack = Struct(endian + b'i 3f 2i')

    key = (3202, 32, 265)
    nfields = 6
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

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
    op2_file.write(spack.pack(*data))
    return nbytes

def write_aeros(model: Union[BDF, OP2Geom], name: str,
                aeros: list[AEROS], ncards: int,
                op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
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
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    data = [aeroi.acsid, aeroi.rcsid,
            aeroi.cref, aeroi.bref, aeroi.sref,
            aeroi.sym_xz, aeroi.sym_xy]
    #print(data)
    assert None not in data, data
    op2_ascii.write(f'  AEROS data={data}\n')
    op2_file.write(spack.pack(*data))
    return nbytes

def write_flfact(model: Union[BDF, OP2Geom], name: str,
                 flfact_ids, ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    (4102, 41, 274)
    NX 2019.2

    data = (1, 0.206, -1,
            2, 1.3, -1,
            3, 14400.0, 15600.0, 16800.0, 18000.0, 19200.0, 20400.0, -1)

    """
    key = (4102, 41, 274)
    #nfields = 3

    fmt = ''
    data = []
    for flfact_id in flfact_ids:
        flfact = model.flfacts[flfact_id]
        #print(flfact.get_stats())
        nfactors = len(flfact.factors)
        fmt += 'i %if i' % nfactors
        if isinstance(flfact.factors, np.ndarray):
            factors = flfact.factors.tolist()
        else:
            factors = flfact.factors
        datai = [flfact.sid] + factors + [-1]
        data.extend(datai)
        #flutter = model.loads[flutter_id]  # type: FLUTTER
        assert None not in datai, datai
        op2_ascii.write(f'  FLFACT data={data}\n')

    nvalues = len(data)
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)
    structi = Struct(fmt)
    op2_file.write(structi.pack(*data))
    return nbytes

def write_set1(model: Union[BDF, OP2Geom], name: str,
               set_ids: list[int], ncards: int,
               op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    SET1: (3502, 35, 268)
    MSC 2018.2
    Word Name Type Description
    1 SID I
    2 G1  I Grid ID or -2 when SKIN is specified
    Word 2 repeats until End of Record
    """
    key = (3502, 35, 268)

    ngrids = 0
    for set_id in set_ids:
        seti = model.sets[set_id] # type: PAERO5
        ngrids += len(seti.ids)
        if seti.is_skin:
            ngrids += 1

    # 2* = sid and the -1 flag
    nvalues = 2 * ncards + ngrids
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for set_id in set_ids:
        seti = model.sets[set_id] # type: PAERO5
        ngridsi = len(seti.ids)

        # +2 = sid and the -1 flag
        dn = 2
        if seti.is_skin:
            dn += 1
        fmt = endian + b'%di' % (ngridsi + dn)
        #print(fmt)
        structi = Struct(fmt)

        data = [
            # 2 NALPHA I
            # 3 LALPHA I
            # 4 NXIS   I
            # 5 LXIS   I
            # 6 NTAUS  I
            # 7 LTAUS  I
            seti.sid,
        ]
        if seti.is_skin:
            data.append(-2)
        data += seti.ids
        data.append(-1)

        assert None not in data, data
        op2_ascii.write(f'  SET1 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
        all_data.extend(data)
    assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'
    return nbytes

def write_set3(model: Union[BDF, OP2Geom], name: str,
               set_ids: list[int], ncards: int,
               op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2
    Word Name Type Description
    1 SID I
    2 DES I Set description:
        1=ELEM
        2=GRID
        3=POINT
        4=PROP
        5=RBEin
        6=RBEex
    3 ID1 I IDs of Grids, Elements, Points or Properties.
    4 "ID1 THRU ID2" format will be EXPANDED into explicit IDs.
    Words 3 through 4 repeat until End of Record

    data = (1, 1, 190, ..., 238, -1,
            2, 1, 71, ..., 189, -1,
            4, 1, 309, ..., ..., 378, -1)
    """
    key = (8001, 80, 511)
    raise NotImplementedError('SET3')

    ngrids = 0
    for set_id in set_ids:
        seti = model.sets[set_id] # type: PAERO5
        ngrids += len(seti.ids)
        if seti.is_skin:
            ngrids += 1

    # 2* = sid and the -1 flag
    nvalues = 2 * ncards + ngrids
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for set_id in set_ids:
        seti = model.sets[set_id] # type: PAERO5
        ngridsi = len(seti.ids)

        # +2 = sid and the -1 flag
        fmt = endian + b'%di' % (ngridsi + 2)
        #print(fmt)
        structi = Struct(fmt)

        data = [
            # 2 NALPHA I
            # 3 LALPHA I
            # 4 NXIS   I
            # 5 LXIS   I
            # 6 NTAUS  I
            # 7 LTAUS  I
            seti.sid,
        ]
        if seti.is_skin:
            data.append(-2)
        data += seti.ids
        data.append(-1)
        #print(data)

        assert None not in data, data
        op2_ascii.write(f'  SET3 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
        all_data.extend(data)
    assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'
    return nbytes

def write_aelink(model: Union[BDF, OP2Geom], name: str,
                 aelink_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 ID I
    2 LABLD(2) CHAR4
    4 LABLI(2) CHAR4
    6 C1 RS
    Words 4 through 6 repeat until (-1,-1,-1) occurs

    """
    aelink_ids = np.unique(aelink_ids)
    # 1, b'AILR_L  '
    # b'AILR_R  ', -1.0
    # (-1, -1, -1)
    #
    # 1, b'ELEV_L  '
    # b'ELEV_R  ', 1.0
    # (-1, -1, -1)
    #
    # 2, b'ELEV_L  '
    # b'ELEV_R  ', 1.0
    # (-1, -1, -1)
    #
    # 2, b'AILR_L  '
    # b'AILR_R  ', -1.0
    # (-1, -1, -1)
    #
    # 3, b'AILR_L  '
    # b'AILR_R  ', -1.0
    # (-1, -1, -1)
    #
    # 3, b'ELEV_L  '
    # b'ELEV_R  ', 1.0
    # (-1, -1, -1)
    #op2 = self.op2
    #struct1 = Struct(op2._endian + b'i8s')
    #struct2 =Struct(op2._endian + b'8sf')
    #struct_end = Struct(op2._endian + b'3i')
    key = (2602, 26, 386)

    ncoeffs = 0
    for aelink_id in aelink_ids:
        aelinks = model.aelinks[aelink_id] # type: PAERO5
        for aelink in aelinks:
            #print(aelink.get_stats())
            ncoeffs += len(aelink.linking_coefficients)

    # 3 (main value)
    # 3 (-1)s
    # -------
    # 6
    #
    # 3 values per coeff
    nvalues = 6 * ncards + ncoeffs * 3
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for aelink_id in aelink_ids:
        for aelink in aelinks:
            aelinks = model.aelinks[aelink_id] # type: AELINK
            ncoeffsi = len(aelink.linking_coefficients)

            fmti = b'8sf' * ncoeffsi
            fmt = endian + b'i8s ' + fmti + b'3i'
            #print(fmt)
            structi = Struct(fmt)
            assert isinstance(aelink_id, integer_types) and aelink_id >= 0, aelink_id
            label = b'%-8s' % aelink.label.encode('latin1')
            data = [
                aelink_id, label,
            ]
            for ind_label, coeff in zip(aelink.independent_labels, aelink.linking_coefficients):
                ind_label_bytes = b'%-8s' % ind_label.encode('latin1')
                data.extend([ind_label_bytes, coeff])
            data.extend([-1, -1, -1])
            op2_ascii.write(f'  AELINK data={data}\n')
            assert None not in data, data
            op2_file.write(structi.pack(*data))
            all_data.extend(data)

    #assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'
    return nbytes

def write_monpnt1(model: Union[BDF, OP2Geom], name: str,
                  monpnt_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2
    NX-92

    Word Name Type Description
    1 NAME(2)   CHAR4
    3 LABEL(14) CHAR4
    17 AXES         I
    18 COMP(2)  CHAR4
    20 CP           I
    21 X           RS
    22 Y           RS
    23 Z           RS

    """
    key = (7601, 76, 577)
    if nastran_format == 'msc':
        nfields = 24
        structi = Struct(endian + b'8s 56s i 8s i 3f i')  # msc
    else:
        nfields = 23
        structi = Struct(endian + b'8s 56s i 8s i 3f')  # nx

    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)
    #op2 = self.op2
    #ntotal = 4 * 24 # 4 * 24
    #ntotal = 92 # 4 * 23
    #ndatai = len(data) - n
    #ncards = ndatai // ntotal
    #assert ndatai % ntotal == 0
    for monpnt_id in monpnt_ids:
        monitor = model.monitor_points[monpnt_id]
        name_bytes = ('%-8s' % monitor.name).encode('latin1')
        label_bytes = ('%-56s' % monitor.label).encode('latin1')
        aecomp_name_bytes = ('%-8s' % monitor.comp).encode('latin1')
        #name_bytes, label_bytes, axes, comp_bytes, cp, x, y, z, cd = out
        axes = int(monitor.axes)
        p1 = monitor.xyz
        data = [name_bytes, label_bytes, axes, aecomp_name_bytes, monitor.cp] + list(p1)
        assert len(data) == 8, f'data={data} ndata={len(data)}'
        if nastran_format == 'msc':
            data.append(monitor.cd)
        assert None not in data, data
        op2_ascii.write(f'  MONPNT1 data={data}\n')
        #print(f'  MONPNT1 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
        #all_data.extend(data)
    return nbytes

def write_monpnt2(model: Union[BDF, OP2Geom], name: str,
                  monpnt_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Record 59 - MONPNT2(8204,82,621)

    MSC
    ===
    Word Name Type Description
    1 NAME(2)   CHAR4
    3 LABEL(14) CHAR4
    17 TABLE(2) CHAR4
    19 ELTYP(2) CHAR4
    21 ITEM(2)  CHAR4
    23 EID      I

    NX
    ==
    Word Name Type Description
    1 NAME(2)   CHAR4
    3 LABEL(14) CHAR4
    17 TABLE(2) CHAR4
    19 ELTYP(2) CHAR4
    21 ITEM(2)  CHAR4
    23 EID      I
    Words 17 thru 23 repeat until -1 occurs
    Words 1 thru 23 repeat until (-2,-2) occurs

    """
    key = (8204, 82, 621)
    nvalues = 0
    #nfields = 23
    #nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)
    #structi = Struct(endian + b'8s 56s 8s8s8s i')  # msc
    fmt = endian
    for monpnt_id in monpnt_ids:
        monitor = model.monitor_points[monpnt_id]
        #print(monitor.get_stats())
        name_bytes = ('%-8s' % monitor.name).encode('latin1')
        label_bytes = ('%-56s' % monitor.label).encode('latin1')
        data = [name_bytes, label_bytes]

        fmt = endian
        if nastran_format == 'msc':
            assert len(monitor.tables) == 1, monitor.tables
            assert len(monitor.element_types) == 1, monitor.element_types
            assert len(monitor.nddl_items) == 1, monitor.nddl_items
            assert len(monitor.eids) == 1, monitor.eids
            for table, element_type, nddl_item, eids in zip(monitor.tables, monitor.element_types, monitor.nddl_items, monitor.eids):
                assert len(eids) == 1, eids
                table_bytes = ('%-8s' % table).encode('latin1')
                eltype_bytes = ('%-8s' % element_type).encode('latin1')
                item_bytes = ('%-8s' % nddl_item).encode('latin1')
                data.extend([table_bytes, eltype_bytes, item_bytes]+eids)
                fmt += b'8s 56s 8s8s8s i'
                nvalues += 23
                Struct(fmt).pack(*data)
            assert len(data) == 6, f'data={data} ndata={len(data)}'
        else:
            # TODO: should have -1 and -2 flags, but I need to update nbytes for the header...
            #assert len(monitor.tables) == 1, monitor.tables
            #assert len(monitor.element_types) == 1, monitor.element_types
            #assert len(monitor.nddl_items) == 1, monitor.nddl_items
            #assert len(monitor.eids) == 1, monitor.eids
            fmt = b'8s 56s'
            nvalues += 16
            for table, element_type, nddl_item, eids in zip(monitor.tables, monitor.element_types, monitor.nddl_items, monitor.eids):
                assert len(eids) == 1, eids
                table_bytes = ('%-8s' % table).encode('latin1')
                eltype_bytes = ('%-8s' % element_type).encode('latin1')
                item_bytes = ('%-8s' % nddl_item).encode('latin1')
                data.extend([table_bytes, eltype_bytes, item_bytes]+eids+[-1])
                neids = len(eids)
                nvalues += (3 * 2) + neids + 1

                fmt += b'8s8s8s i ' + neids * b'i'
                Struct(fmt).pack(*data)
            fmt += b'2i'
            data.extend([-2, -2])
            nvalues += 2
            Struct(fmt).pack(*data)
        #name_bytes, label_bytes, table_bytes, eltype_bytes, item_bytes, eid = out
        #name = reshape_bytes_block_size(name_bytes, self.size)
        #label = reshape_bytes_block_size(label_bytes, self.size)
        #table = reshape_bytes_block_size(table_bytes, self.size)
        #Type = reshape_bytes_block_size(eltype_bytes, self.size)
        #nddl_item = reshape_bytes_block_size(item_bytes, self.size)

        assert nvalues > 0, nvalues
        nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

        assert None not in data, data
        op2_ascii.write(f'  MONPNT2 data={data}\n')
        #print(f'  MONPNT2 data={data}\n')
        structi = Struct(fmt)
        op2_file.write(structi.pack(*data))
    return nbytes

def write_monpnt3(model: Union[BDF, OP2Geom], name: str,
                  monpnt_ids: list[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    Record 60 - MONPNT3(8304,83,622)
    Word Name    Type   Description
    1  NAME(2)   CHAR4
    3  LABEL(14) CHAR4
    17 AXES      I      Axes to compute
    18 GRIDSET   I      GPF Grid Set
    19 ELEMSET   I      GPF Elem Set
    20 CID       I      Coord system x,y,z input in
    21 X         RS
    22 Y         RS
    23 Z         RS
    24 XFLAG     I      Exclude forces from class
    25 CD        I

    """
    key = (8304, 83, 622)
    nfields = 25
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)
    structi = Struct(endian + b'8s 56s 4i 3f 2i')  # msc
    for monpnt_id in monpnt_ids:
        monitor = model.monitor_points[monpnt_id]
        #print(monitor.get_stats())
        name_bytes = ('%-8s' % monitor.name).encode('latin1')
        label_bytes = ('%-56s' % monitor.label).encode('latin1')

        xflag = 4
        axes = int(monitor.axes)
        data = [name_bytes, label_bytes, axes, monitor.grid_set, monitor.elem_set, monitor.cp] + list(monitor.xyz) + [xflag, monitor.cd]
        assert len(data) == 11, f'data={data} ndata={len(data)}'
        assert None not in data, data
        op2_ascii.write(f'  MONPNT3 data={data}\n')
        #print(f'  MONPNT3 data={data}\n')
        #print('npaero2', len(data), data)
        op2_file.write(structi.pack(*data))
        #all_data.extend(data)
    return nbytes


def write_deform(model: Union[BDF, OP2Geom], name: str,
                 loads: list[DEFORM], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    (104, 1, 81)
    NX 2019.2

    Word Name Type Description
    1 SID I Deformation set identification number
    2 EID I Element number
    3 D RS Deformation

    """
    key = (104, 1, 81)
    nfields = 3
    structi = Struct(endian + b'iif')
    nbytes = write_header(name, nfields, ncards, key, op2_file, op2_ascii)

    for load in loads:
        data = [load.sid, load.eid, load.deformation]
        #flutter = model.loads[flutter_id]  # type: FLUTTER
        #print(flutter.get_stats())
        assert None not in data, data
        op2_ascii.write(f'  DEFORM data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def write_aefact(model: Union[BDF, OP2Geom], name: str,
                 aefact_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 SID I
    2 D RS
    Word 2 repeats until End of Record

    (1, 0.0, 0.1, 0.2, 1.0, -1,
     2, 0.0, 0.1, 0.2, 0.5, 1.0, -1,
    )

    """
    key = (4002, 40, 273)

    nvalues = 0
    for aefact_id in aefact_ids:
        aefacti = model.aefacts[aefact_id]
        nvalues += 2 + len(aefacti.factors)

    # 2* = sid and the -1 flag
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for aefact_id in aefact_ids:
        aefacti = model.aefacts[aefact_id]
        nfactors = len(aefacti.factors)

        fmt = endian + b'i %df i' % nfactors
        #print(fmt)
        structi = Struct(fmt)

        data = [aefacti.sid,] + aefacti.factors + [-1]

        assert None not in data, data
        op2_ascii.write(f'  AEFACT data={data}\n')
        op2_file.write(structi.pack(*data))
        all_data.extend(data)
    assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'
    return nbytes

def write_aelist(model: Union[BDF, OP2Geom], name: str,
                 aelist_ids: list[int], ncards: int,
                 op2_file, op2_ascii, endian: bytes, nastran_format: str='nx') -> int:
    """
    MSC 2018.2

    Word Name Type Description
    1 SID I
    2 E I
    Word 2 repeats until End of Record

    """
    key = (2302, 23, 341)

    nvalues = 0
    for aelist_id in aelist_ids:
        aelist = model.aelists[aelist_id]
        nvalues += 2 + len(aelist.elements)

    # 2* = sid and the -1 flag
    nbytes = write_header_nvalues(name, nvalues, key, op2_file, op2_ascii)

    all_data = []
    for aelist_id in aelist_ids:
        aelist = model.aelists[aelist_id]
        nelements = len(aelist.elements)

        fmt = endian + b'i %di i' % nelements
        structi = Struct(fmt)

        data = [aelist.sid,] + aelist.elements + [-1]
        assert None not in data, data
        op2_ascii.write(f'  AELIST data={data}\n')
        op2_file.write(structi.pack(*data))
        all_data.extend(data)
    assert len(all_data) == nvalues, f'ndata={len(all_data)}; nvalues={nvalues}'

    return nbytes

EDT_MAP = {
    'AEFACT': write_aefact,
    'AELIST': write_aelist,

    'MONPNT1': write_monpnt1,
    'MONPNT2': write_monpnt2,
    'MONPNT3': write_monpnt3,
    'SET1': write_set1,
    'DIVERG': write_diverg,
    'CAERO1': write_caero1,
    'CAERO2': write_caero2,
    'CAERO3': write_caero3,
    'CAERO4': write_caero4,
    'CAERO5': write_caero5,
    'PAERO1': write_paero1,
    'PAERO2': write_paero2,
    'PAERO5': write_paero5,
    'MKAERO1': write_mkaero1,
    'AERO': write_aero,
    'AEROS': write_aeros,
    'TRIM': write_trim,
    'FLUTTER': write_flutter,
    'DEFORM': write_deform,
    'FLFACT': write_flfact,
    'SPLINE1': write_spline1,
    'SPLINE2': write_spline2,
    #'SPLINE3': write_spline3,
    'AESURF': write_aesurf,
    'AESURFS': write_aesurfs,
    'AESTAT': write_aestat,
    'AELINK': write_aelink,
    #'GUST': write_gust,  # part of DIT
}
