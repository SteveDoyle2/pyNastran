"""writes the MPT/MPTS table"""
from __future__ import annotations
from collections import defaultdict
from struct import pack, Struct
from typing import List, Union, TYPE_CHECKING

import numpy as np

from .geom1_writer import write_geom_header, close_geom_table
from .geom4_writer import write_header, write_header_nvalues
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.cards.aero.static_loads import AEROS # , AESTAT, CSSCHD, DIVERG, TRIM, TRIM2
    from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, FLUTTER # , FLFACT, MKAERO2
    from pyNastran.op2.op2_geom import OP2Geom, BDF

def write_edt(op2_file, op2_ascii, model: Union[BDF, OP2Geom], endian: bytes=b'<') -> None:
    """writes the EDT/EDTS table"""
    if not hasattr(model, 'loads'):  # OP2
        return
    card_types = [
        'MKAERO1', # 'MKAERO2',
        'AERO', 'AEROS',
        'CAERO1', 'PAERO1',
        #'AELIST', 'AEFACT', 'AESURF', 'AESURFS'
        'TRIM', 'FLUTTER',
        'DEFORM',
        'FLFACT',
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
    for unused_load_id, loads in model.loads.items():
        for load in loads:
            if load.type in ['DEFORM', 'CLOAD']:
                out[load.type].append(load)

    for eid, caero in sorted(model.caeros.items()):
        out[caero.type].append(eid)
    for pid, paero in sorted(model.paeros.items()):
        out[paero.type].append(pid)
    for spline_id, spline in sorted(model.splines.items()):
        out[spline.type].append(spline_id)
    for aesurf_id, aesurf in sorted(model.aesurf.items()):
        out[aesurf.type].append(aesurf_id)
    for aesurfs_id, aesurfs in sorted(model.aesurfs.items()):
        out[aesurfs.type].append(aesurfs_id)
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
        try:
            func = EDT_MAP[name]
        except KeyError:  # pragma: no cover
            raise NotImplementedError(name)

        nbytes = func(model, name, ids, ncards, op2_file, op2_ascii, endian)
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

def remove_unsupported_cards(card_dict: Dict[str, Any],
                             card_types: List[str],
                             log: SimpleLogger):

    for card_type in list(card_dict):
        if card_type not in card_types:
            del card_dict[card_type]
            log.warning(f"removing {card_type} in OP2 writer because it's unsupported")

def _write_trim(model: Union[BDF, OP2Geom], name: str,
                trim_ids: List[int], ncards: int,
                op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_caero1(model: Union[BDF, OP2Geom], name: str,
                  caero_ids: List[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
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
        op2_file.write(structi.pack(*data))
    return nbytes

def _write_paero1(model: Union[BDF, OP2Geom], name: str,
                  paero_ids: List[int], ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
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
        paero = model.paeros[paero_id]
        n0s = 7 - len(paero.caero_body_ids)
        data = [paero.pid, ] + paero.caero_body_ids
        if n0s:
            data += [0] * n0s
        assert None not in data, data
        op2_ascii.write(f'  PAERO1 data={data}\n')
        op2_file.write(structi.pack(*data))
    return nbytes

def _write_flutter(model: Union[BDF, OP2Geom], name: str,
                   flutter_ids: List[int], ncards: int,
                   op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_mkaero1(model: Union[BDF, OP2Geom], name: str,
                   mkaero1s: List[MKAERO1], ncards: int,
                   op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_aero(model: Union[BDF, OP2Geom], name: str,
                aero: List[AERO], ncards: int,
                op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_aeros(model: Union[BDF, OP2Geom], name: str,
                 aeros: List[AEROS], ncards: int,
                 op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_deform(model: Union[BDF, OP2Geom], name: str,
                  loads: List[AEROS], ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
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

def _write_flfact(model: Union[BDF, OP2Geom], name: str,
                  flfact_ids, ncards: int,
                  op2_file, op2_ascii, endian: bytes) -> int:
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

EDT_MAP = {
    'CAERO1': _write_caero1,
    'PAERO1': _write_paero1,
    'MKAERO1': _write_mkaero1,
    'AERO': _write_aero,
    'AEROS': _write_aeros,
    'TRIM': _write_trim,
    'FLUTTER': _write_flutter,
    'DEFORM': _write_deform,
    'FLFACT': _write_flfact,
}
